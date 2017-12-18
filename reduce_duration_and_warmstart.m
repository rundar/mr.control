function [d,q,skip] = reduce_duration_and_warmstart(d,q,N_distribute,m)
%     written by Armin Rund, Dec 18, 2017
%############## Copyright: ##################################
%  Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                Christoph Aigner, christoph.aigner@tugraz.at
% Implementation of the time optimal control method of paper
% "Simultaneous Multislice Refocusing via Time Optimal Control"
% by A. Rund, C. Aigner, K. Kunisch, and R. Stollberger.
%%###########################################################
% cut the point m out of the time grid. distribute the RF an 
% Gs at time m symmetrically to neighboring time instances 
% while holding the amplitude constraints as well as the Gs 
% slew rate constraint.
%%###########################################################
% d (I/O) is the struct of parameters
% q (I/O) is the optimization variable
% N_distribute (I) is the number of neighbors for distribution
%       in each direction
% m (I) is the index of the time point that is cut. 1<m<d.Nu
% skip (O) is true if the algorithm cannot fully distribute the
%       cutted values
%%###########################################################

%  if (m<=1||m>=d.Nu)
%    fprintf('error in reduce_duration_and_warmstart(): cut of initial or terminal time is forbidden! Nu=%d  m=%d\n',d.Nu,m); 
%  end;    
  Told = d.T;  
  q_proj = min(d.box,max(-d.box,q));	% projecting the control
  u=q_proj(1:d.Nu);   % RF
  slew=q_proj(d.Nu+1:end);  % Gs slew rate
  w = zeros(d.Nu,1); %  Gs amplitude
  for n=2:d.Nu-1 
    w(n) = w(n-1) + d.dt*slew(n-1);
  end

  % now cut and reduce time: (working on RF and Gs)
  [u,w,skip] = cut_point(d,u,w,m,N_distribute);

  % finish the cut:
  u(m) = [];
  w(m) = [];
  % get the slew rate again:
  slew=diff(w)/d.dt; 
  slew=slew(1:end-1);
  % adapt optimization variable and parameters:
  q = [u;slew];
  d.Nu = size(u,1); d.T = d.Nu*d.dt; 
  d.box=[d.box(1)*ones(d.Nu,1);d.box(end)*ones(d.Nu-2,1)];  
end


% cut time point m and distribute that u,w to neighboring times
function [u,w,skip] = cut_point(d,u,w,m,N_distribute)
  skip = false;
  % 1. distribute the cutted RF pulse to neighboring time instances:
  box=d.box(1);
  cut = u(m); % this has to be distributed to neighbors
  k=1;
  while abs(cut) >1e-10 && k<=N_distribute  
    % first try half and half:
    half = 0.5*cut;
    if m+k<=d.Nu-1
      pre = u(m+k); u(m+k) = min(box,max(-box,u(m+k) + half)); cut = cut + pre - u(m+k);
    end
    if m-k>=2
      pre = u(m-k); u(m-k) = min(box,max(-box,u(m-k) + half)); cut = cut + pre - u(m-k);
    end
    if abs(cut) >1e-10  % at least one of them is at the box. add the rest to the other:
      if m+k<=d.Nu-1,pre = u(m+k); u(m+k) = min(box,max(-box,u(m+k) + cut)); cut = cut + pre - u(m+k); end;
      if m-k>=2,pre = u(m-k); u(m-k) = min(box,max(-box,u(m-k) + cut)); cut = cut + pre - u(m-k); end;
    end    
    k = k+1;
  end
  if abs(cut) >1e-10
    skip = true; return;
  end
  
  % 2. distribute the cutted gradient field to neighboring time instances:
  cut = w(m); % this has to be distributed to neighbors
  boxw = d.w_box - 1e-3;
  del = d.box(end)*d.dt; Lfirst=false;
  onlyLeft = false; onlyRight=false;
  % a) get the direct neighbors admissible:
  if m==d.Nu-1 % then there is no right neighbour: already admissible!
    onlyLeft=true; % just distribute to the left
  elseif m==2  
    onlyRight=true; % just distribute to the right
    % maximally increase right neighbour:
    pre = w(m+1); w(m+1) = min(min(w(m+2)+del,min(del,w(m+1)+cut)),boxw); cut = cut + pre - w(m+1);
  else % we can change the left and right neighbour:
    if abs(w(m-1)) < abs(w(m+1))-del % correct the left neighbor
      pre = w(m-1); Lfirst = true;
      w(m-1) = min(w(m+1)-sign(cut)*del,boxw);
      cut = cut + pre - w(m-1); % rest that has to be distributed
      [w,skip2,cut] = distribute(d,w,m,-1,N_distribute,del,cut); if skip2,skip=true;end;% go Left      
    elseif abs(w(m-1)) < abs(w(m+1))-0.9*del % barely admissible:
      Lfirst = true;
    elseif  abs(w(m+1)) < abs(w(m-1))-del % correct the right neighbor
      pre = w(m+1); 
      w(m+1) = min(w(m-1)-sign(cut)*del,boxw);
      cut = cut + pre - w(m+1);
      [w,skip2,cut] = distribute(d,w,m,1,N_distribute,del,cut); if skip2,skip=true;end;% go Right
    end
  end    
  % b) distribute rest to L and R:
  k = 1;
  while abs(cut)>1e-10 && k<100
    C = abs(cut)/(2*N_distribute); % conservative choice for admissibility! maybe increase if loop too slow
    if onlyLeft
      pre = w(m-1); 
      w(m-1) = min(w(m-1) +sign(cut)*C,boxw);
      cut = cut + pre - w(m-1);
      [w,skip2,cut] = distribute(d,w,m,-1,N_distribute,del,cut); if skip2,skip=true;end;% go Left           
    elseif onlyRight
      pre = w(m+1); 
      w(m+1) = min(w(m+1)+sign(cut)*C,boxw);
      cut = cut + pre - w(m+1);
      [w,skip2,cut] = distribute(d,w,m,1,N_distribute,del,cut); if skip2,skip=true;end; % go Right          
    else
      if Lfirst % first Left than Right
        pre = w(m-1); 
        w(m-1) = min(w(m-1) +sign(cut)*C,boxw);
        cut = cut + pre - w(m-1);
        [w,skip2,cut] = distribute(d,w,m,-1,N_distribute,del,cut); if skip2,skip=true;end;% go Left
        pre = w(m+1); 
        w(m+1) = min(w(m+1)+sign(cut)*C,boxw);
        cut = cut + pre - w(m+1);
        [w,skip2,cut] = distribute(d,w,m,1,N_distribute,del,cut); if skip2,skip=true;end;% go Right
      else % other way round: first Right then Left
        pre = w(m+1); 
        w(m+1) = min(w(m+1)+sign(cut)*C,boxw);
        cut = cut + pre - w(m+1);
        [w,skip2,cut] = distribute(d,w,m,1,N_distribute,del,cut); if skip2,skip=true;end; % go Right         
        pre = w(m-1); 
        w(m-1) = min(w(m-1) +sign(cut)*C,boxw);
        cut = cut + pre - w(m-1);
        [w,skip2,cut] = distribute(d,w,m,-1,N_distribute,del,cut); if skip2,skip=true;end;% go Left
      end   
    end
    k = k+1;
  end
  if abs(cut) >1e-10
    skip = true; return;
  end
  % finally, fix last slewrate:
  if abs(w(end))>del-1e-10
    cut = w(end); 
    w(end)=0;
    [w,skip2,cut] = distribute(d,w,size(w,1)+1,-1,N_distribute,del,cut); if skip2,skip=true;end;
    if abs(cut) >1e-10
      [w,skip2,cut] = distribute(d,w,size(w,1)+1,-1,N_distribute*2,del,cut);
      skip = true; return;
    end
  end
end

function [w,skip,cut] = distribute(d,w,m,s,N_distribute,del,cut)
  boxw = d.w_box - 1e-3;
  % s is the signum for the indices:   -1=Left, 1=Right
  k=2; skip = false; 
  while k<=N_distribute % consider up to N_distribute neighbors
    if m+s*k<=d.Nu-1 && m+s*k>=2    && abs(w(m+s*k)-w(m+s*(k-1)))>del  % we have to correct here
      pre = w(m+s*k); 
      w(m+s*k) = min(w(m+s*(k-1)) + del*(w(m+s*(k-1))<w(m+s*k)) -del*(w(m+s*(k-1))>w(m+s*k)),boxw); % then it is just admissible
      cut = cut + pre - w(m+s*k);
      k = k+1;      
    else % already admissible, nothing to correct
      k = N_distribute + 1; % end the loop
    end        
  end
  if m+s*N_distribute<=2
    if abs(w(2))>del+1e-12
      skip = true;
    end
  elseif m+s*N_distribute>=d.Nu-1
    if abs(w(end-1))>del+1e-12
      skip = true; 
    end        
  else
    if abs(w(m+s*N_distribute)-w(m+s*(N_distribute+1)))>del+1e-12
      skip = true; % atm this point should not be cutted!
    end
  end
end
