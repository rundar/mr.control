function [d,q] = reduce_PINS_duration(d,q)
%% This routine shortens PINS pulses instantaneously by deleting most of the parts with Gs=0.
%     written by Armin Rund, Dec 18, 2017
%############## Copyright: ##################################
%  Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                Christoph Aigner, christoph.aigner@tugraz.at
% Implementation of the time optimal control method of paper
% "Simultaneous Multislice Refocusing via Time Optimal Control"
% by A. Rund, C. Aigner, K. Kunisch, and R. Stollberger.
%%###########################################################
  u_in = min(d.box,max(-d.box,q));	
  u = u_in(1:d.Nu);
  slew = u_in(d.Nu+1:2*d.Nu-2);  
  w = zeros(d.Nu,1); 
  for n=2:d.Nu-1
     w(n) = w(n-1) + d.dt*slew(n-1);
  end
  id=find(w==0);
  i2=find(diff(id)>1);
  iR=[id(i2);d.Nu];
  iL=[1;id(i2+1)];
  % then w(iL) is first point in a vale, w(iR) last 
  N_vale = length(iR);
  for n=N_vale:-1:1  % backwards to allow for direct deletion of time instances
    RFsum=sum(u(iL(n):iR(n)));
    nt = ceil(abs(RFsum)/d.u_box); % number of needed time instances
    NT = iR(n)-iL(n)+1; % number of time instances in the vale
    if (nt>=NT), continue; end;
    u(iL(n):iL(n)+nt-1)=RFsum/nt;
    % cut the unneccessary right part of the vale:
    u(iL(n)+nt:iR(n))=[];
    w(iL(n)+nt:iR(n))=[];    
  end
  % finish the cuts:
  slew=diff(w)/d.dt; 
  slew=slew(1:end-1);
  q = [u;slew];
  d.Nu = size(u,1); d.T = d.Nu*d.dt; 
  d.box=[d.box(1)*ones(d.Nu,1);d.box(end)*ones(d.Nu-2,1)]; 
end
