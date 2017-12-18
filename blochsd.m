function [a,b] = blochsd(RF,Gs,dt,dis,gamma)
% Bloch solver in the spin domain 
  Nx=size(dis,2); %spatial points
  Nu=size(RF,1);  %temporal points
  %% Bloch simulation in spin domain:
  a = ones(Nx,Nu); b=zeros(Nx,Nu);  
  % assemble alpha, beta vectorwise:
  gadt = gamma*dt;    
  B = repmat(gadt*RF', Nx,1);
  K = gadt*dis'*Gs';
  phi = -sqrt(abs(B).^2+K.^2);
  sinp = sin(0.5*phi)./phi; ind = isinf(sinp)|isnan(sinp); sinp(ind) = 0.5; 
  alpha = cos(0.5*phi)+1i*K.*sinp;
  beta = 1i*B.*sinp;
  % initialized by zero-rotation: (a=1,b=0)
  a(:,1) = alpha(:,1);
  b(:,1) = beta(:,1);
  for m=2:Nu  % time loop
    a(:,m) = alpha(:,m).*a(:,m-1)-conj(beta(:,m)) .*b(:,m-1);
    b(:,m) = beta(:,m) .*a(:,m-1)+conj(alpha(:,m)).*b(:,m-1);
  end
  %save only last timepoint
  a = a(:,end);
  b = b(:,end);
end