function store_solution(q,d,final)
%     written by Armin Rund, Dec 18, 2017
%############## Copyright: ##################################
%  Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                Christoph Aigner, christoph.aigner@tugraz.at
% Implementation of the time optimal control method of paper
% "Simultaneous Multislice Refocusing via Time Optimal Control"
% by A. Rund, C. Aigner, K. Kunisch, and R. Stollberger.
%%###########################################################
  % get RF and Gs:
  o = objfun(q,d);
  RF=o.u;
  Gs=o.w;
  tdis=linspace(0,d.T,d.Nu+1);
  T=d.T;
  if final
      if d.Nu>4000
          filename = strcat('control_T',num2str(d.T,'%2.3f'),'.mat');
      else    
          filename = strcat('control_T',num2str(d.T,'%2.2f'),'.mat');       
      end
  else    
      filename = 'current_best_control.mat';
  end
  save(filename,'q','d','RF','Gs','tdis','T');
end
