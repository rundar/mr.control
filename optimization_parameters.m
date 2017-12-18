%% Parameters of the time-optimal control method:
%     written by Armin Rund, Dec 18, 2017
%############## Copyright: ##################################
%  Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                Christoph Aigner, christoph.aigner@tugraz.at
% Implementation of the time optimal control method of paper
% "Simultaneous Multislice Refocusing via Time Optimal Control"
% by A. Rund, C. Aigner, K. Kunisch, and R. Stollberger.
%%###########################################################

%%##### Parameters of the warmstart: #####
% N_distribute: parameter of the warmstart. At cut and warmstart, the cut 
% RF/Gs is distributed locally to up to N_distribute neighbors to the left 
% and right
N_distribute = min(50,max(3,floor(d.Nu/200)));     % integer number between 3 and 50 depending on the refinement level 

%%##### Parameters of the globalization for phase 1: #####
% N_glob:  prepare_globalization() yields a list of the best N_glob time points for possible cut. 
N_glob = max(20,min(400,floor(d.Nu/5))); % integer between 20 and 400
% N_tries: this list is then used (and reused) for finding a next good cut. The first N_tries entries are tested, compared, deleted from the list, and the best is proceeced and accepted (if admissible)
N_tries = 2;   % integer between 1 and 5. 2-3 is good for fast globalization
% N_skip: At an accept also the N_skip neighboring times are removed from the list of candidate cuts, since they are possibly outdated. If not, they appear again in the next list from new globalization().
N_skip = floor(N_distribute/3);   % integer between 0 and 2*N_distribute+1
N_red = 2;  % per globalization reduce at least N_red time points successfully, otherwise p++, integer between 1 and 5

%%###### Parameters controlling the integration order p: #######
p0 = 8; d.pp=p0; % initial penalty exponent (even). take 2, 4, 8, or 16
p_max = 32; % maximum integration exponent p in phase-1, even integer e.g. 32
p_add = 6;  % increase exponent p by this additively (even integer between 2 and 16) 

%%##### Parameters of the lower level optimization: #####
d.rho = 1e-5;     % initial trust-region radius (1e-5 is good)
%% the following four are not to be changed:
d.xi1=1;          % initial outslice error parameter, self-calibrating
d.xi2=1;          % initial inslice error parameter,  self-calibrating
d.alpha = 1e-4;   % initial RF power parameter, self-calibrating
d.beta = d.alpha*1e-2*init.RF_max/init.Gs_max;  % slew rate, is driven to 0





  

  


