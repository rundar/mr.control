function init = prep_init(example)
%     written by Armin Rund and Christoph Aigner, Dec 18, 2017
%############## Copyright: ##################################
%  Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                Christoph Aigner, christoph.aigner@tugraz.at
% Implementation of the time optimal control method of paper
% "Simultaneous Multislice Refocusing via Time Optimal Control"
% by A. Rund, C. Aigner, K. Kunisch, and R. Stollberger.
%%###########################################################
%% load the initial pulse container or provide the following:
% init.RF, init.Gs   column vector of RF and Gs values
% details: piecewise constant discretization on an equidistant time grid 
% with step size init.dt, identical grid for RF and Gs, Gs(1)=Gs(end)=0 is 
% necessary. In this code, the time step is kept constant throughout the optimization.
% Both RF and Gs should fulfill the amplitude constraints
%         |Gs|<=Gs_max, |Gs'|<=Gs_slew_max, |RF|<=RF_max
% other required parameters are described in case 'example==2' below.

% examples: 11,...,15,21,...,25,31,...,35 are the DIFF examples with PINS
% add 100 to get the same example with a superposition initial
% 1,2,14 are the examples from the paper, Figure 3
% ex 3 for customized parameters (see below)
if example == 1, load('inits/PDIFF_1_4.mat');
elseif example == 2, load('inits/SDIFF_1_4.mat');
elseif example == 11, load('inits/CDIFF_1_1.mat');
elseif example == 12, load('inits/CDIFF_1_2.mat');
elseif example == 13, load('inits/CDIFF_1_3.mat');
elseif example == 14, load('inits/CDIFF_1_4.mat');
elseif example == 15, load('inits/CDIFF_1_5.mat');
elseif example == 21, load('inits/CDIFF_2_1.mat');
elseif example == 22, load('inits/CDIFF_2_2.mat');
elseif example == 23, load('inits/CDIFF_2_3.mat');
elseif example == 24, load('inits/CDIFF_2_4.mat');
elseif example == 25, load('inits/CDIFF_2_5.mat');
elseif example == 31, load('inits/CDIFF_3_1.mat');
elseif example == 32, load('inits/CDIFF_3_2.mat');
elseif example == 33, load('inits/CDIFF_3_3.mat');
elseif example == 34, load('inits/CDIFF_3_4.mat');
elseif example == 35, load('inits/CDIFF_3_5.mat');
%%% same examples with a superposition initialization:
elseif example == 111, load('inits/CDIFF_1_1_SUP.mat');
elseif example == 112, load('inits/CDIFF_1_2_SUP.mat');
elseif example == 113, load('inits/CDIFF_1_3_SUP.mat');
elseif example == 114, load('inits/CDIFF_1_4_SUP.mat');
elseif example == 115, load('inits/CDIFF_1_5_SUP.mat');
elseif example == 121, load('inits/CDIFF_2_1_SUP.mat');
elseif example == 122, load('inits/CDIFF_2_2_SUP.mat');
elseif example == 123, load('inits/CDIFF_2_3_SUP.mat');
elseif example == 124, load('inits/CDIFF_2_4_SUP.mat');
elseif example == 125, load('inits/CDIFF_2_5_SUP.mat');
elseif example == 131, load('inits/CDIFF_3_1_SUP.mat');
elseif example == 132, load('inits/CDIFF_3_2_SUP.mat');
elseif example == 133, load('inits/CDIFF_3_3_SUP.mat');
elseif example == 134, load('inits/CDIFF_3_4_SUP.mat');
elseif example == 135, load('inits/CDIFF_3_5_SUP.mat');

end
% for PINS we may initially reduce many of the times with Gs=0 before starting the optimization (PINS reduction)
init.is_PINS = 1;   % 0=PINS, 1=superposition or other
if example>110 && example < 136 % superposition inits
  init.is_PINS = 0;  
end

if example==3  % create customized example and initialization
  %% set constraints on RF and Gs:
  % the initial RF pulse and Gs have to fulfill the constraints 
  % |Gs|<=Gs_max, |Gs'|<=Gs_slew_max, |RF|<=RF_max
  init.Gs_max      = 80;    % maximum Gs amplitude, >0, e.g. 80 mT/m
  init.Gs_slew_max = 200;   % maximum Gs slew rate, e.g. 180 mT/m/ms
  init.RF_max      = 0.018; % maximum RF amplitude, e.g. 0.018 mT
  init.RF_power_max = 7.68e-4; % maximum RF pulse power d.dt*sum(|RF|.^2), 
      % e.g. for ISMRM Challenge examples: DIFF 7.68e-4, TSE 2.3467e-4 (mT)^2 ms
     
  %% set allowed refocusing Errors:
  % the initial (RF,Gs) should give a good profile nearly fulfilling the 
  % profile constraints of perfectly crushed refocusing profiles abs(b.^2) 
  % using the Cayley-Klein formalism :
  init.maxErr_in  = 0.03; % maximum Error in-slice (normalized to 1 e.g. |b|^2>=1-maxErr_in)
                          % in [0,0.1] e.g. 0.03
  init.maxErr_out = 0.025; % maximum Error out-of-slice (normalized to 1)
                          % e.g. 0.02
                     
  %% set algorithmic parameters:
  init.fac_RFpower = 1;   % relation between minimum duration and minimum 
                          % RF power, in [0.1,1]: 1 is minimum duration 
                          % (RF power constraint loose or missing), try e.g. 0.3 if 
                          % the RF power constraint is tight
  %% set other parameters based on initial file: 
  init.xdis = linspace(-0.06, 0.06, 481);     % range of spatial positions in m
  init.outslice = logical([ones(186,1); zeros(109,1); ones(186,1)]); % spatial mask to define out-slice postions (same length as xdis, zeros and ones)
  init.inslice  = logical([zeros(229,1); ones(23,1); zeros(229,1)]);  % spatial mask to define in-slice positions (same length as xdis, zeros and ones)
  
  %%% for reading in Challenge examples we note:   
  %%% init.inslice = b2d' > 0.5; %b2d is a 1xNz vector that can be generated by gen_mb_eval_prof.m (ISMRM challenge example code, see http://challenge.ismrm.org/node/71)
  %%% init.outslice = (roi-b2d)' > 0.5; %b2d and roi are 1xNz vectors that can be generated by gen_mb_eval_prof.m (ISMRM challenge example code, see http://challenge.ismrm.org/node/71)
  
  init.dt   = 0.01;       % temporal resolution in ms  
  init.RF   = [0; sinc(-2:0.015:2)'*init.RF_max; 0]; % RF shape in mT (piecewise constant on an equidistant time grid)
  init.Gs = [0; ones(267,1)*2; 0];  % vector for Gs shape  (piecewise constant on the same equidistant time grid), mt/m
  init.gamma = 2*pi*42.58; % proton gyromagnetic ratio (rad/s)/uT
  init.is_PINS = 0; 

  % plot the results
%  plot_results(init.RF, init.Gs, init.dt, init);
end
init.example = example;

end 
