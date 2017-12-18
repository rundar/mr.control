%########## Description: #############################################
% This software solves the time-optimal control problem introduced in the
% paper
%  "Simultaneous Multislice Refocusing via Time Optimal Control"
% [submitted, TBA]
% by Armin Rund, Christoph Aigner, Karl Kunisch, and Rudolf Stollberger.
% Here, the code for the diffusion examples is published. 
%
% The optimization is done jointly for (real-valued) RF and Gs amplitude, as well as
% duration T. Included are the inequality constraints
%     |RF|<=RF_max,  Gs<=Gs_max,   |Gs_slew|<= Gs_slew_max,
%              d.dt*sum(|RF|.^2) <= RF_power_max
%
% The lower-level solver tr_semism_quasiNewton() and objfun() are based on 
% the paper
%  "Magnetic Resonance RF pulse design by optimal control with physical
% constraints" [https://doi.org/10.1109/TMI.2017.2758391]
% by the same authors.
%
% The examples are taken from the ISMRM Challenge 2015/2016, see 
% http://challenge.ismrm.org/node/71. Some initial guesses (RF and Gs) were 
% created by using the PINS principle with a SLR subpulse, others by a 
% superposition of the phase shifted SLR subpulses (RF) with a matched const
% Gs shape, similar to the code supplied by the ISMRM challenge organizers.
%
% The code was tested with MATLAB version R2016a (9.0.0.341360) 64-bit (glnxa64).
%############## Copyright: ##########################################
%     Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                         Christoph Aigner, christoph.aigner@tugraz.at
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial 
% 4.0 International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% 
% This work is distributed WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% 
%#####################################################################
% clear all
addpath('inits/');

%% the required input data is described in prep_init.m
%% Decide the example: 
% 11,...,15,21,...,25,31,...,35 are the DIFF examples with PINS
% add 100 to get the same example with a superposition initial
% 1,2,14 are the examples from the paper, Figure 3
% choose 3 for customized parameters (to enter in prep_init.m)
% example = 15;  
if ~exist('example'), 
  example = 15;  % default case
end
%example = example+100 % to switch from PINS to superposition initial

%% prepare initial structure
init = prep_init(example);
init.globalization_type = 1;    % refinement level of globalization: 0 (fastest, less effective) or 1 (more expensive, more effective)

%% start the optimization
[RF,Gs,tdis,T] = time_optimal_control(init);

%% plot the results
%  plot_results(RF, Gs, init.dt, init);













  

  


