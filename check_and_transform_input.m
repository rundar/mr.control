%     written by Armin Rund, Dec 18, 2017
%############## Copyright: ##################################
%  Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                Christoph Aigner, christoph.aigner@tugraz.at
% Implementation of the time optimal control method of paper
% "Simultaneous Multislice Refocusing via Time Optimal Control"
% by A. Rund, C. Aigner, K. Kunisch, and R. Stollberger.
%%###########################################################
%% please do not change this file. All changeable optimization 
%  parameters are in optimization_parameters.m
%#####################################################
%% basic input parameter check: ######################
%#####################################################
if init.Gs_max<=0,       fprintf('error in input parameter Gs_max: must be positive! Taking it to be 80.\n'); init.Gs_max = 80; end;
if init.Gs_slew_max<=0,  fprintf('error in input parameter Gs_slew_max: must be positive! Taking it to be 180.\n'); init.Gs_max = 180; end;
if init.RF_max<=0,       fprintf('error in input parameter RF_max: must be positive! Taking it to be 0.018.\n'); init.RF_max = 0.018; end;
if init.maxErr_in<=0,    fprintf('error in input parameter maxErr_in: must be positive! Taking it to be 0.03.\n'); init.maxErr_in = 0.03; end;
if init.maxErr_out<=0,   fprintf('error in input parameter maxErr_out: must be positive! Taking it to be 0.02.\n'); init.maxErr_out = 0.02; end;
if init.RF_power_max<=0, fprintf('error in input parameter RF_power_max: must be positive! Taking 7.68e-4.\n'); init.RF_power_max = 7.68e-4; end;
if init.fac_RFpower<=0,  fprintf('error in input parameter fac_RFpower: must be >0 and <=1! Taking it to be 1.\n'); init.fac_RFpower = 1; end;
if size(init.RF,1)==1, init.RF=init.RF'; end; if size(init.Gs,1)==1,init.Gs=init.Gs'; end;
if length(init.Gs)~= length(init.RF),fprintf('error in input parameter RF/Gs: RF and Gs must have the same length.\n');return; end;
if length(init.RF)<10 || length(init.RF)>10000, fprintf('too few or too many time samples. please choose between 10 and 10000 samples.\n'); return; end;
if length(init.RF)>5000, fprintf('Warning: >5k time samples. Better use a coarser initial time grid.\n'); end;
if length(init.xdis)>10000, fprintf('too many spatial points. please choose below 10000 points in xdis.\n'); return; end;
if size(init.xdis,1)>1 && size(init.xdis,2)>1, fprintf('error in input xdis: needs to be a vector, not a matrix.\n'); end;
if size(init.xdis,1)>1, init.xdis=init.xdis'; end;
if size(init.inslice,2)==1, init.inslice=init.inslice'; end; if size(init.outslice,2)==1, init.outslice=init.outslice'; end; 
if (size(init.inslice,2)~=size(init.xdis,2))||(size(init.outslice,2)~=size(init.xdis,2)), fprintf('error in input: inslice and outslice must be row vectors of the length of xdis.\n'); end; 
if sum(init.inslice.*init.outslice)>0, fprintf('error in input: inslice and outslice are not disjoint.\n'); return; end;
if init.dt<=0, fprintf('error in input: temporal step size dt needs to be positive.\n'); return; end;
if init.gamma<0, fprintf('error in input: gamma needs to be positive.\n'); return; end; 

%#####################################################
%% transform parameters to struct d: #################
%#####################################################
d.Nu = length(init.RF);    % number of samples
d.T = d.Nu*init.dt;        % current pulse duration
d.u_box = init.RF_max; 
d.w_box = init.Gs_max; 
d.slew_box = init.Gs_slew_max-1e-12;   % allow round-off errors here
d.maxErr_in = init.maxErr_in;
d.maxErr_out = init.maxErr_out;
d.maxErr_in0 = init.maxErr_in;
d.maxErr_out0 = init.maxErr_out;
d.RF_power_max = init.RF_power_max;
d.fac_SAR = init.fac_RFpower;
d.dx = init.xdis(2)-init.xdis(1); 
d.outslice = init.outslice;
d.inslice = init.inslice;
d.Nx = size(init.xdis,2);
d.dt = init.dt;
d.xdis = init.xdis;
d.box = [d.u_box*ones(d.Nu,1); d.slew_box*ones(d.Nu-2,1)]; % symmetric box constraints for controls
d.gamma = init.gamma;
d.is_PINS = init.is_PINS; 
d.adapt_param = 0;
d.example = init.example;
d.globalization_type = init.globalization_type;
     


