function [RF,Gs,tdis,T] = time_optimal_control(init)
%     written by Armin Rund, Dec 18, 2017
%############## Copyright: ##################################
%  Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                Christoph Aigner, christoph.aigner@tugraz.at
% Implementation of the time optimal control method of paper
% "Simultaneous Multislice Refocusing via Time Optimal Control"
% by A. Rund, C. Aigner, K. Kunisch, and R. Stollberger.
%%###########################################################
tdis=0; 
fprintf('Starting time optimal control method for example %d at %s.\n',init.example,datestr(now));
%%###########################################################
% step-0 are initializations (checks, parameter calibration, PINS reduction), the time-optimal control starts from step-1
% output: tdis is the discrete time vector, T the terminal time/pulse duration
% input: see prep_init.m
%%###########################################################

%############################################################
%% 0. transform input parameters: ###########################
%############################################################
check_and_transform_input    % input parameter check and creating struct d
slew=diff(init.Gs)/d.dt;          % Gs slew rate
q=[init.RF;slew(1:end-1)];        % optimization variable 
    % (the last slew rate is not free but given by terminal condition Gs=0)
optimization_parameters      % introduce optimization parameters
T = d.T; T_init=T; 

fprintf('the optimization parameters are: p0=%d PINS=%d N_distribute=%d N_skip=%d N_glob=%d N_tries=%d rho=%1.1e\n',p0,d.is_PINS,N_distribute,N_skip,N_glob,N_tries,d.rho); 

%###########################################
%% 0.1. check initialization:  ###############
%###########################################
q0=q; d0=d;
o1=objfun(q,d); err = o1.err_o+o1.err_i;
blank=' '; if o1.max_g<10 && d.w_box>10, blank=''; end;
fprintf('\n maximum of:  |Gs| |Gs_slew| err_out err_in RF_power|  T\n'); 
fprintf( 'init:          %2.1f   %3.1f   %1.3f   %1.3f   %1.1e | %1.3f\n', ...
     o1.max_g,o1.max_slew,o1.err_o,o1.err_i,o1.ccosts_rf,d.T);
fprintf('constraints: %s %2.1f   %3.1f   %1.3f   %1.3f   %1.1e |\n',...
     blank,d.w_box,d.slew_box,d.maxErr_out,d.maxErr_in,d.RF_power_max);
init_admiss = (o1.max_g<=d.w_box & o1.last_slew<=d.slew_box & o1.err_o<=d.maxErr_out & o1.err_i<=d.maxErr_in & o1.ccosts_rf<=d.RF_power_max);
admissible = init_admiss;

%###########################################
%% 0.2. Instantaneous reduction of PINS INITS:
%###########################################
T_P = 0;
if d.is_PINS
  T0=d.T; tic;
  [d,q]=reduce_PINS_duration(d,q); t1=toc;
  T_P=d.T;
  fprintf('step 0.2: modifying PINS: reducing T from %1.2f to %1.2f  in %1.2f min  %s\n',T0,d.T,t1/60); 
  o2=objfun(q,d); err = o2.err_o+o2.err_i;
  fprintf('PINS reduced:  %2.1f   %3.1f   %1.3f   %1.3f   %1.1e | %1.3f\n', ...
     o2.max_g,o2.max_slew,o2.err_o,o2.err_i,o2.ccosts_rf,d.T);
  admissible = (o2.max_g<=d.w_box & o2.last_slew<=d.slew_box & o2.err_o<=d.maxErr_out & o2.err_i<=d.maxErr_in & o2.ccosts_rf<=d.RF_power_max);   
end

%###########################################
%% 0.3. check init and if needed get it admissible:
%###########################################
if ~admissible
  fprintf('The initial pulse is not admissible to the constraints. Try to fix this by applying the optimal control method...\n');
  maxit = 500; % maximum number of TR iterations
  it_adapt_param = 5; % calibrate the parameters every such TR step
  [d,q]=calibrate_parameters(1.0,d,q,maxit,it_adapt_param); % calibrate the objective parameters
  o=objfun(q,d); err = o.err_o+o.err_i;
  fprintf('..after optim: %2.1f   %3.1f   %1.3f   %1.3f   %1.1e | %1.3f\n', ...
      o.max_g,o.max_slew,o.err_o,o.err_i,o.ccosts_rf,d.T);
  admissible = (o.max_g<=d.w_box & o.last_slew<=d.slew_box & o.err_o<=d.maxErr_out & o.err_i<=d.maxErr_in & o.ccosts_rf<=d.RF_power_max);
  if admissible
    fprintf('successful.\n');
  else
    if d.is_PINS && init_admiss
      fprintf('not successful. Revoking the PINS reduction.\n'); 
      q=q0; d=d0; 
      clear q0 d0
      err = o1.err_o+o1.err_i;
    else  
      fprintf('not successful. Please design a better initial guess or loosen the constraints.\n'); return;
    end
  end
end

%###########################################
%% 0.4. Calibrate the obj parameters:
%###########################################
%% apply the automatic parameter calibration once at the beginning
%% do 200 steps of parameter calibration
  err0 = err; qold=q; dold=d;
  [d,q]=calibrate_parameters(0.85,d,q,200,4); % calibrate the slice error to 85 percent
  % check admissibility:
  o=objfun(q,d); err = o.err_o+o.err_i;
  admissible = (o.max_g<=d.w_box & o.last_slew<=d.slew_box & o.err_o<=d.maxErr_out & o.err_i<=d.maxErr_in & o.ccosts_rf<=d.RF_power_max);
  if admissible
    fprintf('...calibrated: %2.1f   %3.1f   %1.3f   %1.3f   %1.1e | %1.3f\n', ...
      o.max_g,o.max_slew,o.err_o,o.err_i,o.ccosts_rf,d.T);
  else 
    fprintf('warning: parameter calibration not successful. Continuing with time-optimal control.\n')
    q = qold; d=dold; err=err0;
  end
  clear qold dold

%###########################################
%% 1. Starting the time-optimal control method:
%###########################################
% always do the cheap globalization first:
tic; p1=p0; T0=d.T; d.pp=p0;
T00=d.T+1; keep_param=false; 
fprintf('T=%1.3f,  increasing p to %d\n',d.T,p1);
while d.T<T00-(N_red-0.5)*d.dt  % as long as we have enough decrease
  T00=d.T; 
  [d,q] = globalization_fast(d,q,N_glob,N_distribute,N_tries,N_skip,keep_param);
  %% for a description of globalization_fast() see below
end
while d.pp<p_max
  p1=min(p1+p_add,p_max);
  fprintf('T=%1.3f,  increasing p to %d\n',d.T,p1);
  d.pp=p1; 
  d.beta = d.beta*0.5;
  maxit = 150+2*d.pp; % maximum number of TR iterations
  it_adapt_param = 3+ceil(d.pp/10); % calibrate the parameters every such TR step
  [d,q]=calibrate_parameters(0.85,d,q,maxit,it_adapt_param); % calibrate the slice error to 90 percent
  T00=d.T+1; keep_param=false; 
  while d.T<T00-(N_red-0.5)*d.dt   % as long as we have decrease
    T00=d.T; 
    [d,q] = globalization_fast(d,q,N_glob,N_distribute,N_tries,N_skip,keep_param);
  end
  mx = max(d.xi1,d.xi2); if mx>1e20, d.xi1=d.xi1/mx*1e20; d.xi2=d.xi2/mx*1e20; end;
end
save('oc_phase1.mat','d','q');
% store_solution(q,d,1);   % to write out the RF and Gs at this stage
t1=toc/60; 
fprintf('... terminated in %1.2f min\n',t1);
fprintf('phase 1 with p=%d reducing T from %1.2f to %1.2f  in %1.2f min %s\n',d.pp,T0,d.T,t1,datestr(now));


% do a re-print of the table:
o=objfun(q,d); T1=d.T;
fprintf('\n maximum of:  |Gs| |Gs_slew| err_out err_in RF_power|  T\n'); 
fprintf( 'init:          %2.1f   %3.1f   %1.3f   %1.3f   %1.1e | %1.3f\n', ...
     o1.max_g,o1.max_slew,o1.err_o,o1.err_i,o1.ccosts_rf,d.T);
fprintf('constraints: %s %2.1f   %3.1f   %1.3f   %1.3f   %1.1e |\n',...
     blank,d.w_box,d.slew_box,d.maxErr_out,d.maxErr_in,d.RF_power_max);
fprintf('toc phase 1:    %2.1f   %3.1f   %1.3f   %1.3f   %1.1e | %1.3f\n', ...
   o.max_g,o.max_slew,o.err_o,o.err_i,o.ccosts_rf,d.T);



%###########################################
%% 2. Change to a more expensive globalitation
%###########################################
T2=0; t2=0;
if d.globalization_type>0
  fprintf('Phase 2: Changing to a more expensive globalization.\n'); 
  tic; T0=d.T; T00=d.T+1;
  maxit_per_try = 20;  
  N_tries = 3; p_max2=100;   
  maxit_cont = 100;  
  while d.pp<p_max2
    p1=min(p1+p_add,p_max2);
    fprintf('T=%1.3f,  increasing p to %d\n',d.T,p1);
    d.pp=p1;
    d.beta = d.beta*0.5;
    maxit = 150+2*d.pp; % maximum number of TR iterations
    it_adapt_param = 3+ceil(d.pp/10); % calibrate the parameters every such TR step
    [d,q]=calibrate_parameters(0.95,d,q,maxit,it_adapt_param); % calibrate the slice error to 95 percent
    maxit_per_try = min(20,maxit_per_try+1);  
    maxit_cont = min(200,maxit_cont+10);  
    T00=d.T+1; 
    while d.T<T00  % as long as we have decrease
      T00=d.T; 
      [d,q] = globalization(d,q,N_glob,N_distribute,maxit_per_try, N_tries, maxit_cont,N_skip);
      %% for a description of globalization() see below
    end
    mx = max(d.xi1,d.xi2); if mx>1e20, d.xi1=d.xi1/mx*1e20; d.xi2=d.xi2/mx*1e20; end;  
  end 
  t2=toc/60; fprintf('... terminated in %1.2f min\n',t2);
  fprintf('phase 2 with p=%d reducing T from %1.2f to %1.2f  in %1.2f min \n',d.pp, T0,d.T,t2); 
  o2 = objfun(q,d); T2=d.T;
  save('oc_phase2.mat','d','q');
end



%###########################################
%% 7. Finalize: Postprocessing
%###########################################
o = objfun(q,d);
RF=o.u;
Gs=o.w;
tdis=linspace(0,d.T,d.Nu+1);
T=d.T;
save('optimal_control.mat','d','q','RF','Gs','tdis','T');
% store_solution(q,d,1);   % to write out the RF and Gs at this stage
fprintf('\n maximum of:  |Gs| |Gs_slew| err_out err_in RF_power|  T\n'); 
fprintf( 'init:          %2.1f   %3.1f   %1.3f   %1.3f   %1.1e | %1.3f\n', ...
     o1.max_g,o1.max_slew,o1.err_o,o1.err_i,o1.ccosts_rf,d.T);
fprintf('constraints: %s %2.1f   %3.1f   %1.3f   %1.3f   %1.1e |\n',...
     blank,d.w_box,d.slew_box,d.maxErr_out,d.maxErr_in,d.RF_power_max);
fprintf('final opt:      %2.1f   %3.1f   %1.3f   %1.3f   %1.1e | %1.3f\n', ...
   o.max_g,o.max_slew,o.err_o,o.err_i,o.ccosts_rf,d.T);
end

function [d,q]=calibrate_parameters(err_fac,d,q,maxit,it_adapt_parameters)
%  maxit = 200;   % maximum number of trust-region steps
%  it_adapt_parameters = 4;  % parameter calibration every X steps. integer between 4 and 20, e.g. 4 for calibration and 20 for optimization
  % calibrate the parameters again: 
  % err_fac between 0.7 (far from minimum-duration) to  <1
  maxErr_in=d.maxErr_in0; maxErr_out=d.maxErr_out0;
  d.maxErr_in = err_fac*d.maxErr_in0; d.maxErr_out = err_fac*d.maxErr_out0;
  %% optimization call to the lower-level solver (trust-region semismooth Newton/quasi-Newton method from paper https://doi.org/10.1109/TMI.2017.2758391)
  d.adapt_param = 1;
  [q,success,d,err] = tr_semism_quasiNewton(d,@objfun,q,maxit,it_adapt_parameters); 
  d.adapt_param = 0;
  d.maxErr_in = maxErr_in; d.maxErr_out = maxErr_out;
end


%% Descriptions of the globalizaton methods:
%%########################################################################
% function [d,q] = globalization_fast(d,q,N_glob,N_distribute,N_tries,N_skip,keep_param)
% fast variant of the globalization
% upper level method with globalization. First, prepare_globalization()
% yields a list of good candidate time points for cut/reduction. Then, this
% list is used for reduction (accepting/rejecting certain cuts according to
% performance). Here, the leading N_tries cuts are tested, the best one is
% accepted (if admissible), the list is reduced by an accepted cut and its
% N_skip neighbors per side. 
%%%%%%%%%%%  Input/Output:
% d (I/O): struct of data and parameters
% q (I/O): optimization variable (RF amplitude and Gs slew rate)
% N_glob:  prepare_globalization() yields the best N_glob time points for 
%          possible cut
% N_distribute: At cut and warmstart, the cut RF/Gs is distributed locally
%          to up to N_distribute neighbors to the left and right
% N_tries: number of cuts to be tested, the best is accepted (if admissible)
% N_skip:  At an accepted cut, the N_skip neighboring time points are 
%          removed from the list of candidate cuts, since possibly outdated
% keep_param: if false, a new calibration of the parameters is done
%%########################################################################


%%########################################################################
% function [d,q] = globalization(d,q,N_glob,N_distribute,maxit_per_try, N_tries, maxit_cont,N_skip)
% upper level method with globalization. First, prepare_globalization()
% yields a list of good candidate time points for cut/reduction. Then, this
% list is used for reduction (accepting/rejecting certain cuts according to
% performance). Here, the leading N_tries cuts are tested for maxit_per_try 
% trust-region (TR) iterations, the best one is continued for further 
% maxit_cont TR iterations and accepted (if admissible). If accepted, the
% list of best cuts is reduced by it and its N_skip neighbors per side. 
%% Input/Output:
% d (I/O): struct of data and parameters
% q (I/O): optimization variable (RF and Gs slew rate)
% N_glob:  prepare_globalization() yields the best N_glob time points for 
%          possible cut
% N_distribute: At cut and warmstart, the cut RF/Gs is distributed locally
%          to up to N_distribute neighbors to the left and right
% maxit_per_try: number of TR iterations for each N_tries different cuts
% N_tries: number of cuts to be tested, the best is accepted (if admissible)
% maxit_cont: the best cut is continued for this number of TR iterations
% N_skip:  At an accepted cut, the N_skip neighboring time points are 
%          removed from the list of candidate cuts, since possibly outdated
%%########################################################################
  

  


