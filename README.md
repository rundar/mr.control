mr.control
-------
Magnetic resonance RF pulse design by time optimal control

This software solves the time-optimal control problem introduced in the
paper
	"Simultaneous Multislice Refocusing via Time Optimal Control"
	[submitted, TBA]
by Armin Rund, Christoph Aigner, Karl Kunisch, and Rudolf Stollberger.
The code implements the diffusion examples.

Introduction
-------
The optimization is done jointly for (real-valued) RF and Gs amplitude, as well as
duration T. Included are the inequality constraints
     |RF|<=RF_max,  Gs<=Gs_max,   |Gs_slew|<= Gs_slew_max,
              d.dt*sum(|RF|.^2) <= RF_power_max

The lower-level solver tr_semism_quasiNewton() and objfun() are based on 
the paper
  "Magnetic Resonance RF pulse design by optimal control with physical
	constraints" [https://doi.org/10.1109/TMI.2017.2758391]
by the same authors.

The examples are taken from the ISMRM Challenge 2015/2016, see 
http://challenge.ismrm.org/node/71. Some initial guesses (RF and Gs) were 
created by using the PINS principle with a SLR subpulse, others by a 
superposition of the phase shifted SLR subpulses (RF) with a matched const
Gs shape, similar to the code supplied by the ISMRM challenge organizers.

Contents
-------

##### Test script (run this):
        main.m                              test script to start the optimization 

##### Routines called by the test script:
        time_optimal_control.m              implements the time optimal control method
        prep_init.m:                        problem definition
        optimization_parameters.m:          parameters of the time-optimal control method
        check_and_transform_input.m         basic input parameter check
        reduce_PINS_duration.m              shortens PINS pulses instantaneously by deleting most of the parts with Gs=0 
        reduce_duration_and_warmstart.m     cut points out of the time grid
        tr_semism_quasiNewton.p             implements trust region semismooth quasi Newtion
        objfun.p                            computes the objective function
        prepare_globalization.p             find good candidate time points for cut/reduction        
        globalization.p                     upper level method with globalization (further description in time_optimal_control.m)
        globalization_fast.p                fast variant of the globalization (further description in time_optimal_control.m)
        blochsd.m                           spin domain bloch simulation
        store_solution.m                    saves the optimized RF and Gs shape
        plot_results.m                      plots optimized RF, Gs, slew-rate of Gs and the refocusing profile

        
##### Data files used by the test scripts:        
        inits/CDIFF%_%.mat                  binary file containing a DIFF PINS init (challenge constraints)
        inits/CDIFF%_%_SUP.mat              binary file containing a DIFF superposition init (challenge constraints)
        inits/PDIFF1_4.mat                  binary file containing a DIFF PINS init (scanner 1 constraints)
        inits/SDIFF1_4.mat                  binary file containing a DIFF PINS init (scanner 2 constraints)
        
Usage
-------

The code is started with the MATLAB script main.m and was tested with MATLAB version R2016a (9.0.0.341360) 64-bit (glnxa64).
Parameters can be changed in 'prep_init.m' (problem definition) and 'optimization_parameters.m' (parameters of the time-optimal control method).
        
License
-------        
This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. To view a copy of this license, see license_CC_BY-NC.txt, visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
	 
mr.control is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
