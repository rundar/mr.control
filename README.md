# mr.control
Magnetic resonance RF pulse design by time optimal control

Introduction
------------

This software solves the time-optimal control problem introduced in the
paper
	"Simultaneous Multislice Refocusing via Time Optimal Control"
	[submitted, TBA]
by Armin Rund, Christoph Aigner, Karl Kunisch, and Rudolf Stollberger.
The code implements the diffusion examples.

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

The code was tested with MATLAB version R2016a (9.0.0.341360) 64-bit (glnxa64).

Usage
-----

The code is started with the MATLAB script main.m.
Parameters can be changed in 
    prep_init.m: problem definition
		optimization_parameters.m: parameters of the time-optimal control method
