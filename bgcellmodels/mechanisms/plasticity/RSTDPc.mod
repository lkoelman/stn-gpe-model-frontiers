TITLE Reward-modulated STDP weight adjuster, continuous time version

COMMENT

Reward-modulated STDP weight adjuster mechanism with continuously updated 
LTP/LTD traces, Eligibility trace, and DA concentration (as opposed to discrete
updates during spike processing). This is not efficient but useful for plotting 
in continuous time and understanding the learning mechanism.

Original STDP code adapted from:
http://senselab.med.yale.edu/modeldb/showmodel.asp?model=64261&file=\bfstdp\stdwa_songabbott.mod

Algorithm
---------

Adapted to implement reward-modulated STDP using eligibility traces 
as described in Izhikevich (2007) - 'Solving the distal reward problem 
through linkage of STDP and dopamine signaling'.

The STDP rule adjusts the eligibility trace instead of the synaptic weight. 
The synaptic weight update is the product of the reward (DA) signal 
and the eligibility variable.


Example
-------

>>> from neuron import h

Create cells:

>>> dummy = h.Section() # STDP mechanism can be put in any Section
>>> ncells = 2
>>> cells = []
>>> for c in range(ncells): cells.append(h.IntFire4(0,sec=dummy)) # Create the cells

Create synapses:

>>> threshold = 10 # Set voltage threshold
>>> delay = 1 # Set connection delay
>>> singlesyn = h.NetCon(cells[0],cells[1], threshold, delay, 0.5)

Create weight adjuster:

>>> stdpmech = h.RSTDPc(0, sec=dummy)
>>> presyn = h.NetCon(cells[0], stdpmech, threshold, delay, 1)
>>> pstsyn = h.NetCon(cells[1],stdpmech, threshold, delay, -1)
>>> h.setpointer(singlesyn._ref_weight[0], 'wsyn', stdpmech)

ENDCOMMENT

NEURON {
	POINT_PROCESS RSTDPc
	RANGE M, P, e_trace
	RANGE deltae, deltaw
	RANGE RLon
	GLOBAL emax, emin, wmax, wmin, DAmax, DAinf, aLTP, aLTD
	GLOBAL tauLTP, tauLTD, tau_e, tauDA : global to all instances of this mechanism
	POINTER wsyn : Pointer to the weight (in a NetCon object) to be adjusted
	RANGE DA	: NOTE: DA concentration can be  global variable for the mechanism
				: but this is more difficulty to set up and maintain consistent
				: NOTE: In Izhikevich (2007), DA is incremented by 0.5 upon reward
				: and decays exponentially as DA*0.995 every second
	: POINTER DA	: method 2: pointer to hoc var and update in custom proc advance()
	
}

ASSIGNED {
	deltae				: change in eligibility trace
	deltaw				: change in synaptic weight
	wsyn				: weight of the synapse
	t0					: time of last update
	DA					: global dopamine concentration

	: State variables (STATE but without DERIVATIVE block)
	M					: LTD/anti-Hebbian part of eligibility trace (trace for post-synaptic spikes) will be NEGATIVE
	P					: LTP/Hebbian part of eligibility trace (trace for pre-synaptic spikes) will be POSITIVE
	e_trace				: Synaptic eligibility trace (a.k.a. synaptic tag)
	delta				: Time step
}


INITIAL {
	M = 0
	P = 0
	e_trace = 0
	DA = 0
	deltae = 0
	deltaw = 0	
	delta = 0
	t0 = 0
}

PARAMETER {
	tauLTP  = 10	(ms)    : decay time for LTP part ( values from           )
	tauLTD  = 10	(ms)    : decay time for LTD part ( Song and Abbott, 2001 )
	tau_e 	= 100 	(ms)	: decay time constant for eligibility trace (NOTE: this determines how 'distal' the reward can be)
	tauDA	= 200	(ms)	: decay time for DA concentration (200 corresponds to DA*=0.995)
	wmax    = 1		: max value of synaptic weight
	wmin	= 0		: min value of synaptic weight
	emax	= 1		: max of eligibility trace variable
	emin	= -1		: min of eligibility trace variable
	DAmax	= 100 	: max DA concentration
	DAinf	= 0		: equilibrium DA concentration
	aLTP    = 0.001		: amplitude of LTP steps
	aLTD    = 0.00106	: amplitude of LTD steps (CONSTRAINT: area under LTD curve must be larger/LTD part of STDP must be dominant over LTP part)
	RLon	= 1		: turn on reinforcement learning
}

BREAKPOINT {
	delta = t-t0 : dt since last update

	: decay of PRE/PLUS trace and POST/MINUS trace
	P = P - delta*P/tauLTP
	M = M - delta*M/tauLTD

	: decay of eligibility trace
	e_trace = e_trace - delta*e_trace/tau_e

	: decay of DA concentration
	DA = DA - delta*(DA-DAinf)/tauDA

	: update weights based only on eligibility trace and reinforcement signal
	if (RLon) {
		wsyn = wsyn + DA*e_trace
		: saturation constraint
		if (wsyn > wmax) {
			wsyn = wmax
		}
		if (wsyn < wmin) {
			wsyn = wmin
		}
	}

	t0 = t : Reset last time so delta can be calculated in the next time step
}

NET_RECEIVE (w) {
	: see Scholarpedia article on STDP, section 'online implementation'
	: - the weight is increased at the moment of postsynaptic firing 
	:   by an amount that depends on the value of the trace x left by 
	:   the presynaptic spike
	: - Similarly, the weight is depressed at the moment of presynaptic 
	:   spikes by an amount proportional to the trace y left by previous 
	:   postsynaptic spikes.

	deltae = 0

	if (w >= 0) { : PRE SPIKE OCCURRED
		: update trace of pre-synaptic spikes (PLUS)
		P = P + aLTP

		: perform LTD update (MINUS increment)
		: dE negatively proportional to amount of post-synaptic activity left
		deltae = emax * M

	} else {      : POST SPIKE OCCURRED
		: update trace of post-synaptic spikes (MINUS)
		M = M - aLTD

		: perform LTP update (PLUS increment)
		: dE positively proportional to amount of pre-synaptic activity left
		deltae = emax * P

	}

	: eligibility trace adjustment
	e_trace = e_trace + deltae
	: saturation constraint
	if (e_trace > emax) {
		e_trace = emax
	}
	if (e_trace < emin) {
		e_trace = emin
	}

	
}

PROCEDURE reward_punish(reinf) {
	: increment DA concentration
	DA = DA + reinf
	if (DA > DAmax) { DA = DAmax }
}
