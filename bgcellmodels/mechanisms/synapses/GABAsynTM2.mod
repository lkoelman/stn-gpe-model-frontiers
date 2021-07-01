TITLE GABAA and GABAB receptor with independent Tsodyks-Markram short-term dynamics

COMMENT
GABA-A and GABA-B receptor conductance using a dual-exponential profile
and independent Tsodyks-Markram short-term dynamics
ENDCOMMENT

: Definition of variables which may be accessed by the user 
NEURON {
    THREADSAFE

    POINT_PROCESS GABAsynTM2
    
    RANGE tau_r_GABAA, tau_d_GABAA
    RANGE tau_r_GABAB, tau_d_GABAB
    RANGE gmax_GABAA, gmax_GABAB
    RANGE i, i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, Erev_GABAA, Erev_GABAB
    RANGE tau_rec_A, tau_facil_A, U1_A
    RANGE tau_rec_B, tau_facil_B, U1_B
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (umho) = (micromho)
    (mM) = (milli/liter)
}


: Definition of constants which may be set by the user
PARAMETER {

    tau_r_GABAA = 0.2   (ms)  : Dual-exponential conductance profile
    tau_d_GABAA = 1.7   (ms)  : IMPORTANT: tau_r < tau_d
    tau_r_GABAB = 5    (ms)  : Dual-exponential neurotransmitter profile
    tau_d_GABAB = 25   (ms)  : IMPORTANT: tau_r < tau_d

    Erev_GABAA = -80    (mV)  : GABAA reversal potential
    Erev_GABAB = -95    (mV)  : GABAB reversal potential
    
    : NOTE: gmax can be set either via NetCon.weight[0] in [nS] or here in [uS] (assuming that weight is 1)
    gmax_GABAA = .001   (uS)  : Weight conversion factor (from nS to uS)
    gmax_GABAB = .001   (uS)  : Weight conversion factor (from nS to uS)
                              : (typically about 70% of GABA-A conductance)

    : Parameters of TM model
    tau_rec_A = 5        (ms)  : time constant of recovery from depression
    tau_facil_A = 100    (ms)  : time constant of facilitation
    U1_A = 0.05           (1)   : baseline release probability
    tau_rec_B = 5        (ms)  : time constant of recovery from depression
    tau_facil_B = 100    (ms)  : time constant of facilitation
    U1_B = 0.05           (1)   : baseline release probability
}

: Declaration of state variables 
STATE {

    A_GABAA       : GABAA state variable to construct the dual-exponential profile
                 : rise kinetics with tau_r_GABAA

    B_GABAA       : GABAA state variable to construct the dual-exponential profile
                 : decay kinetics with tau_d_GABAA

    A_GABAB       : GABAB state variable to construct the dual-exponential profile
                 : rise kinetics with tau_r_GABAB

    B_GABAB       : GABAB state variable to construct the dual-exponential profile
                 : decay kinetics with tau_d_GABAB

: State variables for the TM model

    Rrp_A          : Readily releaseable pool: fraction of available vesicles => STD
    Use_A          : Synaptic efficacy: fraction of RRP released per spike => STP
    Rrp_B
    Use_B


}

: Declaration of variables that are computed, e.g. in the BREAKPOINT block
ASSIGNED {
    v (mV)
    i (nA)
    i_GABAA (nA)
    i_GABAB (nA)
    g_GABAA (uS)
    g_GABAB (uS)
    g (uS)
    factor_GABAA
    factor_GABAB
}

: Definition of initial conditions
INITIAL{
    LOCAL tp_GABAA, tp_GABAB     : Declaration of some local variables

    : Zero receptor rise and fall kinetics variables
    A_GABAA = 0
    B_GABAA = 0

    A_GABAB = 0
    B_GABAB = 0

    : Compute constants needed to normalize the dual-exponential receptor dynamics

    : Expression for time to peak of the GABAA dual-exponential conductance
    tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA)
    : Expression for time to peak of the GABAB dual-exponential conductance
    tp_GABAB = (tau_r_GABAB*tau_d_GABAB)/(tau_d_GABAB-tau_r_GABAB)*log(tau_d_GABAB/tau_r_GABAB)

    : GABAA Normalization factor - so that when t = tp_GABAA, gsyn = gpeak
    :   (i.e. a weight of 1 will cause a peak in B_GABAA-A_GABAA of 1 at time tp_GABAA)
    factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA) 
    factor_GABAA = 1/factor_GABAA
    : GABAB Normalization factor - so that when t = tp_GABAB, gsyn = gpeak
    :   (i.e. a weight of 1 will cause a peak in B_GABAB-A_GABAB of 1 at time tp_GABAB)
    factor_GABAB = -exp(-tp_GABAB/tau_r_GABAB)+exp(-tp_GABAB/tau_d_GABAB) 
    factor_GABAB = 1/factor_GABAB

    : Synaptic facilitation and depression
    Rrp_A = 1
    Use_A = 0
    Rrp_B = 1
    Use_B = 0

}

: Declare method to propagate the state variables in time
BREAKPOINT {

    : Specify to solve system of equations "odes", declared below (DERIVATIVE block)
    : "cnexp" specifies the intergration method, it is
    : an implicit integration method that is stable even for stiff systems
    SOLVE odes METHOD cnexp

    : Compute and assign quantities which depend on the state variables

    : Compute the time varying GABAA receptor conductance as 
    : the difference of state variables B_GABAA and A_GABAA
    g_GABAA = gmax_GABAA*(B_GABAA-A_GABAA) 
    g_GABAB = gmax_GABAB*(B_GABAB-A_GABAB)
    g = g_GABAA + g_GABAB

    : Compute the GABAA and GABAB specific currents
    i_GABAA = g_GABAA * (v-Erev_GABAA) 
    i_GABAB = g_GABAB * (v-Erev_GABAB) 

    : Compute the total current
    i = i_GABAA + i_GABAB
}

: Declaration of ODEs solved for in the BREAKPOINT block
DERIVATIVE odes {
    : State variables for bi-exponential profile
    A_GABAA' = -A_GABAA/tau_r_GABAA
    B_GABAA' = -B_GABAA/tau_d_GABAA
    A_GABAB' = -A_GABAB/tau_r_GABAB
    B_GABAB' = -B_GABAB/tau_d_GABAB

    : State variables for synaptic facilitation & depression
    Rrp_A' = (1-Rrp_A)/tau_rec_A
    Use_A' = -Use_A/tau_facil_A

    Rrp_B' = (1-Rrp_B)/tau_rec_B
    Use_B' = -Use_B/tau_facil_B
}


: Block to be executed for a pre-synaptic spike event
NET_RECEIVE (weight) {
    LOCAL rel_A, rel_B

    : Use = 'running release probability' = facilitation variable
    :     => The longer (slower) tau_facil, the more easily Use builds up.
    : U1  = baseline release probability
    :     = increment in Use per spike as fraction of what's left of Use.
    :     => The higher U1, the faster Use build up (increment per spike)
    
    : Increment release probability.
    : (spike-triggered delta-term in diff. eq. for Use)
    Use_A = Use_A + U1_A*(1-Use_A)
    Use_B = Use_B + U1_B*(1-Use_B) 
    
    : A is product of (size of RRP, i.e. available vesicles) * (their release probability).
    : This corresponds to the amount of vesicles released.
    : That's why A functions both as the decrement for R and the increment for gsyn.
    rel_A = Use_A*Rrp_A
    rel_B = Use_B*Rrp_B
  
    : Decrease size of RRP (spike-triggered delta-term in diff. eq. for R)
    Rrp_A = Rrp_A - rel_A
    Rrp_B = Rrp_B - rel_B
    
    A_GABAA = A_GABAA + rel_A*weight*factor_GABAA
    B_GABAA = B_GABAA + rel_A*weight*factor_GABAA

    A_GABAB = A_GABAB + rel_B*weight*factor_GABAB
    B_GABAB = B_GABAB + rel_B*weight*factor_GABAB
}
