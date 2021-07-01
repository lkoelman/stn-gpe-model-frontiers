TITLE GABAA receptor with TM short-term dynamics

COMMENT

Brief Description:

    GABAA receptor conductance using a dual-exponential profile
    and incorporating Tsodyks-Markram short-term dynamics


Authors:

    - implemented by Lucas Koelman
    
    - Tsodyks-Markram model based on `TsodyksMarkram_AMPA_NMDA.mod` by E. Muller at BBP/EPFL

Date:
    
    - first version on 2017-07-20

USAGE
-----

- STP parameters and imax_gabao must be calibrated together, since STP parameters
  determine magnitude of initial peak

ENDCOMMENT

: Definition of variables which may be accessed by the user 
NEURON {
    THREADSAFE
    POINT_PROCESS syn_gabadiff

    USEION gaba WRITE igaba
    NONSPECIFIC_CURRENT i
    
    RANGE tau_r_GABAA, tau_d_GABAA
    RANGE gmax_GABAA
    RANGE i, i_GABAA, g_GABAA, g, Erev_GABAA
    RANGE tau_rec, tau_facil, U1, use_stdp_A
    RANGE imax_gabao
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
    Erev_GABAA = -80    (mV)  : GABAA reversal potential
    gmax_GABAA = .001   (uS)  : Weight conversion factor (from nS to uS)
    : NOTE: gmax can be set either via NetCon.weight[0] in [nS] or here in [uS] (assuming that weight is 1)

    : Parameters of TM model
    tau_rec = 200      (ms)  : time constant of recovery from depression
    tau_facil = 200    (ms)  : time constant of facilitation
    U1 = 0.5           (1)   : baseline release probability

    use_stdp_A = 1     (1)   : flag indicating whether STD/STP is enabled for GABA-A

    : Neurotransmitter release to extracellular pool
    imax_gabao = 1.0   (mM/ms)
}

: Declaration of state variables 
STATE {

    A_GABAA       : GABAA state variable to construct the dual-exponential profile
                 : rise kinetics with tau_r_GABAA

    B_GABAA       : GABAA state variable to construct the dual-exponential profile
                 : decay kinetics with tau_d_GABAA

: State variables for the TM model

    Rrp          : Running fraction of available vesicles (vesicles in RRP)
                 : Accounts for synaptic depression
    
    Use          : Running value of the release probability ('utilization of synaptic efficacy')
                 : Accounts for synaptic facilitation
}

: Declaration of variables that are computed, e.g. in the BREAKPOINT block
ASSIGNED {
    v (mV)
    i (nA)
    i_GABAA (nA)
    g_GABAA (uS)
    g (uS)
    factor_GABAA
    igaba (mM/ms)
}

: Definition of initial conditions
INITIAL{
    LOCAL tp_GABAA     : Declaration of some local variables

    : Zero receptor rise and fall kinetics variables
    A_GABAA = 0
    B_GABAA = 0

    : Compute constants needed to normalize the dual-exponential receptor dynamics
    : Expression for time to peak of the GABAA dual-exponential conductance
    tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA)

    : GABAA Normalization factor - so that when t = tp_GABAA, gsyn = gpeak
    :   (i.e. a weight of 1 will cause a peak in B_GABAA-A_GABAA of 1 at time tp_GABAA)
    factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA) 
    factor_GABAA = 1/factor_GABAA

    : Synaptic facilitation and depression
    Rrp = 1
    Use = 0

    : Neurotransmitter release to extracellular pool
    igaba = 0

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

    : Total conductance
    g = g_GABAA

    : GABA transmitter release
    igaba = imax_gabao*(B_GABAA-A_GABAA) 

    : Compute the GABAA specific currents
    i_GABAA = g_GABAA * (v-Erev_GABAA) 

    : Compute the total current
    i = i_GABAA
}

: Declaration of ODEs solved for in the BREAKPOINT block
DERIVATIVE odes {
    : State variables for bi-exponential profile
    A_GABAA' = -A_GABAA/tau_r_GABAA
    B_GABAA' = -B_GABAA/tau_d_GABAA

    : State variables for synaptic facilitation & depression
    Rrp' = (1-Rrp)/tau_rec
    Use' = -Use/tau_facil
}

: Block to be executed for a pre-synaptic spike event
NET_RECEIVE (weight) {
    LOCAL A
    
    Use = Use + U1*(1-Use) : Increment release probability.
                           : (spike-triggered delta-term in diff. eq. for Use)
    
    A = Use*Rrp : A is product of (size of RRP, i.e. available vesicles) * (their release probability).
              : This corresponds to the amount of vesicles released.
              : That's why A functions both as the decrement for R and the increment for gsyn.
    
    Rrp = Rrp - A : Decrease size of RRP (spike-triggered delta-term in diff. eq. for R)
    
    : Enable STD/STP for GABA-A
    if (use_stdp_A == 0) {
        A = 1
    }
    A_GABAA = A_GABAA + A*weight*factor_GABAA
    B_GABAA = B_GABAA + A*weight*factor_GABAA
}
