TITLE GABAA and GABAB receptor with TM short-term dynamics

COMMENT

Brief Description:

    GABAA and GABAB receptor conductance using a dual-exponential profile
    and incorporating Tsodyks-Markram short-term dynamics


Detailed Description:

    - For details on GABA-B signaling cascade and its parameters,
      see file gagab3.mod obtainable at http://cns.iaf.cnrs-gif.fr/alain_demos.html
      under package "Kinetic Models of Synaptic Transmission (zip format) "

Authors:

    - implemented by Lucas Koelman
    
    - Tsodyks-Markram model based on `TsodyksMarkram_AMPA_NMDA.mod` by E. Muller at BBP/EPFL
    
    - GABA-B signaling cascade model based on `gabab3.mod` by A. Destexhe

Date:
    
    - first version on 2017-07-20

ENDCOMMENT

: Definition of variables which may be accessed by the user 
NEURON {
    THREADSAFE

    POINT_PROCESS GABAsyn
    
    RANGE tau_r_GABAA, tau_d_GABAA
    RANGE tau_r_GABAB, tau_d_GABAB
    RANGE gmax_GABAA, gmax_GABAB
    RANGE i, i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, Erev_GABAA, Erev_GABAB
    RANGE tau_rec, tau_facil, U1, use_stdp_A, use_stdp_B

    RANGE R, D, G
    RANGE K1, K2, K3, K4, KD, d1, d2, Tmax

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
    
    : Parameters tau_r and tau_d have different function for GABA-B receptor:
    :   they represents rise and decay time constants for neurotransmitter
    :   concentration that kicks of the signaling cascade.
    tau_r_GABAB = 0.2   (ms)  : Dual-exponential neurotransmitter profile
    tau_d_GABAB = 1.7   (ms)  : IMPORTANT: tau_r < tau_d

    Erev_GABAA = -80    (mV)  : GABAA reversal potential
    Erev_GABAB = -95    (mV)  : GABAB reversal potential
    
    : NOTE: gmax can be set either via NetCon.weight[0] in [nS] or here in [uS] (assuming that weight is 1)
    gmax_GABAA = .001   (uS)  : Weight conversion factor (from nS to uS)
    gmax_GABAB = .001   (uS)  : Weight conversion factor (from nS to uS)
                              : (typically about 70% of GABA-A conductance)

    : Parameters of TM model
    tau_rec = 200      (ms)  : time constant of recovery from depression
    tau_facil = 200    (ms)  : time constant of facilitation
    U1 = 0.5           (1)   : baseline release probability

    use_stdp_A = 1     (1)   : flag indicating whether STD/STP is enabled for GABA-A
    use_stdp_B = 1     (1)   : flag indicating whether STD/STP is enabled for GABA-B

: Parameters for kinetic model GABA-B signaling cascade

    K1  = 0.66  (/ms mM)    : forward binding rate to receptor
    K2  = 0.020 (/ms)       : backward (unbinding) rate of receptor
    K3  = 0.083 (/ms)       : rate of G-protein production
    K4  = 0.0079 (/ms)      : rate of G-protein decay
    d1  = 0.017 (/ms)       : rate of desensitization
    d2  = 0.0053 (/ms)      : rate of re-sensitization
    KD  = 100               : dissociation constant of K+ channel
    n   = 4                 : nb of binding sites of G-protein on K+
    Tmax = 1.3 (mM)         : peak neurotransmitter concentration
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

    Rrp          : Running fraction of available vesicles (vesicles in RRP)
                 : Accounts for synaptic depression
    
    Use          : Running value of the release probability ('utilization of synaptic efficacy')
                 : Accounts for synaptic facilitation

: State variables for GABA-B signaling cascade
    R               : fraction of activated receptor
    D               : fraction of desensitized receptor
    G               : fraction of activated G-protein
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
    Gn
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
    Rrp = 1
    Use = 0

    : GABA-B signaling cascade
    R = 0
    D = 0
    G = 0

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

    : For GABAB, use signaling cascade to determine receptor activation
    Gn = G^n
    g_GABAB = gmax_GABAB * Gn / (Gn+KD)

    : Total conductance
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
    Rrp' = (1-Rrp)/tau_rec
    Use' = -Use/tau_facil

    : Signaling cascade for GABA-B receptor activation (kinetic scheme)
    : (originally in DERIVATIVE block solved with method euler)
    R' = K1 * Tmax*(B_GABAB-A_GABAB) * (1-R-D) - K2 * R + d2 * D
    D' = d1 * R - d2 * D
    G' = K3 * R - K4 * G
}

: Block to be executed for a pre-synaptic spike event
NET_RECEIVE (weight) {
    LOCAL A, A_GA, A_GB
    
    Use = Use + U1*(1-Use) : Increment release probability.
                           : (spike-triggered delta-term in diff. eq. for Use)
    
    A = Use*Rrp : A is product of (size of RRP, i.e. available vesicles) * (their release probability).
              : This corresponds to the amount of vesicles released.
              : That's why A functions both as the decrement for R and the increment for gsyn.
    
    Rrp = Rrp - A : Decrease size of RRP (spike-triggered delta-term in diff. eq. for R)
    
    : Enable STD/STP for GABA-A
    A_GA = A
    if (use_stdp_A == 0) {
        A_GA = 1
    }
    A_GABAA = A_GABAA + A_GA*weight*factor_GABAA
    B_GABAA = B_GABAA + A_GA*weight*factor_GABAA
    
    : Enable STD/STP for GABA-B
    A_GB = A
    if (use_stdp_B == 0) {
        A_GB = 1
    }
    A_GABAB = A_GABAB + A_GB*weight*factor_GABAB
    B_GABAB = B_GABAB + A_GB*weight*factor_GABAB
}
