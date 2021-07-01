TITLE GABAA and GABAB receptor with TM short-term dynamics

COMMENT

Brief Description:

GABAA and receptor conductance using a dual-exponential profile
and incorporating Tsodyks-Markram short-term dynamics. GABA-B is implemented
based on a signaling cascade proposed by Destexhe et al with the bound receptor
fraction represented by the the Tsodyks-Markram conductance variable.

Signaling cascade:

For the GABA-B signaling cascade, the Tsodyks-Markram conductance variable 
is interpreted as the bound receptor concentration, including effects of 
desensitization (corresponding to state variable 'R' in Destexhe's GABA-B 
signaling cascade model, see Destexhe, Alain, Zachary F. Mainen, and Terrence J.
Sejnowski. "Kinetic models of synaptic transmission." Methods in neuronal 
modeling 2 (1998): 1-25). The G-protein production and decay are modeled in the
same way, and the the final GABA-B conductance is a sigmoid representing Hill
kinetics applied to the G-protein concentration. The parameters of the sigmoid
and G-protein production and decay rates have been adjusted.


CREDITS
-------

- Tsodyks-Markram model based on `TsodyksMarkram_AMPA_NMDA.mod` by E. Muller at BBP/EPFL

- GABA-B signaling cascade model based on `gabab.mod` by A. Destexhe
  obtainable at http://cns.iaf.cnrs-gif.fr/alain_demos.html

- Integration of GABA-B signaling cascade in Tsodyks-Markram synapse + 
  addition of sigmoidal threshold on bound recptor fraction by Lucas Koelman


USAGE
-----


1. Find the peak value of the G-protein concentration for the maximum pre-synaptic
  firing rate

    - record state variable 'G' and plot it

    - this will depend on U1, tau_rec, tau_facil


2. Set KD (half-maximum of sigmoid describing Hill kinetics) appropriately
    
    - Configure sigmoid to map desired range of G-proten concentration to range [0, 1]

    - It makes most sense to map the actual dynamic range of G (during the simulation)
      to the linear (middle) and sublinear part of the sigmoid. E.g. by setting KD
      to slightly above peak G.
    
    - Depending on how many spikes you want to result in a significant rise in g_GABAB


3. Adjust ratio of gmax for GABA-A and GABA-B

    - the actual ratio of g during the simulation will generally differ from the
      values set as parameters due to the different multiplication factors


ENDCOMMENT

: Definition of variables which may be accessed by the user 
NEURON {
    THREADSAFE

    POINT_PROCESS GABAsyn2
    
    RANGE tau_r_GABAA, tau_d_GABAA
    RANGE tau_r_GABAB, tau_d_GABAB
    RANGE gmax_GABAA, gmax_GABAB
    RANGE i, i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, Erev_GABAA, Erev_GABAB
    RANGE tau_rec, tau_facil, U1

    RANGE G : state variables for GABAB-B dynamics
    RANGE K3, K4, KD, n : parameters for GABA-B dynamics

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
    tau_r_GABAB = 5    (ms)  : Dual-exponential neurotransmitter profile
    tau_d_GABAB = 25   (ms)  : IMPORTANT: tau_r < tau_d

    Erev_GABAA = -80    (mV)  : GABAA reversal potential
    Erev_GABAB = -95    (mV)  : GABAB reversal potential
    
    : NOTE: gmax can be set either via NetCon.weight[0] in [nS] or here in [uS] (assuming that weight is 1)
    gmax_GABAA = .001   (uS)  : Weight conversion factor (from nS to uS)
    gmax_GABAB = .001   (uS)  : Weight conversion factor (from nS to uS)
                              : (typically about 70% of GABA-A conductance)

    : Parameters of TM model
    : Good settings seem (tau_r, tau_f) = (5, 100) and (100, 100)
    tau_rec = 5        (ms)  : time constant of recovery from depression
    tau_facil = 100    (ms)  : time constant of facilitation
    U1 = 0.05           (1)   : baseline release probability

: Parameters for kinetic model GABA-B signaling cascade
    : K1  = 0.52  (/ms mM)    : forward binding rate to receptor
    : K2  = 0.0013 (/ms)      : backward (unbinding) rate of receptor
    K3  = 0.098 (/ms)       : rate of G-protein production
    K4  = 0.033 (/ms)       : rate of G-protein decay
    KD  = 100               : dissociation constant of K+ channel
    n   = 4                 : nb of binding sites of G-protein on K+
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

    Gn  : exponentiated G-protein concentration
    Kn  : exponentiated dissociation constant
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
    G = 0
    Kn = KD^n

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

    : For GABAB, use signaling cascade to determine receptor activation.
    : (Hill kinetics function is sigmoid constrained between 0-1
    :  with half-maximum at G = Kd)
    Gn = G^n
    g_GABAB = gmax_GABAB * Gn / (Gn+Kn)

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

    : ==========================================================================
    : Tsodyks-Markram dynamics

    : State variables for bi-exponential profile
    A_GABAA' = -A_GABAA/tau_r_GABAA
    B_GABAA' = -B_GABAA/tau_d_GABAA
    A_GABAB' = -A_GABAB/tau_r_GABAB
    B_GABAB' = -B_GABAB/tau_d_GABAB

    : State variables for synaptic facilitation & depression
    Rrp' = (1-Rrp)/tau_rec
    Use' = -Use/tau_facil

    : ==========================================================================
    : GABA-B signaling cascade

    : Consider Tsodyks-Markram synapse activation (B-A) as bound receptor
    : fraction including desensitization

    : G-protein production rate is proportional to bound receptor fraction.
    G' = K3 * (B_GABAB-A_GABAB) - K4 * G
}


: Block to be executed for a pre-synaptic spike event
NET_RECEIVE (weight) {
    LOCAL A

    : Use = 'running release probability' = facilitation variable
    :     => The longer (slower) tau_facil, the more easily Use builds up.
    : U1  = baseline release probability
    :     = increment in Use per spike as fraction of what's left of Use.
    :     => The higher U1, the faster Use build up (increment per spike)
    
    : Increment release probability.
    : (spike-triggered delta-term in diff. eq. for Use)
    Use = Use + U1*(1-Use) 
    
    : A is product of (size of RRP, i.e. available vesicles) * (their release probability).
    : This corresponds to the amount of vesicles released.
    : That's why A functions both as the decrement for R and the increment for gsyn.
    A = Use*Rrp
  
    : Decrease size of RRP (spike-triggered delta-term in diff. eq. for R)
    Rrp = Rrp - A
    
    A_GABAA = A_GABAA + A*weight*factor_GABAA
    B_GABAA = B_GABAA + A*weight*factor_GABAA
    A_GABAB = A_GABAB + A*weight*factor_GABAB
    B_GABAB = B_GABAB + A*weight*factor_GABAB
}
