COMMENT
/**
 * @file AMPAsynTM.mod
 * @brief An AMPA glutamate receptor model with short-term depression
 * and facilitation
 * @author emuller, lkoelman
 * @date 2017-05-11
 * @remark Copyright © BBP/EPFL 2005-2014; All rights reserved. 
 */
ENDCOMMENT


TITLE AMPA glutamate receptor with TM short-term dynamics

COMMENT
AMPA glutamate receptor conductance using a dual-exponential profile
and incorporating Tsodyks-Markram short-term dynamics
ENDCOMMENT

: Definition of variables which may be accesses by the user 
NEURON {
    THREADSAFE

    POINT_PROCESS AMPAsynTM
    RANGE tau_r_AMPA, tau_d_AMPA
    RANGE gmax_AMPA
    RANGE tau_rec, tau_facil, U1
    RANGE i, i_AMPA, g_AMPA, g, e
    NONSPECIFIC_CURRENT i
}


: Definition of constants which may be set by the user
PARAMETER {

    tau_r_AMPA = 0.2   (ms)  : Dual-exponential conductance profile
    tau_d_AMPA = 1.7   (ms)  : IMPORTANT: tau_r < tau_d

    e = 0              (mV)  : AMPA reversal potential
    
    : NOTE: gmax can be set either via NetCon.weight[0] in [nS] or here in [uS] (assuming that weight is 1)
    gmax_AMPA = .001   (uS)  : Weight conversion factor (from nS to uS)

    : Parameters of TM model
    tau_rec = 200      (ms)  : time constant of recovery from depression
    tau_facil = 200    (ms)  : time constant of facilitation
    U1 = 0.5           (1)   : baseline release probability
}

: Declaration of state variables 
STATE {

    A_AMPA       : AMPA state variable to construct the dual-exponential profile
                 : rise kinetics with tau_r_AMPA

    B_AMPA       : AMPA state variable to construct the dual-exponential profile
                 : decay kinetics with tau_d_AMPA

: State variables for the TM model
    R            : Running fraction of available vesicles
    Use          : Running value of the release probability
}

: Declaration of variables that are computed, e.g. in the BREAKPOINT block
ASSIGNED {
    v (mV)
    i (nA)
    i_AMPA (nA)
    g_AMPA (uS)
    g (uS)
    factor_AMPA
}

: Definition of initial conditions
INITIAL{
    LOCAL tp_AMPA     : Declaration of some local variables

    : Zero receptor rise and fall kinetics variables
    A_AMPA = 0
    B_AMPA = 0

    : Compute constants needed to normalize the dual-exponential receptor dynamics

    : Expression for time to peak of the AMPA dual-exponential conductance
    tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA)

    : AMPA Normalization factor - so that when t = tp_AMPA, gsyn = gpeak
    factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA) 
    factor_AMPA = 1/factor_AMPA

    R = 1
    Use = 0

}

: Declare method to propagate the state variables in time
BREAKPOINT {

    : Specify to solve system of equations "odes", declared below (DERIVATIVE block)
    : "cnexp" specifies the intergration method, it is
    : an implicit integration method that is stable even for stiff systems
    SOLVE odes METHOD cnexp

    : Compute and assign quantities which depend on the state variables

    : Compute the time varying AMPA receptor conductance as 
    : the difference of state variables B_AMPA and A_AMPA
    g_AMPA = gmax_AMPA*(B_AMPA-A_AMPA) 

    : Total conductance
    g = g_AMPA

    i_AMPA = g_AMPA*(v-e) 

    : Compute the total current
    i = i_AMPA
}

: Declaration of ODEs solved for in the BREAKPOINT block
DERIVATIVE odes {
    A_AMPA' = -A_AMPA/tau_r_AMPA
    B_AMPA' = -B_AMPA/tau_d_AMPA

    R' = (1-R)/tau_rec
    Use' = -Use/tau_facil
}

: Block to be executed for a pre-synaptic spike event
NET_RECEIVE (weight) {
    LOCAL A
    
    Use = Use + U1*(1-Use) : spike-triggered (delta) term in diff. eq. for Use
    A = Use*R
    R = R - A : spike-triggered (delta) term in diff. eq. for R
    
    A_AMPA = A_AMPA + A*weight*factor_AMPA
    B_AMPA = B_AMPA + A*weight*factor_AMPA
}

