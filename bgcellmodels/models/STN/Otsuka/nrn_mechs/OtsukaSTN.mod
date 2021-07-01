TITLE Ion channels for Otsuka STN neuron

COMMENT
Next to the HH currents 
  - I_Na - fast sodium current
  - I_K - delayed rectifier K+ current
  - I_l - leak current
, the Otsuka STN neuron contains four additional currents:
  - I_A - A-type K+ current with low activation threshold and fast activation and inactivation time constants
  - I_L - L-type Ca2+ currents, long-lasting, involved in generation of Plateau potentials
  - I_T - T-type Ca2+ currents, low threshold
  - I_Ca-K - Ca2+ activated K+ channels (repolarization component of plateu potentials)

Disclaimer: this is a modified version of the mechanism found in the Hahn & McIntyre model
published in ModelDB at http://senselab.med.yale.edu/modeldb/showModel.cshtml?model=127388&file=\BGnet\pSTN.tem

ENDCOMMENT

: Public interface of the mechanism
NEURON {
    SUFFIX stn

    : All ions and derived ionic variables used in computations
    USEION ca READ cai, cao WRITE ica, cai, eca
    USEION k READ ek WRITE ik
    USEION na READ ena WRITE ina
    : See Otsuka2004, p. 257 (left): We need to compute [Ca]_i because the steady-state
    : value of two of the gating variables depends on this value. We then compute eca
    : ourselves via the Nernst equilibrium potential equation


    :############ Nonspecific (no ion) ############
    NONSPECIFIC_CURRENT ilk
    NONSPECIFIC_CURRENT ibias
    RANGE ilk, gl, el
    RANGE ibias, bias_amp

    :############ Na ion ############
    :Na fast (fast sodium current)
    RANGE ina, ena0
    RANGE gnabar

    :############ K ion ############
    RANGE ik, ek0
    :K delayed rectifier
    RANGE ikDR
    RANGE gkdrbar
    :Ca-activated K current
    RANGE ikCA
    RANGE gkcabar
    :A-type K current
    RANGE ikA
    RANGE gkabar

    :############ Ca ion ############
    RANGE ica
    : L-type Ca current
    RANGE icaL
    RANGE gcalbar
    : T-type Ca current
    RANGE icaT
    RANGE gcatbar
    : Ca ion kinetics equation (p. 257, left)
    : vol is a volume needed to make units of Ca kinetics equation match
    RANGE kca, alphaca

}

: Aliases for units
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (molar) = (1/liter)
    (mM)    = (millimolar)
    FARADAY = (faraday) (coulomb)  :units are really coulombs/mole
}

: Parameters (self-defined or NEURON globals)
PARAMETER {
    R = 8.31441 (Gas constant)
    T       (Absolute temp)
    celsius     (degC)
    bias_amp = 0 (mA/cm2)

:Reversal potentials (fixed)
    ena0 = 60   (mV)
    ek0 = -90   (mV)
    el  = -60   (mV)

:Leakage current
    gl  = 0.35   (S/cm2)

:Ca dynamics
    kca   = 2 (1/ms) : since NEURON uses ms, don't use value of 2e-3 (units 1/s) from paper
    alphaca = 0.000005182136


: ############ Fast Na channel ############
    gnabar   = 49 (S/cm2)
: gating variable m (eq. dm/dt has params m_inf and tau_m of sigmoidal form)
    theta_m = -40 (mV)      : eq. for m_inf
    k_m = -8 (mV)           : eq. for m_inf
    tau_m0 = 0.2 (ms)       : eq. for tau_m
    tau_m1 = 3 (ms)         : eq. for tau_m
    tht_m = -53 (mV)        : eq. for tau_m
    sig_m = -0.7 (mV)       : eq. for tau_m
: gating variable h (param tau_h in dh/dt has equation of bell-shaped form)
    theta_h = -45.5 (mV)    : eq. for h_inf
    k_h = 6.4 (mV)          : eq. for h_inf
    tau_h0 = 0 (ms)         : eq. for tau_h
    tau_h1 = 24.5 (ms)      : eq. for tau_h
    sig_h1 = -15 (mV)       : eq. for tau_h
    sig_h2 = 16 (mV)        : eq. for tau_h
    tht_h1 = -50 (mV)       : eq. for tau_h
    tht_h2 = -50 (mV)       : eq. for tau_h


: ############ Delayed rectifier K ############
    gkdrbar  = 57   (S/cm2)  
: gating variable m (dm/dt has params n_inf and tau_n of bell-shaped form)
    theta_n = -41 (mV)      : eq. for n_inf
    k_n = -14 (mV)          : eq. for n_inf
    tau_n0 = 0 (ms)         : eq. for tau_n
    tau_n1 = 11 (ms)        : eq. for tau_n
    tht_n1 = -40 (mV)       : eq. for tau_n
    tht_n2 = -40 (mV)       : eq. for tau_n
    sig_n1 = -40 (mV)       : eq. for tau_n
    sig_n2 = 50 (mV)        : eq. for tau_n


: ############ T-type Ca current ############
    gcatbar   = 5 (S/cm2)
: gating variable p (dp/dt has params q_inf and tau_p of bell-shaped form)
    theta_p = -56 (mV)
    k_p = -6.7 (mV)    
    tau_p0 = 5 (ms)
    tau_p1 = 0.33 (ms)
    tht_p1 = -27 (mV)
    tht_p2 = -102 (mV)
    sig_p1 = -10 (mV)
    sig_p2 = 15 (mV)
: gating variable q (dq/dt has params q_inf and tau_q of bell-shaped form)
    theta_q = -85 (mV)
    k_q = 5.8 (mV)
    tau_q0 = 0 (ms)
    tau_q1 = 400 (ms)
    tht_q1 = -50 (mV)
    tht_q2 = -50 (mV)
    sig_q1 = -15 (mV)
    sig_q2 = 16 (mV)


: ############ L-type Ca current ############
    gcalbar   = 15 (S/cm2)
: gating variable c (dc/dt has params c_inf and tau_c of bell-shaped form)
    theta_c = -30.6 (mV)
    k_c = -5 (mV)
    tau_c0 = 45 (ms)
    tau_c1 = 10 (ms)
    tht_c1 = -27 (mV)
    tht_c2 = -50 (mV)
    sig_c1 = -20 (mV)
    sig_c2 = 15 (mV)
: gating variable d1 (d1/dt has params d1_inf and tau_d1 of bell-shaped form)
    theta_d1 = -60 (mV)
    k_d1 = 7.5 (mV)
    tau_d10 = 400 (ms)
    tau_d11 = 500 (ms)
    tht_d11 = -40 (mV)
    tht_d12 = -20 (mV)
    sig_d11 = -15 (mV)
    sig_d12 = 20 (mV)
: gating variable d2 (d2/dt has params d2_inf and tau_d2)
:   - d2_inf is of same form as others but [Ca]_i dependent instead of voltage-dependent
:   - tau_d2 is fixed (not voltage dependent)
    theta_d2 = 0.1e-3 (mM)
    k_d2 = 0.02e-3 (mM)
    tau_d2 = 130 (ms)

: ############ A-type K current ############
    gkabar  = 5  (S/cm2)
: gating variable a (da/dt has params a_inf and tau_a of sigmoidal form)  
    theta_a = -45 (mV)
    k_a = -14.7 (mV)    
    k_b = 7.5 (mV)   
    tau_a0 = 1 (ms)
    tau_a1 = 1 (ms)
    tht_a = -40 (mV)
    sig_a = -0.5 (mV)
: gating variable b (db/dt has params b_inf and tau_b of bell-shaped form)
    theta_b = -90 (mV)
    tau_b0 = 0 (ms)
    tau_b1 = 200 (ms)
    tht_b1 = -60 (mV)
    tht_b2 = -40 (mV)
    sig_b1 = -30 (mV)
    sig_b2 = 10 (mV)

: ############ AHP current (Ca dependent K current) ############
    gkcabar   = 1 (S/cm2)
: gating variable 4 (dr/dt has fixed tau_r and [Ca]_i-dependent r_inf)
    theta_r = 0.17e-3 (mM)
    k_r = -0.08e-3 (mM)
    tau_r = 2 (ms)
    power_r = 2         : exponent of r in the equation for ikCA
}

: variables calculated by the mechanism or by Neuron and used by the mechanism
ASSIGNED {
: total membrane potential
    v   (mV)
: total ionic currents
    ina (mA/cm2)
    ik  (mA/cm2) 
    ica (mA/cm2)
: non-specific currents
    ilk (mA/cm2)
    ibias (mA/cm2)
: specific currents
    ikDR (mA/cm2)   
    ikA (mA/cm2) 
    ikCA   (mA/cm2)  
    icaT    (mA/cm2) 
    icaL    (mA/cm2)
: ion concentrations read
    cao
    
: ############ Computed parameters of ionic currents ############
: Fast Na current
    h_inf
    tau_h   (ms)
    m_inf
    tau_m   (ms)
    ena           (mV)   := 60  

: Delayed rectifier K current
    n_inf
    tau_n   (ms)
    ek         (mV) := -90

: T-type Ca current
    p_inf
    q_inf
    tau_p   (ms)
    tau_q   (ms)
    eca           (mV)   :calc from Nernst

: L_type Ca current
    c_inf
    tau_c   (ms)
    d1_inf
    tau_d1  (ms)
    d2_inf
    :tau_d2 (ms)  :in PARAMETERS

: A-type K current
    a_inf
    tau_a   (ms)
    b_inf
    tau_b   (ms)

: Ca-activated K current (AHP current)
    r_inf 
}

: State variables for kinetic schemes
STATE {
    : Gating variables
    m h     : Fast Na current
    n       : Delayed-rectifier K current
    a b     : A-type K current
    c d1 d2 : L-type Ca current
    p q     : T-type Ca current
    r       : Ca-activated K current

    : Ionic concentrations
    cai (mM) <1e-10>
}

: Main routine for the mechanism
BREAKPOINT {
    : Solve the differential equations
    SOLVE states METHOD cnexp

    : Compute Nernst potential only for Ca (rest is fixed)
    : T = 273 + celsius - 9.5
    T = celsius + 273.15
    : eca = 13.05*log(2/cai)
    eca = -1e3*(R*T)/(2*FARADAY)*log(cai/cao) : paper: cao=2, T=29.7+273.15
    
    : Na currents
    ina   = gnabar * m*m*m*h * (v - ena0)
    
    : K currents
    ikDR = gkdrbar * n*n*n*n * (v - ek0)
    ikA = gkabar * a*a*b * (v - ek0)
    ikCA = gkcabar * (v - ek0)*r^(power_r)
    ik = ikDR+ikA+ikCA : total current
    
    : Ca currents
    icaT   = gcatbar * p*p*q * (v - eca)
    icaL   = gcalbar * c*c*d1*d2 * (v - eca)
    ica = icaT+icaL : total current
    
    : Nonspecific currents
    ilk = gl * (v - el)
    ibias = -bias_amp : negative because electrode current
}

: Differential equations to update state variables
DERIVATIVE states {
    : compute instantaneous rate parameters
    set_rates(v,cai)

    : voltage-dependent gating/activation variables
    h' = (h_inf - h)/tau_h
    m' = (m_inf - m)/tau_m
    n' = (n_inf - n)/tau_n
    p' = (p_inf - p)/tau_p
    q' = (q_inf - q)/tau_q
    a' = (a_inf - a)/tau_a
    b' = (b_inf - b)/tau_b
    c' = (c_inf - c)/tau_c
    d1' = (d1_inf - d1)/tau_d1

    : Ca-depedent gating/activation variables
    d2' = (d2_inf - d2)/tau_d2
    r' = (r_inf - r)/tau_r

    : equation for point process (disregarding scaling with cell volume)
    cai' = -alphaca*ica - kca*cai

}

: Initialization of any state variables
INITIAL {
    set_rates(v,cai)

    m = m_inf 
    h = h_inf   
    n = n_inf   
    p = p_inf 
    q = q_inf  

    c = c_inf 
    d1 = d1_inf  
    d2 = d2_inf   

    a = a_inf 
    b = b_inf   

    r = r_inf 
}

:UNITSOFF

FUNCTION taux(v, tau_x0, tau_x1, theta_taux, sigma_taux) {
    : compute voltage-dependent time constant
    taux = tau_x0 + tau_x1/(1.0+exp(-(v-theta_taux)/sigma_taux))
}

FUNCTION tauxbi(v, tau_x0, tau_x1, tht_x1, tht_x2, sig_x1, sig_x2) {
    : compute voltage-dependent time constant
    tauxbi = tau_x0 + tau_x1/(exp(-(v-tht_x1)/sig_x1) + exp(-(v-tht_x2)/sig_x2))
}

FUNCTION xinf(v, theta_x, sigma_x) {
    : compute voltage-depedent steady state rate constant
    xinf = 1.0/(1.0+exp(-(v-theta_x)/sigma_x))
}

PROCEDURE set_rates(v(mV),cai(mM)) {
:Fast Na current
    h_inf = 1/(1+exp((v-theta_h)/k_h))
    m_inf = 1/(1+exp((v-theta_m)/k_m))
    tau_h = tau_h0 + tau_h1/(exp(-(v-tht_h1)/sig_h1) + exp(-(v-tht_h2)/sig_h2)) 
    tau_m = tau_m0 + tau_m1/(1+exp(-(v-tht_m)/sig_m)) 

:Delayed rectifier K
    n_inf = 1/(1+exp((v-theta_n)/k_n))
    tau_n = tau_n0 + tau_n1/(exp(-(v-tht_n1)/sig_n1) + exp(-(v-tht_n2)/sig_n2)) 

:Ca T current
    p_inf = 1/(1+exp((v-theta_p)/k_p))
    q_inf = 1/(1+exp((v-theta_q)/k_q))
    tau_p = tau_p0 + tau_p1/(exp(-(v-tht_p1)/sig_p1) + exp(-(v-tht_p2)/sig_p2)) 
    tau_q = tau_q0 + tau_q1/(exp(-(v-tht_q1)/sig_q1) + exp(-(v-tht_q2)/sig_q2))

:Ca L current
    c_inf = 1/(1+exp((v-theta_c)/k_c))
    d1_inf = 1/(1+exp((v-theta_d1)/k_d1))
    tau_c = tau_c0 + tau_c1/(exp(-(v-tht_c1)/sig_c1) + exp(-(v-tht_c2)/sig_c2))  
    tau_d1 = tau_d10 + tau_d11/(exp(-(v-tht_d11)/sig_d11) + exp(-(v-tht_d12)/sig_d12))  

:A current
    a_inf = 1/(1+exp((v-theta_a)/k_a))
    b_inf = 1/(1+exp((v-theta_b)/k_b))
    tau_a = tau_a0 + tau_a1/(1+exp(-(v-tht_a)/sig_a))
    tau_b = tau_b0 + tau_b1/(exp(-(v-tht_b1)/sig_b1) + exp(-(v-tht_b2)/sig_b2))

:Ca L current
    d2_inf = 1/(1+exp((cai-theta_d2)/k_d2))

:AHP current
    r_inf = 1/(1+exp((cai-theta_r)/k_r))
}

COMMENT
: Procedure for calculating v-dependent kinetic parameters
: This makes use of look-up table that is only set once for speed gains
: See:
: http://www.anc.ed.ac.uk/school/neuron/tutorial/tutD.html
: http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/mechanisms/nmodl.html?#table
: Note: the flag `usetable_suffix = 1/0` can be used to enable/disable the use of look-up tables
PROCEDURE set_rates_vdep(v(mV)) { 

:Make look-up table
: TODO: look at v range
: TABLE h_inf, m_inf, tau_h, tau_m, n_inf, tau_n, p_inf, q_inf, tau_p, tau_q, c_inf, d1_inf, tau_c, tau_d1, a_inf, b_inf, tau_a, tau_b FROM -100 TO 100 WITH 400

:Fast Na current
    h_inf = 1/(1+exp((v-theta_h)/k_h))
    m_inf = 1/(1+exp((v-theta_m)/k_m))
    tau_h = tau_h0 + tau_h1/(exp(-(v-tht_h1)/sig_h1) + exp(-(v-tht_h2)/sig_h2)) 
    tau_m = tau_m0 + tau_m1/(1+exp(-(v-tht_m)/sig_m)) 

:Delayed rectifier K
    n_inf = 1/(1+exp((v-theta_n)/k_n))
    tau_n = tau_n0 + tau_n1/(exp(-(v-tht_n1)/sig_n1) + exp(-(v-tht_n2)/sig_n2)) 

:Ca T current
    p_inf = 1/(1+exp((v-theta_p)/k_p))
    q_inf = 1/(1+exp((v-theta_q)/k_q))
    tau_p = tau_p0 + tau_p1/(exp(-(v-tht_p1)/sig_p1) + exp(-(v-tht_p2)/sig_p2)) 
    tau_q = tau_q0 + tau_q1/(exp(-(v-tht_q1)/sig_q1) + exp(-(v-tht_q2)/sig_q2))

:Ca L current
    c_inf = 1/(1+exp((v-theta_c)/k_c))
    d1_inf = 1/(1+exp((v-theta_d1)/k_d1))
    tau_c = tau_c0 + tau_c1/(exp(-(v-tht_c1)/sig_c1) + exp(-(v-tht_c2)/sig_c2))  
    tau_d1 = tau_d10 + tau_d11/(exp(-(v-tht_d11)/sig_d11) + exp(-(v-tht_d12)/sig_d12))  

:A current
    a_inf = 1/(1+exp((v-theta_a)/k_a))
    b_inf = 1/(1+exp((v-theta_b)/k_b))
    tau_a = tau_a0 + tau_a1/(1+exp(-(v-tht_a)/sig_a))
    tau_b = tau_b0 + tau_b1/(exp(-(v-tht_b1)/sig_b1) + exp(-(v-tht_b2)/sig_b2))  

}

: Procedure for calculating cai-dependent kinetic parameters
PROCEDURE set_rates_icadep(cai(mM)) {

:Make look-up table
: TODO: look at cai range for table
: TABLE d2_inf, r_inf FROM -100 TO 100 WITH 400

:Ca L current
    d2_inf = 1/(1+exp((cai-theta_d2)/k_d2))

:AHP current
    r_inf = 1/(1+exp((cai-theta_r)/k_r))
}
ENDCOMMENT

:UNITSON