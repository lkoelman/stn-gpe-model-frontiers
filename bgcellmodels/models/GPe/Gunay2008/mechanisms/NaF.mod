TITLE fast sodium channel NaF for GPe neuron

COMMENT

DESCRIPTION
-----------

Fast Na channel

Activation and fast inactivation made to replicate resurgent
sodium current from Raman & Bean.

Slow inactivation gate added by J. Edgerton, 2004.
Support for voltage-dependent Z-gate by Cengiz Gunay, 2004


CREDITS
-------

modeled in GENESIS by Gunay et al., 2008
implemented in NEURON by Kitano, 2011
modified in NEURON by Lucas Koelman, 2018 to reflect model Hendrickson et al., 2010
ENDCOMMENT

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
}

NEURON {
 SUFFIX NaF
 USEION na READ ena WRITE ina
 RANGE gmax, iNa
}

PARAMETER {
    gmax  = 0.001 (mho/cm2)

    : m-gate - activation & deactivation
    k_m = 5.0 (mV)          : @aliases Km
    tau_m0 = 2.8e-02 (ms)   : @aliases taummin
    tau_m1 = 2.8e-02 (ms)   : @aliases taummax
    k_taum0 = 20.0 (mV)     : @aliases Ktau1
    k_taum1 = 20.0 (mV)     : @aliases Ktau2
    theta_m0 = -32.4 (mV)   : @aliases Vhalfm

    : h-gate - fast inactivaton
    theta_h = -48.0 (mV)    : @aliases V0h
    k_h = -2.8 (mV)         : @aliases Kh
    tau_h0 = 0.25 (ms)      : @aliases tauhmin
    tau_h1 = 4.0 (ms)       : @aliases tauhmax
    phi_h = -43.0 (mV)      : @aliases tauhV0
    sigma_h0 = 5.0 (mV)     : @aliases Ktau1h
    sigma_h1 = 10.0 (mV)    : @aliases Ktau2h

    : s-gate - slow inactivaton
    s0 = 0.15               : @aliases mins
    theta_s = -40.0 (mV)    : @aliases V0s
    k_s = -5.4 (mV)         : @aliases Ks
    tau_s0 = 10.0 (ms)      : @aliases tausmin
    tau_s1 = 1000.0 (ms)    : @aliases tausmax
    phi_s = -40.0 (mV)      : @aliases tausV0
    sigma_s0 = 18.3 (mV)    : @aliases Ktaus1
    sigma_s1 = 10.0 (mV)    : @aliases Ktaus2
}

STATE {
    m h s
}

ASSIGNED {
    : read built-in variables
    v (mV)
    ena (mV)

    : assigned built-in variables
    ina (mA/cm2)

    : assigned mechanism variables
    iNa (mA/cm2)

    minf
    taum (ms)

    hinf
    tauh (ms)

    sinf
    taus (ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina  = gmax*m*m*m*h*s*(v-ena)
    iNa = ina
}

UNITSOFF

INITIAL {
    settables(v)
    m = minf
    h = hinf
    s = sinf
}

DERIVATIVE states {  
    settables(v)
    m' = (minf - m)/taum
    h' = (hinf - h)/tauh
    s' = (sinf - s)/taus
}

PROCEDURE settables(v) {
    LOCAL theta_m
    TABLE minf, taum, hinf, tauh, sinf, taus FROM -100 TO 100 WITH 400

    : derived parameters cannot go in INITIAL (uninitialized when table made)
    theta_m = theta_m0 + (k_m * (log((1 / pow(0.5, 1/3)) - 1)))

    : m-gate - activation & deactivation
    minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
    taum = tau_m0 + ((tau_m1 - tau_m0) / ( exp((v - theta_m)/k_taum0) 
                                         + exp((theta_m - v)/k_taum1) ) )

    : h-gate - fast inactivaton
    hinf = 1.0 / (1.0 + exp((theta_h - v)/k_h))
    tauh = tau_h0 + ((tau_h1 - tau_h0) / ( exp((v - phi_h)/sigma_h0) 
                                         + exp((phi_h - v)/sigma_h1) ) )

    : s-gate - slow inactivaton
    sinf = s0 + (1.0 - s0)/(1.0 + exp((theta_s - v)/k_s))
    taus = tau_s0 + (tau_s1 - tau_s0)/(exp((phi_s - v)/sigma_s0) + exp((v - phi_s)/sigma_s1))
}

UNITSON
