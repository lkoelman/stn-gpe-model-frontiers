TITLE slowly activated potassium Kv2 (Kv2.1) channel for GPe neuron

COMMENT
DESCRIPTION
-----------

Kdr Kv2 (Kv2.1) slow activating

Created based on GP data:
Baranuskas, Tkatch and Surmeier
1999, J Neurosci 19(15):6394-6404


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
    SUFFIX Kv2
    USEION k READ ek WRITE ik
    RANGE gmax, iKv2
}

PARAMETER {
    gmax  = 0.001 (mho/cm2)

    : m-gate
    theta_m0 = -18.0 (mV)
    k_m = 9.1 (mV)
    tau_m0 = 0.1 (ms)
    tau_m1 = 3.0 (ms)
    sigma_m0 = -13.91 (mV)
    sigma_m1 = -21.74 (mV)

    : h-gate
    h0 = 0.2
    theta_h = -20.0 (mV)
    k_h = 10.0 (mV)
    phi_h = 0.0 (mV)
    sigma_h = 10.0 (mV)
    tau_h0 = 3400.0 (ms)
    tau_h1 = 3400.0 (ms)
}

STATE {
    m h
}

ASSIGNED {
    : read simulator variables
    v (mV)
    ek (mV)

    : assigned simulator variables
    ik (mA/cm2)

    : assigned mechanism variables
    iKv2 (mA/cm2)

    minf
    taum (ms)

    hinf
    tauh (ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik  = gmax*m*m*m*m*h*(v-ek)
    iKv2 = ik
}

UNITSOFF

INITIAL {
    settables(v)

    m = minf
    h = hinf
}

DERIVATIVE states {  
    settables(v)
    m' = (minf - m)/taum
    h' = (hinf - h)/tauh
}

PROCEDURE settables(v) {
    LOCAL theta_m
    TABLE minf, taum, hinf, tauh FROM -100 TO 100 WITH 400

    : derived parameters cannot go in INITIAL (uninitialized when table made)
    theta_m = theta_m0 + (k_m * (log((1 / pow(0.5, 1/4)) - 1)))

    : m-gate (also called n-gate)
    minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
    taum = tau_m0 + (tau_m1 - tau_m0)/(exp((theta_m - v)/sigma_m0) + exp(-(theta_m - v)/sigma_m1))

    : h-gate
    hinf = h0 + (1.0 - h0) / (1.0 + exp((v - theta_h)/k_h))
    tauh = (tau_h1-tau_h0) / (1 + exp((v - phi_h) / sigma_h ))  + tau_h0
}

UNITSON
