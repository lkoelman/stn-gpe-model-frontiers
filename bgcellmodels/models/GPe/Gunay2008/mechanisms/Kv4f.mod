TITLE fast A-type potassium Kv4 channel for GPe neuron

COMMENT
DESCRIPTION
-----------

KA Kv4 fast
Low voltage-activated

Created based on GP data:
Tkatch, Baranauskas and Surmeier 2000, J Neurosci 20(2):579-588
Modified by J. R. Edgerton 02/2004


CREDITS
-------

modeled in GENESIS by J.R. Edgerton, 2004
implemented in NEURON by Kitano, 2011
modified in NEURON by Lucas Koelman, 2018 to reflect model Hendrickson et al., 2010
ENDCOMMENT

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX Kv4f
    USEION k READ ek WRITE ik
    RANGE gmax, iKv4f
}

PARAMETER {
    gmax  = 0.001 (mho/cm2)

    : m-gate
    theta_m = -49.0 (mV)
    k_m = 12.5 (mV)
    tau_m0 = 0.25 (ms)
    tau_m1 = 7.0 (ms)
    phi_m = -49.0 (mV)
    sigma_m0 = 29.0 (mV)
    sigma_m1 = 29.0 (mV)

    : h-gate
    theta_h = -83.0 (mV)
    k_h = -10.0 (mV)
    tau_h0 = 7.0 (ms)
    tau_h1 = 21.0 (ms)
    phi_h = -83.0 (mV)
    sigma_h0 = 10.0 (mV)
    sigma_h1 = 10.0 (mV)
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
    iKv4f (mA/cm2)

    minf
    taum (ms)
    hinf
    tauh (ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik  = gmax*m*m*m*m*h*(v-ek)
    iKv4f = ik
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
    TABLE minf, taum, hinf, tauh FROM -100 TO 100 WITH 400

    : m-gate (n-gate in Hendrickson et al., 2010)
    minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
    taum = tau_m0 + (tau_m1 - tau_m0)/(exp((phi_m - v)/sigma_m0) + exp((v - phi_m)/sigma_m1))

    : h-gate
    hinf = 1.0 / (1.0 + exp((theta_h - v)/k_h))
    tauh = tau_h0 + (tau_h1 - tau_h0)/(exp((phi_h - v)/sigma_h0) + exp((v - phi_h)/sigma_h1))
}

UNITSON
