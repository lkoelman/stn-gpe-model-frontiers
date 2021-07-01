TITLE HCN channel for GPe neuron

COMMENT
DESCRIPTION
-----------

H current  (Anomalous rectifier--mixed Na and K current)

HCN2 homomeric channel, GP-specific Channel model from 
Chan et al (2004), J Neurosci 24: 9921-32.

Original model from Siegelbaum lab. Wang et al (2002), Neuron 36:
451-62. Chen et al (2001), JGP 117: 491-504.


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
    SUFFIX HCN2
    NONSPECIFIC_CURRENT ih
    RANGE gmax, iHCN
}

PARAMETER {
    gmax  = 0.5e-4 (mho/cm2)
    eh = -30 (mV)

    : m-gate
    theta_m = -87.5 (mV)
    k_m = -4.0 (mV)
    tau_m0 = 0.0 (ms)
    tau_m1 = 25200 (ms)
    sigma_m0 = 8.2 (mV)
    sigma_m1 = 8.9 (mV)
    dQ10 = 4
}

STATE {
    m
}

ASSIGNED {
    : read simulator variables
    v (mV)
    
    : assigned mechanism variables
    ih (mA/cm2)
    iHCN (mA/cm2)

    minf
    taum (ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ih  = gmax*m*(v-eh)
    iHCN = ih
}

UNITSOFF

INITIAL {
    settables(v)
    m = minf
}

DERIVATIVE states {  
    settables(v)
    m' = (minf - m)/taum
}

PROCEDURE settables(v) {
    TABLE minf, taum FROM -100 TO 100 WITH 400

    : m-gate
    minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
    taum = (tau_m0 + tau_m1 / (exp((v - theta_m)/sigma_m0) + exp((theta_m - v)/sigma_m1))) / dQ10
}

UNITSON
