TITLE calcium-dependent potassium (SK) channel for GPe neuron

COMMENT
DESCRIPTION
-----------

SK channel from Volker Steuber's DCN neuron model, modified to reflect Hill fits
in the following:

    - Hirschberg et al, 1999: Biophys J. 77: 1905-13. 
    - Keen et al, 1999: J. Neurosci 19: 8830-38.

Tau-Ca equation made by Volker based on Hirschberg et al, 1998: JGP 111: 565-581.


CREDITS
-------

modeled in GENESIS by J.R. Edgerton, 2004
implemented in NEURON by Kitano, 2011
modified in NEURON by Lucas Koelman, 2018 to reflect model Hendrickson et al., 2010
ENDCOMMENT

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
}

NEURON {
    SUFFIX SK
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gmax, iSK
}

PARAMETER {
    gmax  = 0.005 (mho/cm2)

    ECh = 0.00035 (mM)
    HCoeff = 4.6
    Ca_sat = 0.005 (mM)
    tau_m0 = 4.0 (ms)
    tau_m1 = 76.0 (ms)
    dQ10_SK = 2
}

ASSIGNED {
    : read simulator variables
    v (mV)
    cai (mM)
    ek (mV)

    : assigned simulator variables
    ik (mA/cm2)

    : assigned mechanism variables
    iSK (mA/cm2)

    minf
    taum (ms)
}

STATE {
    m
}

INITIAL {
    settables(cai)
    m = minf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik  = gmax*m*(v-ek)
    iSK = ik
}

DERIVATIVE states {  
    settables(cai)
    m' = (minf - m)/taum
}

PROCEDURE settables(cai) { LOCAL can, can0
    can = pow(cai*1000, HCoeff)
    can0 = pow(ECh*1000, HCoeff)
    minf = can/(can + can0)
    if(cai < Ca_sat){
        taum = (tau_m1 - cai*(tau_m1 - tau_m0)/Ca_sat) / dQ10_SK
    } else {
        taum = tau_m0 / dQ10_SK
    }
}

