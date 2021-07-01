TITLE high-voltage-activated calcium channel for GPe neuron

COMMENT
DESCRIPTION
-----------

HVA Ca2+ Channels 

Voltage-dependent activation from GP data:  
Surmeier Seno and Kitai 1994, J Neurophysio. 71(3):1272-1280

USAGE
-----

First insert Ca buffering mechanism into compartments where you wish to insert
this mechanism.

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
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

NEURON {
    SUFFIX CaHVA
    USEION ca READ cai,cao,eca WRITE ica
    RANGE gmax, iCaH
    : RANGE theta_m, k_m, tau_m, k_tau
}

PARAMETER {
    gmax  = 0.001 (mho/cm2)

    : m-gate
    theta_m = 20.0 (mV)
    k_m = -7.0 (mV)

    tau_m = 0.2 (ms)
    k_tau = 1e9 (mV)

    : read simulator variables
    v (mV)
    cai (mM)
    cao (mM)
    eca (mV)
    celsius (degC)
}

STATE {
    m
}

ASSIGNED { 
    : assigned simulator variables
    ica (mA/cm2)

    : assigned mechanism variables
    iCaH (mA/cm2)

    taum (ms)
    minf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica  = gmax * m * ghk(v, cai, cao)
    iCaH = ica
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
    minf = 1.0 / (1.0 + exp((v + theta_m) / k_m))
    taum = (tau_m + tau_m*1e-6 * v) / (exp((v) / k_tau))
}

: GHK equation, taken from NEURON examples (nrn/examples/nrniv/nmodl/cachan.mod)
FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    eco = co*efun(z)
    eci = ci*efun(-z)
    :high cao charge moves inward
    :negative potential charge moves inward
    ghk = (.001)*2*FARADAY*(eci - eco)
}

: prevent divide by zero, taken from NEURON examples (nrn/examples/nrniv/nmodl/cachan.mod)
FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
        efun = 1 - z/2
    }else{
        efun = z/(exp(z) - 1)
    }
}

UNITSON
