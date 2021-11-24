TITLE persistent sodium channel NaP for GPe neuron

COMMENT

DESCRIPTION
-----------

Persistent Na channel

Based on Magistretti & Alonso (1999), JGP 114:491-509
and Magistretti & Alonso (2002), JGP 120: 855-873.

Created by J.R. Edgerton, 03/2004
Modified 10/2004 by J.R. Edgerton: add z-gate slow inactivation, improve 
model's y-gate intermediate inactivation.


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
    SUFFIX NaP
    USEION na READ ena WRITE ina
    RANGE gmax, iNa
}

PARAMETER {
    gmax  = 0.001 (mho/cm2)

    celsius (degC)
    activate_Q10 = 1
    Q10 = 3

    : Activation & Deactivation (m-gate)
    theta_m0 = -50.0 (mV)
    k_m = 5.7 (mV)
    phi_m = -41.6 (mV)
    sigma_m = 14.4 (mV)
    alp0 = 2.130 (1/ms)
    bet0 = 2.460 (1/ms)

    : Slow Inactivation (s-gate)
    h0 = 0.154
    theta_h = -57.0 (mV)
    k_h = -4.0 (mV)
    tau_h0 = 30.0 (ms)
    tau_h1 = 51.0 (ms)
    phi_h = -34.0 (mV)
    sigma_h0 = -26.0 (mV)
    sigma_h1 = -31.9 (mV)

    : Slow Inactivation (s-gate)
    theta_s = -10.0 (mV)
    k_s = -4.9 (mV)
    Aa_s = -2.88e-6 (1/ms/mV)
    Ba_s = -4.9e-5 (1/ms)
    Ka_s = 4.63 (mV)
    Ab_s = 6.94e-6 (1/ms/mV)
    Bb_s = 4.47e-4 (1/ms)
    Kb_s = -2.63 (mV)
}

STATE {
    m h s
}

ASSIGNED {
    : read simulator variables
    v (mV)
    ena (mV)

    : assigned simulator variables
    ina (mA/cm2)

    : assigned mechanism variables
    iNa (mA/cm2)

    minf
    taum (ms)
    
    hinf
    tauh (ms)

    sinf
    alphas (1/ms)
    betas (1/ms)
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
    LOCAL alpham, betam, aphas, betas, theta_m, T_Q10
    TABLE minf, taum, hinf, tauh, sinf, taus DEPEND celsius FROM -100 TO 100 WITH 400

    : Temperature adjustment for rates
    if (activate_Q10>0) {
        T_Q10 = Q10^((celsius-22)/10)
    } else {
        T_Q10 = 1.0
    }

    : derived parameters cannot go in INITIAL (uninitialized when table made)
    theta_m = theta_m0 + (k_m * (log((1 / pow(0.5, 1/3)) - 1)))

    : Activation & Deactivation (m-gate)
    minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
    alpham = alp0 * exp((v - phi_m) / sigma_m) * T_Q10
    betam = bet0 * exp((phi_m - v) / sigma_m) * T_Q10
    taum = 1.0 / (alpham + betam)

    : Fast / Intermediate Inactivation (h-gate)
    hinf = h0 + (1.0 - h0)/ (1.0 + exp((theta_h - v)/k_h))
    tauh = tau_h0 + (tau_h1 - tau_h0)/(exp((v - phi_h)/sigma_h0) + exp((phi_h - v)/sigma_h1))

    : Slow Inactivation (s-gate)
    sinf = 1.0 / (1.0 + exp((theta_s - v)/k_s))
    alphas = (Aa_s * v + Ba_s)/(1.0 - exp((v + Ba_s/Aa_s)/Ka_s)) * T_Q10
    betas = (Ab_s * v + Bb_s)/(1.0 - exp((v + Bb_s/Ab_s)/Kb_s)) * T_Q10
    taus = 1.0 / (alphas + betas)
}

UNITSON
