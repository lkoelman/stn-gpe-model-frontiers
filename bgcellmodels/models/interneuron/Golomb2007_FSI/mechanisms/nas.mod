COMMENT
Conceptual model:  Sodium current for a model of a fast-spiking cortical interneuron.

Authors and citation:
    Golomb D, Donner K, Shacham L, Shlosberg D, Amitai Y, Hansel D (2007).
    Mechanisms of Firing Patterns in Fast-Spiking Cortical Interneurons. 
    PLoS Comput Biol 3:e156.

Original implementation and programming language/simulation environment:
    by Golomb et al. for XPP
    Available from http://senselab.med.yale.edu/modeldb/ShowModel.asp?model=97747

This implementation is by N.T. Carnevale and V. Yamini for NEURON.

Revisions:
20130415 NTC introduced tiny first order delay in m 
so that simulations with fixed dt > 0.02 ms would be stable.
With taum = 0.001 ms, fixed dt simulations show slight differences 
in spike timing compared to the original results, 
but adaptive integration with cvode.atol (absolute error tolerance) 1e-4 
and proper tolerance scaling of these states
    statename   cvode.atolscale("statename")
    v           10
    m_nas       1
    a_kd        0.1
    b_kd        "
    n_kdr       "
    h_nas       "
produces results nearly identical to the original published figure.
ENDCOMMENT

NEURON {
    SUFFIX NasFSI
    USEION na READ ena WRITE ina
    RANGE gna, g
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gna = 0.1125 (S/cm2)
    theta_m = -24 (mV)
    sigma_m = 11.5 (mV)
    theta_h = -58.3 (mV)
    sigma_h = -6.7 (mV)
    theta_t_h = -60 (mV)
    sigma_t_h = -12 (mV)
    taum = 0.001 (ms) : for stability with dt>0.01 ms
}

ASSIGNED {
    v (mV)
    ena (mV)
    ina (mA/cm2)
    g (S/cm2)
    hinf
    tauh
    minf
}

STATE {
    m
    h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gna * h * m^3
    ina = g * (v-ena)
}

INITIAL {
    settables(v)
    m = minf
    h = hinf
}

DERIVATIVE states {
    settables(v)
    m' = (minf - m) / taum
    h' = (hinf - h) / tauh
}

PROCEDURE settables(v) {
    TABLE hinf, tauh, minf FROM -100 TO 100 WITH 400

    hinf = 1/(1 + exp(-(v-theta_h)/sigma_h))
    tauh = 0.5 + 14 / ( 1 + exp(-(v-theta_t_h)/sigma_t_h))
    minf = 1/(1 + exp(-(v-theta_m)/sigma_m))
}