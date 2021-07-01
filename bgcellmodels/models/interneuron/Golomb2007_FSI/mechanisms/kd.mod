TITLE Kv1 Current (d-type)

COMMENT
Conceptual model:  D current for a model of a fast-spiking cortical interneuron.

Authors and citation:
  Golomb D, Donner K, Shacham L, Shlosberg D, Amitai Y, Hansel D (2007).
  Mechanisms of Firing Patterns in Fast-Spiking Cortical Interneurons. 
  PLoS Comput Biol 3:e156.

Original implementation and programming language/simulation environment:
  by Golomb et al. for XPP
  Available from http://senselab.med.yale.edu/modeldb/ShowModel.asp?model=97747

This implementation is by N.T. Carnevale and V. Yamini for NEURON.
ENDCOMMENT

NEURON {
    SUFFIX KdFSI
    USEION k READ ek WRITE ik
    RANGE gkd, g
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gkd = 0.00039 (S/cm2)
    theta_a = -50 (mV)
    sigma_a = 20 (mV)
    theta_b = -70 (mV)
    sigma_b = -6 (mV)
    tau_a = 2 (ms)
    tau_b = 150 (ms)
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    g (S/cm2)
    ainf
    binf
}

STATE {
    a
    b
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gkd * a^3 * b
    ik = g * (v-ek)
}

INITIAL {
    settables(v)
    a = ainf
    b = binf
}

DERIVATIVE states {
    settables(v)
    a' = (ainf-a)/tau_a
    b' = (binf-b)/tau_b
}

PROCEDURE settables(v) {
    TABLE ainf, binf FROM -100 TO 100 WITH 400

    ainf = 1/(1 + exp(-(v-theta_a)/sigma_a))
    binf = 1/(1 + exp(-(v-theta_b)/sigma_b))
}