TITLE slowly activated potassium Kv2 (Kv2.1) channel for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
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
 v (mV)
 dt (ms)
 gmax  = 0.001 (mho/cm2)
 iKv2  = 0.0 (mA/cm2)
 ek (mV)

 theta_m = -33.2 (mV)
 k_m = 9.1 (mV)
 tau_m0 = 0.1 (ms)
 tau_m1 = 3.0 (ms)
 phi_m = -33.2 (mV)
 sigma_m0 = 21.7 (mV)
 sigma_m1 = -13.9 (mV)

 h0 = 0.2
 theta_h = -20.0 (mV)
 k_h = -10.0 (mV)
 tauh = 3400.0 (ms)
}

STATE {
 m h
}

ASSIGNED { 
 ik (mA/cm2)
 minf
 taum (ms)
 hinf
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
        TABLE minf, taum, hinf FROM -100 TO 100 WITH 400

	minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
	taum = tau_m0 + (tau_m1 - tau_m0)/(exp((phi_m - v)/sigma_m0) + exp((phi_m - v)/sigma_m1))
	hinf = h0 + (1.0 - h0) / (1.0 + exp((theta_h - v)/k_h))
}

UNITSON
