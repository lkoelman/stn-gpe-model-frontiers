TITLE fast sodium channel NaF for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
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
 v (mV)
 dt (ms)
 gmax  = 0.001 (mho/cm2)
 iNa  = 0.0 (mA/cm2)
 ena (mV)

 theta_m = -39.0 (mV)
 k_m = 5.0 (mV)
 taum = 0.028 (ms)

 theta_h = -48.0 (mV)
 k_h = -2.8 (mV)
 tau_h0 = 0.25 (ms)
 tau_h1 = 4.0 (ms)
 phi_h = -43.0 (mV)
 sigma_h0 = 10.0 (mV)
 sigma_h1 = -5.0 (mV)

 s0 = 0.15
 theta_s = -40.0 (mV)
 k_s = -5.4 (mV)
 tau_s0 = 10.0 (ms)
 tau_s1 = 1000.0 (ms)
 phi_s = -40.0 (mV)
 sigma_s0 = 18.3 (mV)
 sigma_s1 = -10.0 (mV)
}

STATE {
 m h s
}

ASSIGNED { 
 ina (mA/cm2)
 minf
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
        TABLE minf, hinf, tauh, sinf, taus FROM -100 TO 100 WITH 400

	minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
	hinf = 1.0 / (1.0 + exp((theta_h - v)/k_h))
	tauh = tau_h0 + (tau_h1 - tau_h0)/(exp((phi_h - v)/sigma_h0) + exp((phi_h - v)/sigma_h1))
	sinf = s0 + (1.0 - s0)/(1.0 + exp((theta_s - v)/k_s))
	taus = tau_s0 + (tau_s1 - tau_s0)/(exp((phi_s - v)/sigma_s0) + exp((phi_s - v)/sigma_s1))
}

UNITSON
