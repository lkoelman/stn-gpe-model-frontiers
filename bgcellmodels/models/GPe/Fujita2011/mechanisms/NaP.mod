TITLE persistent sodium channel NaP for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
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
 v (mV)
 dt (ms)
 gmax  = 0.001 (mho/cm2)
 iNa  = 0.0 (mA/cm2)
 ena (mV)

 theta_m = -57.7 (mV)
 k_m = 5.7 (mV)
 tau_m0 = 0.03 (ms)
 tau_m1 = 0.146 (ms)
 phi_m = -42.6 (mV)
 sigma_m0 = 14.4 (mV)
 sigma_m1 = -14.4 (mV)

 h0 = 0.154
 theta_h = -57.0 (mV)
 k_h = -4.0 (mV)
 tau_h0 = 10.0 (ms)
 tau_h1 = 17.0 (ms)
 phi_h = -34.0 (mV)
 sigma_h0 = 26.0 (mV)
 sigma_h1 = -31.9 (mV)

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
 ina (mA/cm2)
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
: SOLVE states METHOD cnexp
 SOLVE states
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
        TABLE minf, taum, hinf, tauh, sinf, alphas, betas, taus FROM -100 TO 100 WITH 400

	minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
	taum = tau_m0 + (tau_m1 - tau_m0)/(exp((phi_m - v)/sigma_m0) + exp((phi_m - v)/sigma_m1))
	
    hinf = h0 + (1.0 - h0)/ (1.0 + exp((theta_h - v)/k_h))
	tauh = tau_h0 + (tau_h1 - tau_h0)/(exp((phi_h - v)/sigma_h0) + exp((phi_h - v)/sigma_h1))
	
    sinf = 1.0 / (1.0 + exp((theta_s - v)/k_s))
	alphas = (Aa_s * v + Ba_s)/(1.0 - exp((v + Ba_s/Aa_s)/Ka_s))
	betas = (Ab_s * v + Bb_s)/(1.0 - exp((v + Bb_s/Ab_s)/Kb_s))
	taus = 1.0 / (alphas + betas)
}

UNITSON
