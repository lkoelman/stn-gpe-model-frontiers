TITLE HCN channel for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
ENDCOMMENT

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
}

NEURON {
 SUFFIX HCN
 NONSPECIFIC_CURRENT ih
 RANGE gmax, iHCN
}

PARAMETER {
 v (mV)
 dt (ms)
 gmax  = 0.001 (mho/cm2)
 iHCN  = 0.0 (mA/cm2)
 e (mV)

 theta_m = -76.4 (mV)
 k_m = -3.3 (mV)
 tau_m0 = 0.0 (ms)
 tau_m1 = 3625.0 (ms)
 phi_m = -76.4 (mV)
 sigma_m0 = 6.56 (mV)
 sigma_m1 = -7.48 (mV)
}

STATE {
 m
}

ASSIGNED { 
 ih (mA/cm2)
 minf
 taum (ms)
}

BREAKPOINT {
 SOLVE states METHOD cnexp
 ih  = gmax*m*(v-e)
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

PROCEDURE settables(v) { LOCAL tau
        TABLE minf, taum FROM -100 TO 100 WITH 400

	minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
	tau = tau_m0 + (tau_m1 - tau_m0)/(exp((phi_m - v)/sigma_m0) + exp((phi_m - v)/sigma_m1))
	if(tau < 0.01) {
	 taum = 0.01
	} else {
	 taum = tau
	}
}

UNITSON
