TITLE calcium-dependent potassium (SK) channel for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
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
 v (mV)
 dt (ms)
 gmax  = 0.001 (mho/cm2)
 iSK  = 0.0 (mA/cm2)
 cai (mM)
 ek (mV)

 ECh = 0.00035 (mM)
 HCoeff = 4.6
 Ca_sat = 0.005 (mM)
 tau_m0 = 4.0 (ms)
 tau_m1 = 76.0 (ms)
}

STATE {
 m
}

ASSIGNED { 
 ik (mA/cm2)
 minf
 taum (ms)
}

BREAKPOINT {
 SOLVE states METHOD cnexp
 ik  = gmax*m*(v-ek)
 iSK = ik
}

UNITSOFF

INITIAL {
 settables(cai)
 m = minf
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
	 taum = tau_m1 - cai*(tau_m1 - tau_m0)/Ca_sat
	} else {
	 taum = tau_m0
	}
}

UNITSON
