TITLE calcium accumulation for STh

COMMENT 

 Calcium accumulation into a volume of area*depth next to the
 membrane with an exponential decay (time constant tau) to resting
 level (given by the global calcium variable cai0_ca_ion).

 How the q10 works: There is a q10 for the rates (alpha and beta's)
 called Q10 and a Q10 for the maximum conductance called gmaxQ10.  The
 q10s should have been measured at specific temperatures temp1 and
 temp2 (that are 10degC apart). Ideally, as Q10 is temperature
 dependant, we should know these two temperatures.  We used to
 follow the more formal Arrhenius derived Q10 approach.  The
 temperature at which this channel's kinetics were recorded is tempb
 (base temperature).  What we then need to calculate is the desired
 rate scale for now working at temperature celsius (rate_k).  This was
 given by the empirical Arrhenius equation, using the Q10, but now is 
 using the quick Q10 approximation. 
ENDCOMMENT

NEURON {
	SUFFIX Cacum
	USEION ca READ ica WRITE cai
	GLOBAL con,cai0,buftau,activate_Q10,Q10,rate_k,temp1,temp2,tempb,depth
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)	: Faradays constant 
}

PARAMETER {
        v (mV)
	dt (ms)
	con   = 0.0			: conversion constant (see INITIAL block)
        Avo   = 6.02e23			: Avogadro's number
	elc   = 1.602e-19 (coulombs)	: elementrary charge
	depth = 200.0 (nm)		: assume volume = area*depth
	cai0  = 0.0001(mM)		: replace cai0_ca_ion 
	buftau = 1.857456645e+02 (ms)
	cai0_ca_ion
	celsius

	activate_Q10 = 1
	Q10 = 1.2
	temp1 = 19.0 (degC)
	temp2 = 29.0 (degC)
	tempb = 23.0 (degC)
}

ASSIGNED {
	ica (mA/cm2)
        tau (ms)
	rate_k
}

STATE {
	cai (mM)
}

BREAKPOINT {
	SOLVE integrate METHOD cnexp
}

UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  rate_k = Q10^((celsius-tempb)/10)
	}else{
	  rate_k = 1.0
	}

	con=1e7/(depth*2.0*Avo*elc)

	tau=buftau/rate_k
	cai=cai0
}

DERIVATIVE integrate {
	cai' = -ica*con + (cai0 - cai)/tau
}

UNITSON
