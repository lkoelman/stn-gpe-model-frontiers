: Ih current
: Created 8/6/02 - nwg

NEURON {
	SUFFIX hpkj
	NONSPECIFIC_CURRENT i
	RANGE ghbar, eh
	GLOBAL ninf, ntau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	v	 	(mV)
	
	ghbar = .0001	(S/cm2)

	eh = -30	(mV)
}

ASSIGNED {
	i (mA/cm2)
	ninf
	ntau
}

STATE {
	n
}

INITIAL {
	rates(v)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar*n*(v - eh)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/ntau
}

PROCEDURE rates(v (mV)) {
	ninf = 1/(1+exp((v+90.1)/9.9))
	ntau = 1000 * (.19 + .72*exp(-((v-(-81.5))/11.9)^2))
}