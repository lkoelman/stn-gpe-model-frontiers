: HH Low TEA-sensitive Purkinje potassium current
: Created 8/7/02 - nwg

NEURON {
	SUFFIX kpkj2
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	GLOBAL ninf, ntau
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	v		(mV)
	gkbar = .002	(mho/cm2)
	
	nivh = -24	(mV)
	nik = 20.4
	
	ek
}

ASSIGNED {
	ik
	ninf
	ntau		(ms)
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
	ik = gkbar * n^4 * (v - ek)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n) / ntau
}

PROCEDURE rates(Vm (mV)) {
	LOCAL v
	v = Vm + 11	: Account for Junction Potential
	ninf = 1/(1+exp(-(v-nivh)/nik))
	ntau = 1000 * ntau_func(v)
}

FUNCTION ntau_func(v (mV)) {
	if (v < -20) {
		ntau_func = .000688 + 1/(exp((v+64.2)/6.5)+exp((v-141.5)/-34.8))
	} else {
		ntau_func = .00016 + .0008*exp(-.0267 * v)
	}
}
