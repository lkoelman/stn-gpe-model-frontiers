TITLE passive membrane properties for STh

COMMENT
Temperature is going to effect the passive membrane properties.  How to
model this is difficult.  Kita do not recorded their temperature at
which they measure the membrane time constant, and Nakanishi who
measure input resistance perform their experiments at 36degC.
At present I will not implement a temperature dependence of this
component. 

ENDCOMMENT

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX STh
	RANGE gpas, epas
	NONSPECIFIC_CURRENT ipas
}

PARAMETER {
        v (mV)
	dt (ms)
        gpas  = 7.84112e-05 (mho/cm2) <0,1e9>
	epas  = -58.4477 (mV)
	celsius
}

ASSIGNED { 
	ipas (mA/cm2)
}

BREAKPOINT {
	ipas = gpas*(v - epas)
}
