TITLE leak current for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
ENDCOMMENT

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
}

NEURON {
 SUFFIX leak
 NONSPECIFIC_CURRENT i
 RANGE gmax, ileak
}

PARAMETER {
 v (mV)
 dt (ms)
 gmax  = 0.001 (mho/cm2)
 i  = 0.0 (mA/cm2)
 e (mV)
}

ASSIGNED { 
 ileak (mA/cm2)
}

BREAKPOINT {
 i  = gmax*(v-e)
 ileak = i
}
