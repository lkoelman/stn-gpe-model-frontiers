TITLE calcium dynamics for GPe neuron

COMMENT

ENDCOMMENT

NEURON {
 SUFFIX Calcium
 USEION ca READ ica WRITE cai
 RANGE a2v
}

UNITS {
 (mM) = (milli/liter)
 (mA) = (milliamp)
 (um) = (microm)
 F = (faraday)(coulomb)
}

PARAMETER {
 dt (ms)
 z = 2.0
 kca = 0.4 (1/ms)
 cai0 = 0.00001 (mM)
 a2v = 3000.0 (1/cm)
}

ASSIGNED {
 ica (mA/cm2)
 alpha
}

STATE {
 cai (mM)
}

BREAKPOINT {
 SOLVE integrate METHOD cnexp
}

UNITSOFF

INITIAL{
 alpha = a2v/(z*F)
 cai = cai0
}

DERIVATIVE integrate {
 cai' = -alpha*ica - kca*(cai - cai0)
}

UNITSON
