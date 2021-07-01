COMMENT
Calcium accumulation into a volume of area*depth next to the
membrane with a decay (time constant tau) to resting level
given by the global calcium variable cai0_ca_ion

NOTES
-----

Depth parameter was calculated by equating the 'B' parameter for
calcium buffering in Gunay (2008) model to the multiplicative factor for ica
in equations below, specifically:

1/depth = 5.2e-12/6.283e-8 * 1e-6 * (2.0*Faraday/1e7)

The factor 1e-6 is obtained by adjusting the B factor in
http://genesis-sim.org/GENESIS/Hyperdoc/Manual-26.html#ss26.1
for NEURON units (concentration rate and current are in SI units in GENESIS)


CREDITS
-------

This is the example mechanism included with NEURON under 
/share/examples/nrniv/nmodl/examples/cacum.mod
ENDCOMMENT

NEURON {
    SUFFIX Calcium
    USEION ca READ ica WRITE cai
    RANGE depth, tau, cai0
}

UNITS {
    (mM) = (milli/liter)
    (mA) = (milliamp)
    F = (faraday) (coulombs)
}

PARAMETER {
    depth = 1.6e12 (nm)  : assume volume = area*depth
    tau = 1 (ms)
    cai0 = 50e-6 (mM)   : Requires explicit use in INITIAL
            : block for it to take precedence over cai0_ca_ion
            : Do not forget to initialize in hoc if different
            : from this default.
}

ASSIGNED {
    ica (mA/cm2)
}

STATE {
    cai (mM)
}

INITIAL {
    cai = cai0
}

BREAKPOINT {
    SOLVE integrate METHOD derivimplicit
}

DERIVATIVE integrate {
    cai' = -ica/depth/F/2 * (1e7) + (cai0 - cai)/tau
}
