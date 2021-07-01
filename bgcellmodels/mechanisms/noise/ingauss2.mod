COMMENT
Noise current characterized by gaussian distribution with mean 'mean' and 
standard deviation 'stdev'.

Borrows from NetStim's code so it can be linked with an external instance 
of the Random class in order to generate output that is independent of 
other instances of InGau.

This mechanism is intended for use with fixed time steps only, since it generates a new current value at each advance of the simulation.

If you need to inject noise current into more than one segment of a model, attach a separate instance of InGauss to each segment and scale its standard deviation according to the compartment's surface area (assuming you want uniform noise current per unit of surface area).

To ensure that each noise source is independent of every other noise source, create a separate Hoc Random instance for each instance of InGauss.


EXAMPLE
-------

Python example:

>>> sec = h.Section()

>>> rng = h.Random()
>>> rng.MCellRan4(high_index, low_index) # high_index can also be set using rng.seq()
>>> rng.normal(0, 1)

>>> stim = h.ingauss2(sec(0.5))
>>> noise_amp = 0.01 # [mA/cm2], values of 0.005 - 0.01 seem reasonable
>>> stim.mean = 0.0
>>> stim.stdev = noise_amp * 1e-2 * sum((seg.area() for seg in sec))
>>> stim.noiseFromRandom(rng)


CREDITS
-------

Adapted by Lucas Koelman from the mechanism created by Ted Carnevale,
posted at https://www.neuron.yale.edu/phpBB/viewtopic.php?t=2986

ENDCOMMENT

NEURON {
    POINT_PROCESS ingauss2
    NONSPECIFIC_CURRENT i
    RANGE mean, stdev
    THREADSAFE : true only if every instance has its own distinct Random
    POINTER donotuse
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    mean = 0 (nA)
    stdev = 1 (nA)
}

ASSIGNED {
    i (nA)
    donotuse
}

INITIAL {
    i = 0
}

PROCEDURE seed(x) {
    : for backward compatibility with old method normrand()
    set_seed(x)
}

: WARNING: This method does not seem to work in NEURON 7.5
: BEFORE BREAKPOINT {
:     : dummy BREAKPOINT block must be present if this block is used
:     i = stdev*grand() + mean
: }

BREAKPOINT {
    SOLVE dummy   :has to be in a proc otherwise there is an error
}

PROCEDURE dummy() {
    i = stdev*grand() + mean 
}


VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM


FUNCTION grand() {
VERBATIM
    if (_p_donotuse) {
        /*
         : Supports separate independent but reproducible streams for
         : each instance. However, the corresponding hoc Random
         : distribution MUST be set to Random.normal(0,1)
         */
            _lgrand = nrn_random_pick(_p_donotuse);
    }else{
        /* only can be used in main thread */
        if (_nt != nrn_threads) {
hoc_execerror("multithread random in InUnif"," only via hoc Random");
        }
ENDVERBATIM
        : the old standby. Cannot use if reproducible parallel sim
        : independent of nhost or which host this instance is on
        : is desired, since each instance on this cpu draws from
        : the same stream
:        erand = exprand(1)
:        urand = scop_random()
        grand = normrand(0,1)
: printf("%f\n", grand)
VERBATIM
    }
ENDVERBATIM
}


PROCEDURE noiseFromRandom() {
VERBATIM
 {
    void** pv = (void**)(&_p_donotuse);
    if (ifarg(1)) {
        *pv = nrn_random_arg(1);
    }else{
        *pv = (void*)0;
    }
 }
ENDVERBATIM
}
