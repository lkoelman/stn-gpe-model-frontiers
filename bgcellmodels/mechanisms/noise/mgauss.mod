COMMENT
Specific membrane noise current characterized by gaussian distribution 
with mean 'mean' and standard deviation 'stdev'.

It is sufficient to use a distinct Random instance per cell (different cells 
can be in different threads but all the compartments of any given cell are 
in the same thread, assuming no cell splitting).


EXAMPLE
-------

Python example:

>>> sec = h.Section()
>>> sec.nseg = 5

>>> rng = h.Random()    # one instance per cell is sufficient
>>> rng.MCellRan4(high_index, low_index)
>>> rng.uniform(0, 1)   # distribution MUST be set to uniform(0,1)

>>> sec.insert('mgauss')
>>> for seg in sec:
>>>     seg.noiseFromRandom_noise(rng)
>>>     seg.mean_noise = 0
>>>     seg.stdev_noise = 1

CREDITS
-------

Implemented by Lucas Koelman

Adapted from ingauss.mod posted by Ted Carnevale at
Posted by Ted Carnevale at https://www.neuron.yale.edu/phpBB/viewtopic.php?t=2986

ENDCOMMENT

TITLE Membrane noise

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

UNITS {
    (mA) = (milliamp)
}

NEURON {
    SUFFIX mgauss
    NONSPECIFIC_CURRENT i
    RANGE mean, stdev
    THREADSAFE          : true only if every instance has its own distinct Random
    POINTER donotuse    : pointer to random number generator
}

PARAMETER {
    mean    = 0 (mA/cm2)
    stdev   = 1 (mA/cm2)
    : area         (um2) : (not assigned) segment area in micron, used in POINT_PROCESS current
}

ASSIGNED {
    i           (mA/cm2)
    donotuse
}


:backward compatibility
PROCEDURE seed(x) {
    set_seed(x)
}

INITIAL {
    i = 0
}


BREAKPOINT {
    : NOTE: ingauss.mod used BEFORE BREAKPOINT block
    SOLVE dum   :has to be in a proc otherwise there is an error
}

PROCEDURE dum() {
    i = stdev*grand() + mean 
}


VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

: sample Gaussian standard normal distribution
FUNCTION grand() {
VERBATIM
    if (_p_donotuse) {
        /*
         : Supports separate independent but reproducible streams for
         : each instance. However, the corresponding hoc Random
         : distribution MUST be set to Random.uniform(0,1)
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