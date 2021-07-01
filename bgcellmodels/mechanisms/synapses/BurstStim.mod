: $Id: pregen.mod 3 2003-02-11 18:36:05Z hines $

COMMENT
Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a fast spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed), and 
either tonic or bursting.

Parameters;
   noise:   between 0 (no noise-periodic) and 1 (fully noisy)
   fast_invl:   fast interval, mean time between spikes (ms)
   slow_invl:   slow interval, mean burst silent period (ms), 0=tonic train
   burst_len:   mean burst length (nb. spikes)

Usage:
   see examples in nrn/examples/nrniv/netcon/ -> search for 'SpikeGenerator'

Written by Z. Mainen, modified by A. Destexhe, The Salk Institute

Modified by Michael Hines for use with CVode

Modified by Michael Hines to use logical event style with NET_RECEIVE

Modified by Lucas Koelman: 
    - renamed SpikeGenerator to BurstStim, restructure comments
    - changed mechanism type from POINT_PROCESS to ARTIFICIAL_CELL
    - changed exprand(1) call to custom function erand() copied from netstim.mod
ENDCOMMENT


NEURON  { 
  ARTIFICIAL_CELL BurstStim
  RANGE y
  RANGE fast_invl, slow_invl, burst_len, start, end
  RANGE noise
  POINTER donotuse
}


PARAMETER {
    fast_invl   = 1 (ms)    : time between spikes in a burst (msec)
    slow_invl   = 50 (ms)   : burst period (msec)
: actually, above is interburst period in conformity with original version
: see
    burst_len   = 10        : burst length (# spikes)
    start       = 50 (ms)   : start of first interburst interval
    end     = 1e10 (ms) : time to stop bursting
    noise       = 0     : amount of randomeaness (0.0 - 1.0)
}


ASSIGNED {
    y
    burst
    event (ms)
    burst_off (ms)
    burst_on (ms)
    toff (ms)
    on
    donotuse
}


:backward compatibility
PROCEDURE seed(x) {
    set_seed(x)
}

INITIAL {

    if (noise < 0) {
        noise = 0
    }
    if (noise > 1) {
        noise = 1
    }

    on = 1
    toff = 1e9
    y = -90
    burst = 0
    event = start - slow_invl
    event_time()
    while (on == 1 && event < 0) {
        event_time()
    }
    if (on == 1) {
        net_send(event, 1)
    }
}   

FUNCTION interval(mean (ms)) (ms) {
    : sample an ISI from distribution with given mean
    :
    : @return   float : sample from distribution with given mean
    
    if (mean <= 0.) {
        mean = .01 (ms) : I would worry if it were 0.
        : since mean is a local variable, if the number it is set
        : to is dimensionless, mean will be dimensionless.
    }
    if (noise == 0) {
        interval = mean
    }else{
        interval = (1. - noise)*mean + noise*mean*erand()
    }
}

VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM


: sample negexp distribution from specified random stream

FUNCTION erand() {
VERBATIM
    if (_p_donotuse) {
        /*
        :Supports separate independent but reproducible streams for
        : each instance. However, the corresponding hoc Random
        : distribution MUST be set to Random.negexp(1)
        */
        _lerand = nrn_random_pick(_p_donotuse);
    }else{
        /* only can be used in main thread */
ENDVERBATIM
        : the old standby. Cannot use if reproducible parallel sim
        : independent of nhost or which host this instance is on
        : is desired, since each instance on this cpu draws from
        : the same stream
        erand = exprand(1)
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


PROCEDURE event_time() {
    : Update state variables
    :
    : @post     'burst' state var is updated
    :
    : @post     'event' state var is updated
    :
    : @post     'burst_on' state var is updated
    :
    : @post     'burst_off' state var is updated

    if (slow_invl == 0 || (burst != 0. && burst_len > 1)) {
        : We are in a burst ('spiking mode')
        event = event + interval(fast_invl)
        if (event > burst_on + burst_off) {
            burst = 0.
        }
    }else{
        : We are not in burst ('silent mode')
        burst = 1.
        
        : Use following code if slow_invl represents beginning of burst to beginning of burst
        : > event = event + interval(slow_invl - (burst_len-1)*fast_invl)
        : Use following if slow_invl represents interburst interval (end to beginning)
        event = event + interval(slow_invl)
        
        burst_on = event
        burst_off = interval((burst_len - 1)*fast_invl)-1e-6
    }
    if (event > end) {
        on = 0
    }
}

NET_RECEIVE (w) {
:printf("Pregen receive t=%g flag=%g\n", t, flag) 
    if (flag == 1 && on == 1) {
        y = 20
        net_event(t)
        event_time()
        net_send(event - t, 1)
        net_send(.1, 2)
    }
    if (flag == 2) {
        y = -90
    }
}

