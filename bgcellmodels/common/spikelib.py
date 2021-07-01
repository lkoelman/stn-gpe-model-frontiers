"""
Analysis and generation of spike trains.

"""

import numpy as np

def make_oscillatory_bursts(
    T_burst, dur_burst, f_intra, f_inter, 
    max_dur=10e3, rng=None):
    """
    Make oscillatory bursting spike times with given burst period, 
    burst duration, and intra- and inter-burst firing rates.

    The inter-burst interval (IBI) and burst duration are fixed so bursts
    will occur regularly with fixed intervals.

    Parameters
    ----------

    @param      T_burst : float
                Burst period (ms)

    @param      dur_burst : float
                Burst duration (ms)

    @param      f_intra : float
                Intra-burst firing rate (Hz)

    @param      f_inter : float
                Inter-burst firing rate (Hz)

    @param      rng : numpy.Random
                Random object (optional)

    Returns
    -------

    @return     Generator object that generates spike times


    Examples
    --------

    >>> import numpy
    >>> burst_gen = make_oscillatory_bursts(3500.0, 545.0, 50.0, 5.0)
    >>> bursty_spikes = numpy.fromiter(burst_gen, float)


    EXPERIMENTAL DATA
    -----------------

    Li 2012 - "Therapeutic DBS...":

        - T_burst = 3500 ms (f_burst = 0.286 Hz)
        - dur_burst = 545 ms
        - f_intra = 50 Hz
        - f_inter = 5 hz
    """
    if rng is None:
        rng = np.random
    # taken from example specific_network.py
    # Blend ISIs from two negexp distributions centered at intra- and
    # inter-burst firing rate, respectively
    T_intra, T_inter = 1e3/f_intra, 1e3/f_inter
    intra_ISIs = rng.exponential(T_intra, size=max(10, int(5*max_dur/T_intra)))
    inter_ISIs = rng.exponential(T_inter, size=max(10, int(5*max_dur/T_inter)))
    t = 0.0
    i_inter, i_intra = 0, 0
    while t < max_dur:
        j = t // T_burst # we are in the j-th cycle
        t0_cycle = j * T_burst
        in_burst = (t0_cycle <= t < (t0_cycle + dur_burst))
        if in_burst:
            t += intra_ISIs[i_intra]
            i_intra += 1
            yield t
        else:
            # Make sure we don't overshoot next burst period
            end_cycle = t0_cycle + T_burst 
            if t + inter_ISIs[i_inter] < (end_cycle + 0.5*dur_burst):
                t += inter_ISIs[i_inter]
                i_inter += 1
            else:
                # if ISI overshot burst period: pick time in first half of next burst
                t = end_cycle + 0.5*dur_burst*rng.uniform(0,1)
            yield t


def make_variable_bursts(
    T_burst, dur_burst, f_intra, f_inter, 
    max_dur=10e3, rng=None):
    """
    Make bursty spike train with variable inter-burst interval
    and burst duration. 

    The inter-burst interval (IBI) and burs duration are drawn
    from a negative exponential distribution so the burst rate is 
    Poisson-distributed.

    Arguments
    ---------

    @param      T_burst : float
                Burst period (ms), mean of negexp distribution

    @param      dur_burst : float
                Burst duration (ms), mean of negexp distribution

    @param      f_intra : float
                Intra-burst firing rate (Hz)

    @param      f_inter : float
                Inter-burst firing rate (Hz)

    @param      rng : numpy.Random
                Random object (optional)
    
    Returns
    -------
    
    @return     Generator object that generates spike times


    EXAMPLES
    --------

    >>> import numpy
    >>> burst_gen = make_variable_bursts(3500.0, 545.0, 50.0, 5.0)
    >>> bursty_spikes = numpy.fromiter(burst_gen, float)

    """
    if rng is None:
        rng = np.random
    # Blend ISIs from two negexp distributions centered at intra- and
    # inter-burst firing rate, respectively
    T_intra, T_inter = 1e3/f_intra, 1e3/f_inter
    
    # Pre-sample with safety margin of 5 (times faster than mean)
    max_num_burst = int(10*max_dur/T_burst)
    burst_IBIs = rng.exponential(T_burst, size=max(10, max_num_burst))
    # burst_durs = rng.exponential(dur_burst, size=max(10, max_num_burst))
    burst_durs = rng.uniform(max(5.0, dur_burst-10), dur_burst+10, size=max(10, max_num_burst))
    max_num_intra = int(10 * np.sum(burst_durs) / T_intra)
    max_num_inter = int(10 * np.sum(burst_IBIs) / T_inter)
    intra_ISIs = rng.exponential(T_intra, size=max(10, max_num_intra))
    inter_ISIs = rng.exponential(T_inter, size=max(10, max_num_inter))
    
    t_start_burst = burst_IBIs[0]
    t_end_burst = t_start_burst + burst_durs[0]

    # Indices into samples of interval distributions
    i_IBI, i_dur = 1, 1
    i_inter, i_intra = 0, 0

    in_burst = False
    t = 0.0    
    while t < max_dur:
        if t < t_start_burst:
            if in_burst and t < t_end_burst:
                # In burst: fire at f_intra
                t += intra_ISIs[i_intra]
                i_intra += 1
                yield t
            elif in_burst and t >= t_end_burst:
                # Reached end of burst: turn flag off
                in_burst = False
                i_dur += 1
                t_end_burst = t_start_burst + burst_durs[i_IBI]
                continue
            else:
                # Not in burst: fire at f_inter
                t += inter_ISIs[i_inter]
                i_inter += 1
                yield t
        elif t >= t_start_burst:
            # Transition to inside burst: turn flag on
            in_burst = True
            i_IBI += 1
            t_start_burst += burst_IBIs[i_IBI]
            continue


def generate_bursts_during(intervals, T_burst, dur_burst, f_intra, f_inter, 
                           f_background, duration, max_overshoot=0.25, rng=None):
    """
    Make spiketrain that bursts during specific intervals.
    """
    T_intra, T_inter, T_background = 1e3/f_intra, 1e3/f_inter, 1e3/f_background
    # TODO: check multiplicative properties of exponential distribution
    # normalized_ISI = rng.exponential(1.0, 
    #     size=max(10, int(2*duration/min(T_intra, T_inter, T_background))))
    i_isi = 0

    t = 0.0
    intervals_fifo = list(sorted(intervals)) + [None]
    ival = intervals_fifo.pop(0) # next or current interval
    while t < duration:
        while (ival is not None) and (t > ival[1]):
            ival = intervals_fifo.pop(0)
        if (ival is not None) and (ival[0] <= t < ival[1]):
            # we are inside bursty interval
            j = (t-ival[0]) // T_burst # we are in the j-th cycle
            t0_cycle = ival[0] + j * T_burst
            in_burst = (t0_cycle <= t < (t0_cycle + dur_burst))
            if in_burst:
                # We are intra-burst, during a bursty period
                t += rng.normal(T_intra, 0.1*T_intra)
                i_isi += 1
                yield t
            else:
                # We are inter-burst, during a bursty period
                end_cycle = t0_cycle + T_burst
                dt = T_inter * rng.exponential(T_inter)
                if t + dt < (end_cycle + max_overshoot*dur_burst):
                    t += dt
                    i_isi += 1
                else:
                    # if overshoot: randomize time inside burst period
                    t = end_cycle + max_overshoot*dur_burst*rng.uniform(0,1)
                yield t
        else:
            # We are outside bursty interval
            dt = rng.exponential(T_background)
            if (ival is not None) and (ival[0] <= t+dt < ival[1]):
                t = ival[0] + max_overshoot*dur_burst*rng.uniform(0,1)
            else:
                t += dt
                i_isi += 1
            yield t


def generate_modulated_spiketimes(
        sig_mod, Ts_mod, rate_max, rate_min, 
        dist='normal', rng=None, **kwargs):
    """
    Make modulated poisson spike train.

    The ISIs are drawn from a negative exponential distribution with mean
    equal to the modulator signal amplitude.

    @param  dist : str
            Distribution name, must be attribute of numpy.random

    @param  **kwargs
            Additional keyword arguments for numpy.random distribution function
            besides the 'loc' argument, e.g. 'scale' for normal distribution.
            WARNING: these must be values for loc=1.0 and will be scaled by
            the relative firing rate. Check scaling properties of distribution.

    Algorithm
    ---------

    Interprets the modulator signal as an amplitude-modulated sinusoid.
    The instantaneous firing rate is the slow rate + the fraction of the
    maximum amplitude reached times the difference of slow and fast rate during
    the peaks of the sinusoid. During the throughs the instantaneous firing
    rate is equal to the slow rate. I.e. the instantaneous amplitude during the
    negative peaks is discarded.

    Example
    -------


    """
    if rng is None:
        rng = np.random

    # Get distribution to sample from
    dist_func = getattr(rng, dist)

    mod_max = max(sig_mod)
    mod_min = min(sig_mod)
    mod_amp = (mod_max - mod_min) / 2.0
    mod_mid = (mod_min + mod_max) / 2.0

    Tfast = 1./rate_max * 1e3 # fast period (ms)
    Tslow = 1./rate_min * 1e3 # slow period (ms)
    Tdelta = Tfast - Tslow
    assert Tfast < Tslow

    # Pre-sample ISIs to guarantee statistical independence
    dur = len(sig_mod) * Ts_mod
    max_spikes = int(5*dur/Tfast)
    if dist == 'normal':
        kwargs['loc'] = 1.0
    if dist == 'exponential':
        kwargs['scale'] = 1.0
    ISIs_unscaled = dist_func(size=max(10, max_spikes), **kwargs)

    t = 0.0
    i = 0
    while t <= dur:
        yield t
        # fire at rate t_slow + frac_amp_pos * t_fast
        # - i.e. fire at Tslow when amp < 0
        factor = max(0, (sig_mod[int(t//Ts_mod)] - mod_mid) / mod_amp)
        Tinst = Tslow + factor*Tdelta
        t += ISIs_unscaled[i] * Tinst
        i += 1


def test_spiketime_generator(gen_func, num_spiketrains, *args, **kwargs):
    """
    Test case for generate_modulated_spiketimes()
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for i in range(num_spiketrains):

        generator = gen_func(*args, **kwargs)
        spiketimes = np.fromiter(generator, float)

        y_vec = np.ones_like(spiketimes) * i
        ax.plot(spiketimes, y_vec, marker='|', linestyle='', color='red', snap=True)

    ax.set_title("Modulated spiketimes")
    ax.grid(True)
    plt.show(block=False)


if __name__ == '__main__':
    T_burst, dur_burst, f_intra, f_inter = 50.0, 10.0, 150.0, 5.0
    num_spiketrains = 10
    test_spiketime_generator(make_oscillatory_bursts, num_spiketrains, 
                             T_burst, dur_burst, f_intra, f_inter, max_dur=5e3)