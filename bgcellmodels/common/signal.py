"""
Signal analysis tools for electrophysiology data.

@author     Lucas Koelman

Notes
-----

See features in PyElectro: https://pyelectro.readthedocs.io/en/latest/pyelectro.html

See features in eFEL: https://efel.readthedocs.io/en/latest/eFeatures.html
"""

from __future__ import division # always float division, never trunc to int

import neuron
h = neuron.h
import numpy as np
import scipy.signal
import scipy.stats


def spike_indices(v_rec, v_th, loc='onset'):
    """
    Get spike indices in voltage trace.

    @param  v_rec : np.array
            membrane voltage

    @param  v_th : float
            voltage threshold for spike

    @param  loc : str
            'onset' : time of spike onset, i.e. where threshold is first reached
            'offset' : time of spike offset, i.e. where v dips below threshold
    """
    thresholded = np.array((v_rec - v_th) >= 0, dtype=float)
    i_onset = np.where(np.diff(thresholded) == 1)[0] + 1
    i_offset = np.where(np.diff(thresholded) == -1)[0] + 1
    if loc == 'onset':
        return i_onset
    elif loc == 'offset':
        return i_offset
    else:
        raise ValueError("Unrecognized value for argument 'loc' : {}".format(loc))


def spike_times(v_rec, t_rec, v_th, loc='onset'):
    """
    Get spike times in voltage trace.

    @see    spike_indices(...)

    @param  t_rec : np.array
            Sample times of v_rec.
    """
    spike_idx = spike_indices(v_rec, v_th, loc=loc)
    return t_rec[spike_idx]


def coefficient_of_variation(v_rec, t_rec, v_th):
    """
    Calculate coefficient of variation of spikes.
    This is the ratio of the standard deviation and the mean of the ISIs.
    """
    t_spikes = spike_times(v_rec, t_rec, v_th, loc='onset')
    isi_vals = np.diff(t_spikes)
    return np.std(isi_vals) / np.mean(isi_vals)


def burst_metrics(
        v_rec, t_rec, threshold=-20.0,
        onset_isi=20.0, offset_isi=20.0, min_spk=4):
    """
    Identify bursts using simple algorithm where a burst consists of a minimum
    of <min_spk> spikes where the first ISI is smaller than <onset_isi> and
    susbequent ISIs are smaller than <offset_isi>.

    After burst identification, calculate burst rate, inter-burst intervals,
    intra-burst rates, number of spikes per burst.

    @return     dict[str, <float or list[float]>]

                Dictionary containing burst metrics:

                'spikes_per_burst':         number of spikes in each burst,
                'intra_burst_rates':        firing rate inside each burst,
                'inter_burst_intervals':    time intervals between bursts,
                'burst_rate':               number of bursts per second,
                'ISI_CV':                   coefficient of variation,

    """
    # Using eFEL: use peak indices
    # FIXME: indices are half the values of our own function ???
    # features = ['ISI_values', 'peak_indices']
    # feat_vals = get_efel_features(v_rec, t_rec, interval, features, threshold=threshold)
    # isi_vals = feat_vals['ISI_values']
    # peak_idx = feat_vals['peak_indices']

    # Using threshold crossing indices
    peak_idx = spike_indices(v_rec, v_th=threshold, loc='onset')
    peak_times = t_rec[peak_idx]
    isi_vals = np.diff(peak_times)

    burst_lengths = []
    burst_intra_rates = []
    onset_times = []
    offset_times = []
    inter_burst_intervals = []

    # Find bursts according to absolute criteria
    num_isi = len(isi_vals)
    isi_burst = np.zeros_like(isi_vals, dtype=bool) # flag which ISIs are in a burst
    isi_below_offset = np.array(isi_vals <= offset_isi, dtype=int)
    candidate_offsets = np.diff(isi_below_offset) == -1
    i = 0
    while True:
        if i > num_isi-1:
            break
        isi = isi_vals[i]
        if isi <= onset_isi:
            # i_offset == 0 if i itself is the last ISI of a burst
            offset_dists, = np.where(candidate_offsets[i:])
            if len(offset_dists) == 0:
                offset_dist = num_isi - i - 1
            else:
                offset_dist = offset_dists[0]

            if offset_dist+2 < min_spk:
                i += 1
                continue

            # Now still have to check if all ISI until candidate offset are < offset_ISI
            num_to_offset = np.sum(isi_below_offset[i:i+offset_dist+1])
            if (num_to_offset+1 < min_spk) or (num_to_offset != offset_dist+1):
                i += 1
                continue
            num_follow = num_to_offset

            # Burst properties
            burst_lengths.append(num_follow+1)
            f_intra = 1e3/np.mean(isi_vals[i:i+num_follow])
            if np.isnan(f_intra):
                raise Exception('NaN!')
            burst_intra_rates.append(f_intra)
            t_onset = t_rec[peak_idx[i]]
            t_offset = t_rec[peak_idx[min(i+num_follow, num_isi-1)]]
            # print("Burst at t={} with ISIs: {}".format(t_onset, isi_vals[i:i+num_follow]))
            if len(offset_times) > 0:
                inter_burst_intervals.append(t_onset - offset_times[-1])
            onset_times.append(t_onset)
            offset_times.append(t_offset)
            # Skip to end of burst
            isi_burst[i:i+num_follow+1] = True
            i += num_follow + 1
        else:
            i += 1

    burst_rate = len(burst_lengths) * 1e3 / (t_rec[-1] - t_rec[0])
    metrics = {
        'spikes_per_burst': burst_lengths,
        'intra_burst_rates': burst_intra_rates,
        'inter_burst_intervals': inter_burst_intervals,
        'burst_rate': burst_rate,
        'ISI_CV': np.std(isi_vals) / np.mean(isi_vals),
    }
    return metrics


def burst_metrics_surprise(
        ISI,
        sampling_rate=1000.0,
        min_rate_factor=2.0,
        grow_rate_factor=0.5,
        min_burst_length=3,
        local_length=0,
        surprise_cutoff=3.0):
    """
    Identify bursts using Poisson 'surprise' method (Legendy and Salcman 1985)
    and calculate various burst metrics.

    CONTRACT
    --------

    @pre        Code relies on float division, i.e. in Python2 this function
                definition must be preceded by 'from __future__ import division'

    @param      sampling_rate : float
                Unit conversion factor for ISI times : <time units>/seconds.
                E.g. 1000.0 if ISI values are in units of ms.

    @param      min_rate_factor : float
                As the initial criterion for burst detection, spikes have to occur at a
                frequency which is higher than the baseline frequency by a factor
                'min_rate_factor'.

    @param      min_burst_length : int
                Minimum number of spikes that may be considered as a burst

    @param      surprise_cutoff : float
                Minimum value of the Poisson surprise index S = - log(P)
                for a burst.

    @param      local_length : float
                parameter for calculating the baseline firing rate for comparison
                with burst firing rate. A value of 0 indicates that the entire spike train should be used to yield the firing rate (can be used in case of
                homogeneous firing of units). A value > 0 idicates the length of
                the time interval prior to the burst.

    @return     dict[str, float / list[float]]
                Dictionary containing burst metrics:

        Per-burst metrics:

        'begin'         -> list[int]            : start index for each burst
        'num_spikes'    -> list[int]            : number of spikes in each burst
        'surprise'      -> list[float]          : surprise value for each burst
        'rate'          -> list[float]          : mean firing rate in each burst
        'max_rate'      -> list[float]          : max firing rate in each burst
        'baseline_rate' -> list[float]          : baseline rate used to threshold
                                                  burst firing rate
        'inter_burst_intervals' -> list[float]  : times between bursts

        Summary metrics:

        'num_bursts' -> int                     : ~
        'mean_spikes_per_burst' -> float        : ~
        'median_spikes_per_burst' -> float      : ~
        'total_spikes_in_bursts' -> int         : ~
        'mean_intra_burst_frequency' -> float   : ~
        'median_intra_burst_frequency' -> float : ~
        'proportion_time_in_bursts' -> float    : ~
        'mean_intra_burst_frequency' -> float   : ~
        'proportion_spikes_in_bursts' -> float  : ~
        'burst_rate' -> float                   : bursts per second (Hz)


    ALGORITHM
    ---------

    (Description from Legend & Salcman 1985)

    The measure used here is an evaluation of how improbable it is that the burst is
    a chance occurrence and is computed, for any given burst that contains n spikes
    in a time interval T, as S = -log P where P is the probability that, in a random
    (Poisson) spike train having the same average spike rate r as the spike train
    studied, a given time interval of length T contains n or more spikes. P is given
    by Poisson's formula, as P = exp(-rT) * sum_{i=n}^{i=inf}((rT)^i / i!). S will
    be reffered to as the Poisson surprise of the burst.

    Burst analysis of a spike train began with a quick pass through the [spike
    train], to evaluate the buffer-wide average spike rate [...].  Then, in a second
    pass, the program advanced down the spike train, spike by spike, until it found
    several closely spaced spikes; when it did, it called them a burst. (The
    required number of such spikes had been chosen to be 3, and their largest
    allowed average spacing was chosen to be half the buffer- wide average.) Next,
    the program tested whether including an additional spike at the end of the burst
    increased the Poisson surprise of the burst, and if so, included the additional
    spike. It then repeated the process until either a relatively long interval was
    encountered or inclusion of one or more (up to 10) additional spikes failed to
    increase the Poisson surprise. (The length of such a relatively long interval
    had been chosen to be twice the buffer-wide average.) Next, the program tested
    whether the Poisson surprise could be further increased by removing one or more
    spikes from the beginning of the burst; if this was possible, the program
    removed whatever number of spikes would maximize the Poisson surprise. This con-
    cluded the burst-defining process after

    REFERENCES
    ----------

    - Legendy & Salcman 1985
        - min length of burst = 3
        - surprise value threshold = 10
        - min firing rate for burst detection = 2 x baseline rate
        - min firing rate for burst growing = 0.5 x baselinerate

    Example usage of method, with different algorithm parameters:

    - Wichmann and Soares 2006
        - surprise value threshold = 3
        - min firing rate for burst growing = 2 x baseline rate ? (in MATLAB code)

    - Hahn et al. (2008) - Experimental neurology, 211(1), 243-251.
        - surprise value threshold = 3
        - min firing rate for burst detectio = < 1 x baseline rate
        - burst growing in both directions, instead of shrinking from front

    - Sanders et al. 2013, Sharott et al. 2014
        - surprise value threshold = 3


    CREDITS
    -------

    Algorithm from MATLAB version implemented by Thomas Wichmann found at
    https://github.com/brian-lau/matutils/blob/master/+spk/burst.m

    """
    burst_metrics = {
        'begin': [],
        'num_spikes': [],
        'surprise': [],
        'rate': [],
        'max_rate': [],
        'baseline_rate': [],
        'num_bursts': [],
        'mean_spikes_per_burst':[],
        'median_spikes_per_burst':[],
        'total_spikes_in_bursts':[],
        'mean_intra_burst_frequency':[],
        'median_intra_burst_frequency':[],
        'proportion_time_in_bursts':[],
        'proportion_spikes_in_bursts':[]
    }

    CA = np.cumsum(ISI)

    if local_length == 0:
        # entire data stream should be used to yield the firing rate for comparison
        mean_FR = len(ISI) / (np.sum(ISI) / sampling_rate)
        beg_idx = -1
    else:
        # finds last index within the 'local length' - incremented by one, this will result in the first index that can be evaluate.
        beg_idx = np.where(CA < local_length * sampling_rate)[0][-1]

    n = beg_idx

    while n < (len(ISI) - min_burst_length):

        # n is a running parameter that points to ISIs
        n += 1

        if local_length > 0:
            # find the ISI data segment I that is fully contained within the local_length
            start_idx = np.where(CA > (CA[n] - local_length * sampling_rate))[0][0]
            I = ISI[start_idx:n]

            # calculation of frequency threshold
            mean_FR = len(I) / (np.sum(I) / sampling_rate)

        # Firing rate thresholds in units of time (e.g. ms or s)
        fr_thr = sampling_rate / (min_rate_factor * mean_FR) # units of time (e.g. ms or s)
        grow_fr_thr = sampling_rate / (grow_rate_factor * mean_FR)

        # 1. DETECTION PHASE ===================================================
        # find runs of short ISIs
        if ISI[n] < fr_thr:
            # find areas of the spike train which fulfill the length_of_burst criterion

            q = 0 # running parameter that points to the number of spikes to be added

            while (n + q < len(ISI)) and (ISI[n+q] < fr_thr):
                q += 1

            if q > 0:
                q -= 1

            # at this point, the provisional burst starts at n and ends at n+q;
            # it has q + 1 spikes in it

            # 2. REFINEMENT PHASE ==============================================
            # adjust length of burst to maximize surprise value
            if q + 1 >= min_burst_length:
                m = min_burst_length

                while ((n + m < len(ISI)) and
                       (ISI[n + m] < fr_thr) and
                       (surprise(mean_FR, ISI[n:n+m+1], sampling_rate) >=
                        surprise(mean_FR, ISI[n:n+m], sampling_rate))):

                    m += 1

                if m > min_burst_length:
                    m -= 1

                # at this point, the beginning of the burst is provisionally settled to be n, the end n+m
                # the burst has m+1 spikes in it.

                # 3. BURST GROWING PHASE ==========================================
                # test whether adding up to 10 more ISIs will enhance surprise value
                if n + m + 10 < len(ISI):
                    # mmax is set to 10 unless one couldn't add 10 to n before reaching the end of FR
                    mmax = 10
                else:
                    mmax = len(ISI) - (n + m)

                # looking for 'slow spikes' within the next 10 spikes after the current burst end
                ind_long_ISI = np.where(ISI[n+m+1:n+m+mmax+1] > grow_fr_thr)[0]

                # pmax is set to be the index of the slow spike
                if len(ind_long_ISI) > 0:
                    pmax = ind_long_ISI[0] # max(0, ind_long_ISI[0] - 1)
                else:
                    pmax = mmax

                # formation of an array that will contain surprise values.
                # The first one is that of the burst defined by the ISIs between n and n+m.
                # Additional entries into S are surprise values that would result
                # from adding up to pmax additional spikes (S2 will be the surprise
                # values for ISI(n:n+m+1), S3 the one for ISI(n:n+m+2) etc.)
                S = np.zeros(pmax + 1)
                S[0] = surprise(mean_FR, ISI[n:n+m+1], sampling_rate)

                #  forms array of surprise values for this burst, starting from the end of the burst to pmax values later
                for p in range(1, pmax + 1):
                    S[p] = surprise(mean_FR, ISI[n:n+m+p+1], sampling_rate)

                ind_max_S = np.argmax(S)

                if n+m < len(ISI):
                    # this will set the m-value to the previous m, if the first
                    # entry into the S array is the highest (i.e., if adding additional
                    # ISIs didn't increase the surprise value), or, it will correct m
                    # to the highest value
                    m += ind_max_S
                else:
                    m = len(ISI) - n

                # at this point, the end of the index of the end of the burst
                # is settled to be n+m


                # 4. BURST SHRINKING PHASE =====================================
                # test whether adjusting the front end of the burst enhances the surprise value
                if n > 1:
                    o = 1

                    while ((m - o > min_burst_length) and
                           (surprise(mean_FR, ISI[n+o:n+m+1], sampling_rate) >=
                            surprise(mean_FR, ISI[n+o-1:n+m+1], sampling_rate))):
                        o += 1

                    if o > 1:
                        o -= 1 # reducing o by one to correct for the addition that resulted in the end of the while loop
                        n += o # adjust the beginning of the burst
                        m -= o # adjust the length of the burst

                # at this point, the beginning of the burst is settled to be n, and the length is m+1

                if ((m + 1 >= min_burst_length) and
                    (surprise(mean_FR, ISI[n:n+m+1], sampling_rate) > surprise_cutoff)):

                    burst_ISIs = ISI[n:n+m+1]

                    burst_metrics['begin'].append(n)
                    burst_metrics['num_spikes'].append(m+1)
                    burst_metrics['surprise'].append(
                        surprise(mean_FR, burst_ISIs, sampling_rate))
                    burst_metrics['rate'].append(
                        len(burst_ISIs) / (np.sum(burst_ISIs) / sampling_rate))
                    burst_metrics['max_rate'].append(sampling_rate / min(burst_ISIs))
                    burst_metrics['baseline_rate'].append(mean_FR)

                # adjust ISI pointer to the ISI following the burst
                n += (m + 1)

        # end if ISI[n] < fr_thr:
    # end while n < (len(ISI) - min_burst_length)

    # Calculate aggregate metrics for all bursts
    burst_metrics['num_bursts'] = len(burst_metrics['begin'])
    burst_metrics['inter_burst_intervals']= []

    if burst_metrics['num_bursts'] > 0:
        all_num_spikes = burst_metrics['num_spikes']
        burst_metrics['mean_spikes_per_burst'] = np.mean(all_num_spikes)
        burst_metrics['median_spikes_per_burst'] = np.median(all_num_spikes)
        burst_metrics['total_spikes_in_bursts'] = np.sum(all_num_spikes)

        all_rates = burst_metrics['rate']
        burst_metrics['mean_intra_burst_frequency'] = np.mean(all_rates)
        burst_metrics['median_intra_burst_frequency'] = np.median(all_rates)
        burst_metrics['proportion_time_in_bursts'] = burst_metrics['total_spikes_in_bursts'] / burst_metrics['mean_intra_burst_frequency'] / (np.sum(ISI[beg_idx+1:len(ISI)]) / sampling_rate)
        burst_metrics['proportion_spikes_in_bursts'] = burst_metrics['total_spikes_in_bursts'] / (len(ISI) - beg_idx)

        # burst rate is sampling_rate / mean(inter_burst_intervals)
        if burst_metrics['num_bursts'] >= 2:
            for i_burst, i_begin in enumerate(burst_metrics['begin']):
                if i_burst == burst_metrics['num_bursts']-1:
                    break
                len_burst = burst_metrics['num_spikes'][i_burst]
                i_next = burst_metrics['begin'][i_burst+1]
                sum_inter_ISIs = sum(ISI[i_burst+len_burst:i_next])
                burst_metrics['inter_burst_intervals'].append(sum_inter_ISIs)

            # burst rate is num_bursts / (interval[1] - interval[0])
            # burst_metrics['burst_rate'] = sampling_rate / np.mean(sum_inter_ISIs)
    else:
        burst_metrics['mean_spikes_per_burst'] = 0.
        burst_metrics['median_spikes_per_burst'] = 0.
        burst_metrics['total_spikes_in_bursts'] = 0.
        burst_metrics['mean_intra_burst_frequency'] = 0.
        burst_metrics['median_intra_burst_frequency'] = 0.
        burst_metrics['proportion_time_in_bursts'] = 0.
        burst_metrics['proportion_spikes_in_bursts'] = 0.

    return burst_metrics


def surprise(r, data, sampling_rate, with_deceleration=False):
    """
    Calculate surprise index.

    @param  r : float
            comparison firing rate (spikes per second)

    @param  data : np.array[float]
            ISI values to be included in the burst

    @param  sampling_rate : float
            Sampling rate of ISI measurements
    """

    T = np.sum(data) / sampling_rate
    num_spikes = len(data)

    p = scipy.stats.poisson.cdf(num_spikes, r*T)

    if p == 0:
        burst = 0
        deceleration = 100
    elif p == 1:
        burst = 100
        deceleration = 0
    else:
        burst = -np.log(1-p)
        deceleration = -np.log10(p)

    if with_deceleration:
        return burst, deceleration
    else:
        return burst



def numpy_sum_psth(spiketrains, tstart, tstop, binwidth=10.0, average=False):
    """
    Sum peri-stimulus spike histograms (PSTH) of spike trains

    @param spiketimes   list of Hoc.Vector() objects or other iterables,
                        each containing spike times of one cell.

    @param tstart       stimulus time/start time for histogram

    @param tstop        stop time for histogram

    @param binwidth     bin width (ms)

    @return             hoc.Vector() containing binned spikes
    """
    # Compute bin edges
    num_bins = int((tstop-tstart)/binwidth)
    edges = np.linspace(tstart, tstop, num_bins+1)

    # Concatenate spike trains and compute histogram
    sts_concat = np.concatenate(spiketrains)
    bin_vals, bin_edges = np.histogram(sts_concat, edges, density=False)

    if average:
        bin_vals /= len(spiketrains)

    return bin_vals


def numpy_psth_triggered(spiketrains, trigger_times, bins, interval=None):
    """
    Peri-stimulus time histogram with stimulus times

    @param  spike_trains : iterable[np.array[float]]
            List of spike trains

    @param  trigger_times : np.array[float]
            Stimulus times (= trigger times)

    @param  binds : np.array[float] or int
            Bin edges or number of bins (see numpy.histogram arg 'bins')
    """
    num_trigger = len(trigger_times)

    if interval is not None:
        mask = (trigger_times >= interval[0]) & (trigger_times <= interval[1])
        trigger_times = trigger_times[mask]

    individual_counts = [] # PSTH (counts) for each spike train

    for spike_times in spiketrains:
        mask = (spike_times > trigger_times[0])
        if interval is not None:
            mask = mask & (spike_times >= interval[0]) & (spike_times <= interval[1])
        times = np.array(spike_times[mask]) # copy for in-place modification

        # Normalize all spike times to preceding trigger
        for i, t in enumerate(trigger_times):
            mask = times > t
            if i < (num_trigger - 1):
                mask = mask & (times <= trigger_times[i+1])
            times[mask] = times[mask] - t

        counts, edges = np.histogram(times, bins=bins)
        individual_counts.append(counts)

    summed_counts = sum(individual_counts) # sum of numpy arrays

    return summed_counts, individual_counts, edges



def numpy_avg_rate_simple(spiketrains, tstart, tstop, binwidth):
    """
    Simple algorithm for calculating the running firing rate
    (average firing rate in each bin of the psth, averaged over spike trains)
    """
    # first compute the histogram
    avghist = numpy_sum_psth(spiketrains, tstart, tstop, binwidth)

    # divide by nb. of cells/trials and by binwidth in ms to get rate
    return avghist / (binwidth*1e-3*len(spiketrains))


def moving_average_rate(spiketrains, tstart, tstop, dt, binwidth):
    """
    Moving average firing rate with resolution dt, in interval
    [t-bin/w, t+bin/2]
    """
    all_spiketimes = np.concatenate(spiketrains)
    num_spiketrains = len(spiketrains)

    time = np.arange(tstart, tstop+dt, dt)
    moving_rate = np.zeros_like(time)
    for i, t in enumerate(time):
        t0 = t - binwidth / 2.0
        t1 = t + binwidth / 2.0
        if t0 < tstart:
            t0 = tstart
        if t1 > tstop:
            t1 = tstop
        true_binwidth = t1 - t0
        bin_num_spikes = np.sum((all_spiketimes > t0) & (all_spiketimes < t1))
        moving_rate[i] = float(bin_num_spikes) / true_binwidth / num_spiketrains

    return moving_rate


def nrn_sum_psth(spiketrains, tstart, tstop, binwidth=10.0, average=False):
    """
    Sum peri-stimulus spike histograms (PSTH) of spike trains

    @param spiketimes   list of Hoc.Vector() objects or other iterables,
                        each containing spike times of one cell.

    @param tstart       stimulus time/start time for histogram

    @param tstop        stop time for histogram

    @param binwidth     bin width (ms)

    @return             hoc.Vector() containing binned spikes
    """

    # Create histogram with empty bins
    # avghist = h.Vector(int((tstop-tstart)/binwidth) + 2, 0.0)
    st1 = h.Vector(spiketrains[0])
    avghist = st1.histogram(tstart, tstop, binwidth)

    # add histogram for each cell (add bins element-wise)
    for st in spiketrains[1:-1]:
        if not isinstance(st, neuron.hoc.HocObject):
            st = h.Vector(st)

        spkhist = st.histogram(tstart, tstop, binwidth)
        avghist.add(spkhist) # add element-wise

    # divide by nb. of cells/trials and by binwidth in ms to get rate
    if average:
        avghist.div(len(spiketrains))
    return avghist


def nrn_avg_rate_simple(spiketrains, tstart, tstop, binwidth):
    """
    Simple algorithm for calculating the running firing rate
    (average firing rate in each bin of the psth, averaged over spike trains)

    @param      binwidth: float
                Bin width in (ms)

    @return     h.Vector() containing, the firing rate in each bin.
                The firing rate is just the number of spikes in a bin
                divided by the bin width in seconds.
    """
    # first compute the histogram
    avghist = nrn_sum_psth(spiketrains, tstart, tstop, binwidth)

    # divide by nb. of cells/trials and by binwidth in ms to get rate
    avghist.div(len(spiketrains)*binwidth*1e-3)
    return avghist


def nrn_avg_rate_adaptive(spiketrains, tstart, tstop, binwidth=10.0, minsum=15):
    """
    Compute running average firing rate of population.

    @param spiketrains  list of Hoc.Vector() objects or other iterables,

    @param minsum       minimum number of spikes required in adaptive
                        window used to compute per-bin firing rate

    @return             h.Vector() meanfreq, computed as follows:

                        For bin i, the corresponding mean frequency f_mean[i] is determined by centering an adaptive square window on i and widening the window until the number of spikes under the window equals <minsum>. Then f_mean[i] is calculated as:

                            f_mean[i] = N[i] / (m * binwidth * trials)

                        where m is the number of
                        bins included after widening the adaptive window.

    @see                https://neuron.yale.edu/neuron/static/new_doc/programming/math/vector.html#Vector.psth
    """
    # first compute the histogram
    psth = nrn_sum_psth(spiketrains, tstart, tstop, binwidth)

    # convert to firing rates
    ntrials = len(spiketrains)
    minsum = min(minsum, psth.sum())
    vmeanfreq = h.Vector()
    vmeanfreq.psth(psth, binwidth, ntrials, minsum)
    return vmeanfreq


def composite_spiketrain(spike_times, duration, dt_out=1.0, select=None):
    """
    Construct composite spiketrain for population.
    Baseline is subtracted.

    @param  spike_times : list(numpy.array)
            List of spike time vectors.

    @param  duration : float
            Total duration of simulation (ms).

    @param  select : int / enumerable[int]
            Number of spike trains to select or their indices.
    """
    num_samples = int(np.round(duration, 3)) / dt_out
    # Select subset of spike trains
    if select is None:
        selected = range(len(spike_times))
    elif isinstance(select, int):
        selected = np.random.choice(len(spike_times), select, replace=False)
    else:
        selected = select # assume list/tuple/numpy array
    # Insert ones at spike indices
    binary_trains = np.zeros((num_samples, len(selected)))
    for i_select, i_cell in enumerate(selected):
        st = spike_times[i_cell]
        spike_indices = (st / dt_out).astype(int)
        binary_trains[spike_indices, i_select] = 1.0
    return binary_trains.mean(axis=1)


def composite_spiketrain_coherence(trains_a, trains_b, duration, freq_res=1.0, overlap=0.75):
    """
    Coherence between composite (pooled) spike trains as binary signals, averaged.

    See Terry and Griffin 2008 : https://doi.org/10.1016/j.jneumeth.2007.09.014

    Notes
    -----

    - coherence is the cross-spectral density, normalized by the individual
      spectral densities
    - the cross-spectral density is the Fourrier transform of the cross-correlation
    """
    Ts = 1.0 # ms
    nperseg = int(1e3 / Ts / freq_res) # fs / dF
    noverlap = int(nperseg * overlap)
    comp_trains = [composite_spiketrain(trains, duration, dt_out=Ts)
                    for trains in trains_a, trains_b]
    f, Cxy = scipy.signal.coherence(*comp_trains, fs=1e3/Ts, nperseg=nperseg, noverlap=noverlap)
    return f, Cxy


def combinatorial_spiketrain_coherence(
        spike_trains, duration,
        freq_res=1.0, overlap=0.75, max_comb=100):
    """
    Spiketrain coherence as average over all combinations of two groups.

    As described in McManus et al 2016, https://doi.org/10.1152/jn.00097.2016

    WARNING: scales badly with population size, e.g. for a population of 50
    spike trains, there are 126e12 ways of dividing them in two groups.
    """
    import scipy.special
    import itertools

    Ts = 1.0 # ms
    nperseg = int(1e3 / Ts / freq_res) # fs / dF
    noverlap = int(nperseg * overlap)

    N = len(spike_trains)
    k = N / 2
    inds_N = set(range(N))
    ncomb_tot = scipy.special.comb(N, k) # binomial coefficient
    if ncomb_tot > max_comb:
        grp_a_inds = [set(np.random.choice(N, k, replace=False)) for i in range(max_comb)]
    else:
        # Generate all possible combinations of choosing k out of N
        grp_a_inds = [set(inds) for inds in itertools.combinations(inds_N, k)]
    Cxy_sum = 0.0
    for grp_a in grp_a_inds:
        grp_b = inds_N - grp_a
        comp_spk_a = composite_spiketrain(spike_trains, dt_out=Ts, select=grp_a)
        comp_spk_b = composite_spiketrain(spike_trains, dt_out=Ts, select=grp_b)
        f, Cxy = scipy.signal.coherence(comp_spk_a, comp_spk_b,
                                        fs=1e3/Ts, nperseg=nperseg, noverlap=noverlap)
        Cxy_sum += Cxy
    return f, Cxy_sum / len(grp_a_inds)


def morgera_covariance_complexity(
    v_rec, t_rec, interval,
    t_window=1000.0, t_overlap=800.0):
    """
    Morgera covariance complexity ('synchronization index').

    Works on the raw signal so also measures synchronization of sub-threshold
    voltage oscillations. I.e. not only spike synchronization.

    Explained succinctly in this review article: https://arxiv.org/abs/q-bio/0603035

    Arguments
    ---------

    @param  v_rec : numpy.array
            Array with signals as columns (along first axis).

    References
    ----------

    S. S. Morgera, "Information Theoretic Complexity and Relation to Pattern
    Recognition", IEEE Transactions on Systems, Man, and Cybernetics,
    15 (1985) 608-619
    """
    M_values = []
    delta_t = t_window - t_overlap
    Ts = t_rec[1] - t_rec[0]
    t_start, t_stop = interval
    t0_values = np.arange(t_start, t_stop, delta_t)
    t1_values = []
    for t0 in t0_values:
        window = [t0, t0+t_window]
        if window[1] > t_stop:
            break
        irange = [int((t - t_rec[0])/Ts) for t in window]
        islice = np.s_[irange[0]:irange[1]] # slice object

        # SVD and singular values
        u, s, vh = np.linalg.svd(v_rec[islice, :], full_matrices=True)
        lambas = s**2
        sigmas = lambas / lambas.sum()
        C = - 1./np.log(len(sigmas)) * np.sum(sigmas * np.log(sigmas))

        t1_values.append(window[1])
        M_values.append(1 - C)

    return zip(t0_values, t1_values), M_values


def get_efel_features(
        vm, t, interval, features,
        threshold=None, raise_warnings=True, **kwargs):
    """
    Calculate electrophysiology features using eFEL.

    https://efel.readthedocs.io/en/latest/eFeatures.html
    http://bluebrain.github.io/eFEL/efeature-documentation.pdf

    @param  features : list[str]
            List of eFeature names described in eFEL documentation.

    @param  kwargs : dict[str, double/int]
            Extra eFEL parameters like 'interp_step' or required parameters
            listed under a feature.

    Example
    -------

    >>> from bgcellmodels.common import signal
    >>> features = ['ISI_values', 'burst_ISI_indices', 'burst_mean_freq', 'burst_number']
    >>> feat_vals = signal.get_efel_features(vm, t, [0.5e3, 3e3], features)
    >>> print("\n".join(("{} : {}".format(name, val) for name,val in feat_vals)))

    """

    import efel
    efel.reset()

    if threshold is not None:
        efel.setThreshold(threshold)

    if kwargs is not None:
        for pname, pval in kwargs.iteritems():
            if isinstance(pval, float):
                efel.setDoubleSetting(pname, pval)
            elif isinstance(pval, int):
                efel.setIntSetting(pname, pval)
            else:
                raise ValueError("Unknown eFEL parameter or unexpected type:"
                                 "{} : {}".format(pname, pval))

    # Put trace in eFEL compatible format
    efel_trace = {
        'T': t,
        'V': vm,
        'stim_start': [interval[0]],
        'stim_end': [interval[1]],
    }

    # Calculate spike times from response
    values = efel.getFeatureValues(
        [efel_trace],
        features,
        raise_warnings=raise_warnings
    )
    efeat_values = {
        feat_name: values[0][feat_name] for feat_name in features
    }

    return efeat_values


def get_all_pyelectro_features(
        vm, t, interval, analysis_params=None):
    """
    Calculate electrophysiology features using PyElectro.

    https://pyelectro.readthedocs.io/en/latest/pyelectro.html

    @param  analysis_params : dict[str, float]
            Parameters for calculation of voltage trace metrics.

    @return features : dict[str, float]
            Dictionary containing all PyElectro features.

    Notes
    -----

    See IClampAnalysis.analyse() method for calculating individual features:
    https://github.com/lkoelman/pyelectro/blob/master/pyelectro/analysis.py#L1204
    """

    from pyelectro import analysis

    if analysis_params is None:
        analysis_params = {
            'peak_delta':       1e-4, # the value by which a peak or trough has to exceed its neighbours to be considered outside of the noise
            'baseline':         -10., # voltage at which AP width is measured
            'dvdt_threshold':   0, # used in PPTD method described by Van Geit 2007
        }

    trace_analysis = analysis.IClampAnalysis(
                            vm,
                            t,
                            analysis_params,
                            start_analysis = interval[0],
                            end_analysis = interval[1],
                            show_smoothed_data = False
                        )

    trace_analysis.analyse()

    return trace_analysis.analysis_results


def compute_PRC(t_spikes, t_stim, sort_by='phi',
                exclude_multi_stim_isis=True):
    """
    Compute phase response curve (PRC) using traditional method.

    I.e. compute the phases of incoming spikes (t_stim) in the spike period
    (phi) and the phase advance/delay of each recorded spike (delta_phi).

    @author     Lucas Koelman, ported from Matlab implementation by M. Giugliano
                and J. Couto at https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=155735&file=/pkj_prc/matlab

    @param      t_spikes  :  1 x N array
                times of the recorded spikes

    @param      t_stim : 1 x M array
                times of the delivered pulses

    @param      sort_by : str
                Order of data points: 'spike_times' or 'phi'

    @return     (phi, delta_phi) : tuple[<1 x M array>, <1 x M array>]
                Phases of t_stim and phase shifts of spike_times. Phase
                shifts are negative for delays and positive for advances.
    """
    # Make sure datatype is numpy array
    t_spikes    = np.asarray(t_spikes)
    t_stim      = np.asarray(t_stim)

    mean_ISI    = np.mean(np.diff(t_spikes))
    M           = len(t_stim)
    K           = len(t_spikes)
    phi         = np.zeros(M)
    delta_phi   = np.zeros(M)

    excluded_stim_mask = np.zeros(M, dtype=bool) # np.array([False] * M)

    for i, t_pulse in enumerate(t_stim):

        # Find preceding and following spike
        idx_preceding = np.where(t_spikes < t_pulse)[0][-1]
        idx_following = idx_preceding + 1
        if idx_following >= K:
            break # pulse after last spike : no ISI

        # Check if there is more than one pulse in ISI
        isi_pulses_mask = ((t_stim > t_spikes[idx_preceding]) &
                           (t_stim < t_spikes[idx_following]))

        if isi_pulses_mask.astype(int).sum() > 1 and exclude_multi_stim_isis:
            excluded_stim_mask[i] = True
            continue
            # TODO: remove p

        tau = t_pulse - t_spikes[idx_preceding]
        Ti = t_spikes[idx_following] - t_spikes[idx_preceding]

        phi[i] = tau / mean_ISI
        delta_phi[i] = 1.0 - (Ti / mean_ISI)

    if sort_by == 'phi':
        phi       = phi[~excluded_stim_mask]
        delta_phi = delta_phi[~excluded_stim_mask]
        isort     = np.argsort(phi)
        phi       = phi[isort]
        delta_phi = delta_phi[isort]

    elif sort_by == 'spike_times':
        phi[excluded_stim_mask]         = np.nan
        delta_phi[excluded_stim_mask]   = np.nan

    else:
        raise ValueError(sort_by)

    return phi, delta_phi


def compute_PRC_corrected(t_spikes, t_stim):
    """
    Compute phase response curve (PRC) using corrected method
    (Phoka et al., 2010).

    The corrected method solves the dead zone ('bermuda triangle') problem in
    the PRC where it is impossible to get a phase advance (delta phi > 0) for a
    stimulus arriving late in the phase (phi ~= 1). See Phoka et al. (2010),
    Fig. 2.N vs 3.E.

    @author     Lucas Koelman, ported from Matlab implementation by M. Giugliano
                and J. Couto at https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=155735&file=/pkj_prc/matlab

    @param      t_spikes  :  1 x N array
                times of the recorded spikes

    @param      t_stim : 1 x M array
                times of the delivered pulses

    @return     (phi, delta_phi) : two numpy.array[float] of shape (M*3,)
                Phases of t_stim and phase shifts of spike_times. Phase
                shifts are negative for delays and positive for advances.
    """
    mean_ISI    = np.mean(np.diff(t_spikes))
    M           = len(t_stim)
    K           = len(t_spikes)
    phi         = np.zeros((M, 3))
    delta_phi   = np.zeros((M, 3))


    for i, t_pulse in enumerate(t_stim):
        idx_preceding = np.where(t_spikes < t_pulse)[0][-1]
        idx_following = idx_preceding + 1
        if idx_following >= K:
            break # pulse after last spike : no ISI

        tau     = t_pulse - t_spikes[idx_preceding]
        Ti      = t_spikes[idx_following] - t_spikes[idx_preceding]
        Tim1    = t_spikes[idx_preceding] - t_spikes[idx_preceding-1]
        Tip1    = t_spikes[idx_following+1] - t_spikes[idx_following]

        phi[i, 0]  = tau / mean_ISI
        phi[i, 1] = (Tim1 + tau) / mean_ISI
        phi[i, 2] = (tau - Ti) / mean_ISI

        delta_phi[i, 0] = 1.0 - (Ti / mean_ISI)
        delta_phi[i, 1] = 1.0 - (Tim1 / mean_ISI)
        delta_phi[i, 2] = 1.0 - (Tip1 / mean_ISI)


    phi         = phi.reshape((-1,))
    delta_phi   = delta_phi.reshape((-1,))
    return phi, delta_phi