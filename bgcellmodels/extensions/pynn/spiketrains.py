"""
Spike train generation for PyNN

Notes
-----

The 'spike_times' argument of a PyNN SpikeSourceArray accepts
- a list of Sequence or numpy array
- a function returning such a list
"""

import numpy as np
from pyNN.parameters import Sequence
from bgcellmodels.common import spikelib


def continuous_bursts(bursting_fraction, synchronous, rng,
                      T_burst, dur_burst, f_intra, f_inter,
                      f_background, duration, min_spk_dt=1.0):
    """
    Make generator for continuous regularly bursting spike trains.
    """
    if synchronous:
        make_bursts = spikelib.make_oscillatory_bursts
    else:
        make_bursts = spikelib.make_variable_bursts

    def spike_seq_gen(cell_indices):
        """
        Spike sequence generator

        @param  cell_indices : list(int)
                Local indices of cells (NOT index in entire population)

        @return spiketimes_for_cell : list(Sequence)
                Sequencie of spike times for each cell index
        """
        # Choose cell indices that will emit bursty spike trains
        num_bursting = int(bursting_fraction * len(cell_indices))
        bursting_cells = rng.choice(cell_indices, 
                                    num_bursting, replace=False)

        spiketimes_for_index = []
        for i in cell_indices:
            if i in bursting_cells:
                # Spiketimes for bursting cells
                burst_gen = make_bursts(T_burst, dur_burst, f_intra, f_inter,
                                        rng=rng, max_dur=duration)
                # Get rid of spikes that are spaced too closely
                spk_cell = np.fromiter(burst_gen, float)
                spk_clean = spk_cell[:-1][~(np.diff(spk_cell) < min_spk_dt)]
                spiketimes = Sequence(spk_clean)
            else:
                # Spiketimes for background activity
                number = int(2 * duration * f_background / 1e3)
                if number == 0:
                    spiketimes = Sequence([])
                else:
                    spk_cell = np.add.accumulate(
                        rng.exponential(1e3/f_background, size=number))
                    spk_clean = spk_cell[:-1][~(np.diff(spk_cell) < min_spk_dt)]
                    spiketimes = Sequence(spk_clean)
            spiketimes_for_index.append(spiketimes)
        return spiketimes_for_index

    return spike_seq_gen


def synchronous_bursts_during(intervals, bursting_fraction, 
                          T_burst, dur_burst, f_intra, f_inter, f_background, 
                          duration, randomize_bursting, rng):
    """
    Generator for a given fraction of spiketrains firing synchronized bursts
    during given time intervals.

    @param      intervals : iterable(tuple[float, float])
                Sequence of time intervals in which cells should burst.

    @return     function(iterable(index)) -> iterable(Sequence)
                Function that returns a sequence of spike times for each
                cell index.
    """

    def spike_seq_gen(cell_indices):
        """
        Spike sequence generator
        """
        # First pick bursting cells during each bursty interval
        num_bursting = int(bursting_fraction * len(cell_indices))
        num_intervals = len(intervals)
        if randomize_bursting:
            # pick new bursting cells in each interval
            bursting_ids = [rng.choice(cell_indices, num_bursting, replace=False) 
                                for i in range(num_intervals)]
        else:
            bursting_ids = [rng.choice(cell_indices, num_bursting, replace=False)] * num_intervals

        spiketimes_for_index = []
        for i in cell_indices:
            # Get intervals where this cell is bursting
            burst_intervals = [intervals[j] for j in range(num_intervals) if i in bursting_ids[j]]
            
            if len(burst_intervals) > 0:
                # Spiketimes for bursting cells
                spikegen = spikelib.generate_bursts_during(burst_intervals, 
                                T_burst, dur_burst, f_intra, f_inter, 
                                f_background, duration, max_overshoot=0.25, rng=rng)
                spiketimes = Sequence(np.fromiter(spikegen, float))
            else:
                # Spiketimes for background activity
                number = int(2 * duration * f_background / 1e3)
                if number == 0:
                    spiketimes = Sequence([])
                else:
                    spiketimes = Sequence(np.add.accumulate(
                        rng.exponential(1e3/f_background, size=number)))
            spiketimes_for_index.append(spiketimes)
        return spiketimes_for_index
    return spike_seq_gen


def synchronous_permuted_bursts(
        T_burst         = 20.0,
        phi_burst       = 0.0,
        num_spk_burst   = 4,
        f_intra         = 180.0,
        f_background    = 5.0,
        max_dt_spk      = 1.0,
        t_refrac_pre    = 5.0,
        t_refrac_post   = 10.0,
        bursting_fraction = 0.25,
        intervals       = None,
        duration        = 10e3,
        rng             = None,
        min_spk_dt      = 4.0):
    """
    Make spiketrains that burst semi-synchronously in each cycle of a
    reference sinusoid. In each cycle only a given fraction of spiketrains has
    a burst and the indices of bursting spiketrains are permuted in each cycle.

    Arguments
    ---------

    @param      T_burst : float
                Burst period (ms)

    @param      phi_burst : float
                Phase angle (degrees) where burst starts

    @param      num_spk_burst : tuple[int, int]
                Number of spikes per burst, lower and upper bound

    @param      f_intra : float
                Intra-burst firing rate (Hz)

    @param      f_background : float
                Firing rate when cell is not bursting.

    @param      rng : numpy.Random
                Random generator (optional)

    @param      intervals : iterable(tuple[float, float])
                Sequence of time intervals in which cells should burst.

    @param      max_dt_spk : float
                Max dt (ms) added to spiketimes inside a burst. The interval
                between spikes in a burst is 1/f_intra + uniformly sampled
                value in (0, max_dt_spk). I.e. f_intra is the max firing rate.

    @return     function(iterable(index)) -> iterable(Sequence)
                Function that returns a sequence of spike times for each
                cell index.
    """
    # Make a mask that has the bursting cell indices for each cycle in its rows
    # i.e. element mask[i][j] is the j-th bursting cell index in cycle i
    def spike_seq_gen(cell_indices):
        """
        Spike sequence generator
        """
        # First pick bursting cells during each bursty interval
        num_bursting = int(bursting_fraction * len(cell_indices))
        num_cycles = int(duration / T_burst) + 1
        # element i are cell ids that burst during cycle i
        cycle_burst_ids = np.array(
            [rng.choice(cell_indices, num_bursting, replace=False) 
                for i in range(num_cycles)])

        T_intra = 1e3 / f_intra
        spk_burst_centered = np.arange(0, num_spk_burst[1]*T_intra, T_intra)

        # For each spiketrains
        # - Generate bursty spikes in cycles that it's active and overlap with interval
        # - Generate background spikes in remaining cycles.
        spiketimes_for_index = []
        for cell_idx in cell_indices:
            # Spiketimes for background activity
            number = int(2 * duration * f_background / 1e3)
            spk_bg = np.add.accumulate(
                rng.exponential(1e3/f_background, size=number))
            bg_del_mask = spk_bg > duration # np.zeros_like(spk_bg, dtype=bool)

            # If cell never bursts, skip
            bursting_cycles = np.where(cycle_burst_ids == cell_idx)[0]
            if len(bursting_cycles) == 0:
                spiketimes_for_index.append(Sequence(spk_bg))
                continue

            # Add bursty spikes in burst cycles
            all_spk = []
            for i_cycle in bursting_cycles:
                t0_cycle = i_cycle * T_burst
                t1_cycle = t0_cycle + T_burst
                t0_burst = t0_cycle + (phi_burst / 360.0 * T_burst)
                if intervals is not None:
                    # only add burst in cycle if cycle falls in bursty interval
                    if not any((((ival[0] <= t0_cycle) and (t1_cycle <= ival[1])) for ival in intervals)):
                        continue
                num_spk = rng.randint(num_spk_burst[0], num_spk_burst[1]+1)
                spk_var_dt = rng.uniform(0.0, max_dt_spk, num_spk)
                spk_burst = t0_burst + spk_burst_centered[0:num_spk] + spk_var_dt
                all_spk.append(spk_burst)

                # Delete background spikes in tspk0-refrac, tspk[-1]+refrac
                # - build a deletion mask and apply at end
                cycle_mask = ((spk_bg >= (spk_burst[0] - t_refrac_pre)) & 
                              (spk_bg <= (spk_burst[-1] + t_refrac_post)))
                bg_del_mask = bg_del_mask | cycle_mask

            # Delete background spikes that fall in refractory period around burst
            all_spk.append(spk_bg[~bg_del_mask])

            # Sort the spikes
            spk_merged = np.concatenate(all_spk)
            spk_merged.sort()

            # Get rid of spikes that are spaced too closely
            spk_clean = spk_merged[:-1][~(np.diff(spk_merged) < min_spk_dt)]

            spiketimes = Sequence(spk_clean)
            spiketimes_for_index.append(spiketimes)
        return spiketimes_for_index
    return spike_seq_gen


def synchronous_permuted_bursts_varfreq(
        T_burst         = 20.0,
        phi_burst       = 0.0,
        num_spk_burst   = (4, 6),
        f_intra         = 180.0,
        f_background    = 5.0,
        max_dt_spk      = 1.0,
        t_refrac_pre    = 5.0,
        t_refrac_post   = 10.0,
        bursting_fraction = 0.25,
        intervals       = None,
        duration        = 10e3,
        rng             = None,
        min_spk_dt      = 4.0):
    """
    Same as synchronous_permuted_bursts except each interval can have its own
    bursting parameters.

    Spiketrains that burst synchronously during each interval and fire at a
    background rate between those intervals. During the bursting intervals, 
    in each period of the bursting frequency only a given fraction of cells
    (spike trains) emit a burst. The cells that burst are selected randomly
    in each period.

    Arguments
    ---------

    @see    synchronous_permuted_bursts()
            
            Parameters that can be set per interval are:
            T_burst, num_spk_burst, bursting_fraction
    """
    # Preprocessing arguments : one for each interval
    if isinstance(T_burst, (float, int)):
        T_burst = [T_burst] * len(intervals)
    if isinstance(bursting_fraction, (float, int)):
        bursting_fraction = [bursting_fraction] * len(intervals)
    if len(num_spk_burst) != len(intervals):
        num_spk_burst = [num_spk_burst] * len(intervals)
    T_intra = 1e3 / f_intra

    # Make a mask that has the bursting cell indices for each cycle in its rows
    # i.e. element mask[i][j] is the j-th bursting cell index in cycle i
    def spike_seq_gen(cell_indices):
        """
        Spike sequence generator
        """
        # For each interval select cells that will burst in each period
        ival_burst_ids = []
        for i, interval in enumerate(intervals):
            num_cycles = int((interval[1] - interval[0]) / T_burst[i]) + 1
            num_bursting = int(bursting_fraction[i] * len(cell_indices))
            # matrix: row[i] is indices of bursting cells in period i
            cycle_burst_ids = np.array(
                [rng.choice(cell_indices, num_bursting, replace=False) 
                    for i in range(num_cycles)])
            ival_burst_ids.append(cycle_burst_ids)

        # For each spiketrain
        # - Generate bursty spikes in cycles that it's active and overlap with interval
        # - Generate background spikes in remaining cycles.
        spiketimes_for_index = []
        for cell_idx in cell_indices:
            # First generate all background spikes, then delete as necessary
            number = int(2 * duration * f_background / 1e3)
            spk_bg = np.add.accumulate(
                rng.exponential(1e3/f_background, size=number))

            # Mask which background spikes should be deleted, continued later
            bg_del_mask = spk_bg > duration

            # Create periodic bursts in each interval
            all_spk = []
            for i_ival, interval in enumerate(intervals):
                # Cycles in current interval where this cell emits a burst
                cell_bursting_cycles = np.where(ival_burst_ids[i_ival] == cell_idx)[0]
                if len(cell_bursting_cycles) == 0:
                    continue # cell doesn't burst during this interval

                # Add spikes in each bursting period
                for i_cycle in cell_bursting_cycles:
                    t0_cycle = interval[0] + i_cycle * T_burst[i_ival]
                    t0_burst = t0_cycle + (phi_burst / 360.0 * T_burst[i_ival])

                    num_spk = rng.randint(num_spk_burst[i_ival][0], 
                                          num_spk_burst[i_ival][1] + 1)
                    spk_var_dt = rng.uniform(0.0, max_dt_spk, num_spk)
                    spk_burst = t0_burst + np.arange(0, num_spk*T_intra, T_intra) + spk_var_dt
                    all_spk.append(spk_burst)

                    # Delete background spikes in tspk0-refrac, tspk[-1]+refrac
                    cycle_mask = ((spk_bg >= (spk_burst[0] - t_refrac_pre)) & 
                                  (spk_bg <= (spk_burst[-1] + t_refrac_post)))
                    bg_del_mask = bg_del_mask | cycle_mask

            # Delete background spikes that fall in refractory period around burst
            all_spk.append(spk_bg[~bg_del_mask])

            # Sort the spikes
            spk_merged = np.concatenate(all_spk)
            spk_merged.sort()

            # Get rid of spikes that are spaced too closely
            spk_clean = spk_merged[:-1][~(np.diff(spk_merged) < min_spk_dt)]

            spiketimes = Sequence(spk_clean)
            spiketimes_for_index.append(spiketimes)
        return spiketimes_for_index
    return spike_seq_gen

# only add burst in cycle if cycle falls in bursty interval
# i_cycle = next((i for i, ival in enumerate(intervals) if 
#                             ((ival[0] <= t0_cycle) and (t1_cycle <= ival[1]))), 
#                            None)

# Aliases for backward compatibility
make_bursty_spike_generator = continuous_bursts
bursty_spiketrains_during = synchronous_bursts_during


def test_spiketime_generator(gen_func, num_spiketrains, *args, **kwargs):
    """
    Test case for generate_modulated_spiketimes()
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    generator = gen_func(*args, **kwargs)
    cells_spiketimes = generator(np.arange(num_spiketrains))

    f_avg = 0.0

    for i, sequence in enumerate(cells_spiketimes):
        spiketimes = sequence.value
        print(min(np.diff(spiketimes)))

        y_vec = np.ones_like(spiketimes) * i
        ax.plot(spiketimes, y_vec, marker='|', linestyle='', color='red', snap=True)

        st_dt = np.diff(spiketimes)
        err_idx = np.where(st_dt < 1.0)[0]
        err_spk = spiketimes[err_idx]
        ax.plot(err_spk, np.ones_like(err_spk) * i, 'go', markerfacecolor=None)

        f_avg += len(spiketimes) / kwargs['duration'] * 1e3

    f_avg /= num_spiketrains
    print("Mean population spike rate is %f Hz" % f_avg)

    # Plot zero-phase times
    all_t_phi0 = []
    for i, interval in enumerate(kwargs['intervals']):
        T_burst = kwargs['T_burst']
        if not isinstance(T_burst, float):
            T_burst = T_burst[i]
        num_periods = int((interval[1] - interval[0]) // T_burst) + 1
        all_t_phi0.extend([interval[0] + j*T_burst for j in range(num_periods)])
    y0, y1 = ax.get_ylim()
    ax.vlines(all_t_phi0, y0, y1, color='black', linestyle='dashed', linewidths=0.5)

    ax.set_title("Result of '{}'".format(gen_func.__name__))
    ax.grid(True)
    plt.show(block=False)


if __name__ == '__main__':
    
    rng = np.random
    duration = 12e3
    T_burst = [75., 50., 30.]
    intervals = [(2e3, 4e3), (6e3, 8e3), (10e3, 12e3)]
    phi_burst = 0.0
    f_intra = 200.
    f_background = 10.
    num_spiketrains = 1000
    num_spk_burst = (3, 5)
    max_dt_spk = 1.0
    refrac = 20.0, 10.0
    frac_bursting = 0.50

    # Test for synchronous_permuted_burst
    T_burst = 16.666666666666668 # 18.51851851851852
    f_background = 2.0
    frac_bursting = 0.0641025641025641
    intervals = [(1.0, 12e3)]
    test_spiketime_generator(
        # First two arguments are for test functions
        synchronous_permuted_bursts, num_spiketrains, 
        # Generator function args and kwargs:
        T_burst=T_burst, phi_burst=phi_burst, num_spk_burst=num_spk_burst, 
        f_intra=f_intra, f_background=f_background, max_dt_spk=max_dt_spk,
        t_refrac_pre=refrac[0], t_refrac_post=refrac[1], 
        bursting_fraction=frac_bursting, intervals=intervals, 
        duration=duration, rng=np.random)

    # Test for 'bursty_permuted_spiketrains'
    # T_burst = [75., 50., 30.]
    # test_spiketime_generator(
    #     # First two arguments are for test functions
    #     synchronous_permuted_bursts_varfreq, num_spiketrains, 
    #     # Generator function args and kwargs:
    #     T_burst=T_burst, phi_burst=phi_burst, num_spk_burst=num_spk_burst, 
    #     f_intra=f_intra, f_background=f_background, max_dt_spk=max_dt_spk,
    #     t_refrac_pre=refrac[0], t_refrac_post=refrac[1], 
    #     bursting_fraction=frac_bursting, intervals=intervals, 
    #     duration=duration, rng=np.random)
