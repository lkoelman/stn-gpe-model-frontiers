"""
Population signal analysis tools.

@author     Lucas Koelman
@date       17/07/2019

For use with Jupyter notebooks, making use of the following globar variables:

@global     pops_segments : dict[str, neo.Segment]
            Mapping of recorded population labels to data segments


JUPYTER EXAMPLE
---------------

>>> %load_ext autoreload
>>> %autoreload 1
>>> %aimport bgcellmodels.common.popsignal
>>> popsig = bgcellmodels.common.popsignal

"""

import os, re
import numpy as np
import scipy.signal
import matplotlib
import matplotlib.pyplot as plt
import elephant
# import numba

from bgcellmodels.common import analysis
from bgcellmodels.extensions.neo import signal as neoutil
from bgcellmodels.common.config_global import analysis_data as _data


################################################################################
# Plotting tools
################################################################################

def save_figure(fname, fig=None, **kwargs):
    """
    Save given or current figure.
    For LaTeX embedding, use extension pdf/pgf/eps in the figure name.

    kwargs: see https://matplotlib.org/api/_as_gen/matplotlib.pyplot.savefig.html
    """
    kwargs.setdefault('bbox_inches', 'tight') # prevents cropping
    kwargs.setdefault('transparent', True) # transparent background, see also 'frameon'
    try:
        fname += '_sweep-val-{}'.format(sweep_var_value)
    except NameError:
        pass
    fname += '.' + kwargs.setdefault('format', 'pdf')
    fig_dir = kwargs.get('dir', _data.save_fig_path)
    fig_filepath = os.path.join(fig_dir, fname)
    if fig is None:
        plt.savefig(fig_filepath, **kwargs) # save current figure
    else:
        fig.savefig(fig_filepath, **kwargs) # save specific figure

    print("Figure saved to file {}".format(fig_filepath))
    return fig_filepath


def offset_show_twin_yax(ax, offset=1.2):
    """
    Having been created by twinx, ax has its frame off, so the line of its
    detached spine is invisible.  First, activate the frame but make the patch
    and spines invisible. Then make the detached spine visible.
    """
    ax.spines["right"].set_position(("axes", offset))
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

    ax.spines["right"].set_visible(True)


def hide_axis(ax, frame=True, x=True, y=True):
    """
    Hide axis lines and frame box of figure.
    """
    ax.set_frame_on(not frame)
    ax.get_xaxis().set_visible(not x)
    ax.get_yaxis().set_visible(not y)


def set_axes_size(w, h, ax=None):
    """
    Set the size of the axes (plotting area) instead of the whole figure
    which is the default.

    w, h: width, height in inches
    """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)


def get_pop_color(pop_label):
    """ Colors from https://xkcd.com/color/rgb/ """
    if pop_label.startswith('CTX'):
        if 'axon' in pop_label.lower():
            return 'xkcd:lightblue'
        else:
            return 'c'

    elif pop_label.startswith('STR'):
        return 'm'

    elif pop_label.startswith('GPE'):
        if 'axon' in pop_label.lower():
            return 'xkcd:rose pink'
        else:
            return 'r'

    elif pop_label.startswith('STN'):
        if 'axon' in pop_label.lower():
            return 'xkcd:sage'
        else:
            return 'g'

    else:
        return 'b'

def get_pop_order(pop_label):
    if pop_label.startswith('CTX'):
        return 0
    elif pop_label.startswith('STR'):
        return 1
    elif pop_label.startswith('GPE'):
        return 2
    elif pop_label.startswith('STN'):
        return 3
    else:
        return 100

def plot_signal(signal, interval, channels=None, pop_label=None, title=None,
                ylabel=None, xlabel=None, ylim=None, figsize=None,
                export=False, **plot_kwargs):
    """
    Plot neo.AnalogSignal in given interval.
    """
    if figsize is None:
        figsize = (_data.fig_width, _data.fig_height)

    fig, ax = plt.subplots(figsize=figsize)
    if title == True:
        title = '{} ({})'.format(signal.name, pop_label)
    if title:
        fig.suptitle(title)
    if pop_label:
        sig_label = '{}-{}'.format(signal.name, pop_label)
    else:
        sig_label = signal.name

    islice = neoutil.make_slice(signal, interval)
    if isinstance(channels, slice):
        cslice = channels
    else:
        cslice = slice(channels)

    # sig_data = signal.asarray()
    # if sig_data.ndim == 2:
    #     sig_plotted = sig_data[islice, cslice]
    # else:
    #     sig_plotted = sig_data[islice]
    ax.plot(signal.times[islice], signal[islice, cslice],
            label=sig_label, **plot_kwargs)

    if ylim:
        ax.set_ylim(ylim)

    ax.set_xlim(interval)

    if ylabel is None:
        ylabel = '{} ({})'.format(signal.name, signal.units)

    ax.set_ylabel(ylabel)
    ax.set_xlabel('time ({})'.format(signal.times.units))
    ax.grid(True)

    fig.subplots_adjust(bottom=0.15) # prevent clipping of xlabel

    if _data.export_figs and export:
        fname = '{}_t-{:.1f}-{:.1f}'.format(
            sig_label, interval[0]*1e-3, interval[1]*1e-3)
        save_figure(fname, fig=fig)

    return fig, ax


################################################################################
# Continuous-time Signals
################################################################################


def plot_vm_signals(signal, cell_indices, interval, interval_only=True,
                    as_gids=False, title=None, trigger=None,
                    plot_labels=False, ylim=(-90, 25),
                    figsize=None, export=False, multiple_axes=True,
                    dy_major=50.0, dy_minor=10.0, **plot_kwargs):
    """
    Plot membrane voltage signals on vertically stacked axes.
    """
    if figsize is None:
        figsize = (_data.fig_width, _data.ax_height)

    if as_gids:
        pop_gids = list(signal.annotations['source_ids'])
        cell_indices = [pop_gids.index(gid) for gid in cell_indices]

    rec_dt = signal.sampling_period.magnitude
    tstart = signal.t_start.magnitude
    irange = [0, signal.shape[0]-1] if interval is None else [int((t-tstart)/rec_dt) for t in interval]
    times = signal.times[irange[0]:irange[1]]
    ydiff = ylim[1] - ylim[0]

    num_axes = len(cell_indices) if multiple_axes else 1
    fig, axes = plt.subplots(num_axes, 1,
                             figsize=(figsize[0], 0.5*figsize[1]*len(cell_indices)),
                             sharex=True, sharey=True)
    if title is None:
        title = "{} membrane voltage".format(signal.annotations['source_population'])
        fig.suptitle(title)

    if trigger is not None:
        if interval_only:
            trigger_times = trigger[(trigger >= interval[0]) & trigger <= interval[1]]
        else:
            trigger_times = trigger
            trigger_level = np.zeros_like(trigger_times) - 40.0

    # Plot each Vm on separate axis
    cell_gids = []
    trace_y_offsets = []
    trace_yticks_major = []
    trace_yticks_labels = []
    trace_yticks_minor = []

    for i_ax, i_cell in enumerate(cell_indices):
        try:
            ax = axes[i_ax]
        except TypeError:
            ax = axes

        if 'source_ids' in signal.annotations:
            gid = signal.annotations['source_ids'][i_cell]
            label = "id {}".format(gid)
            cell_gids.append(gid)
            # ax.text(.98, .05, signal.annotations['source_ids'][i_cell] , transform=ax.transAxes, ha='right')
        else:
            label = "cell {}".format(i_cell)

        # Offset for trace if plotted on single axis
        y0 = 0.0 if multiple_axes else i_ax * (ylim[1] - ylim[0])
        trace_y_offsets.append(y0)

        sigdata = signal.as_array()
        if interval_only:
            ax.plot(times, sigdata[irange[0]:irange[1], i_cell] + y0,
                    label=label, **plot_kwargs)
        else:
            ax.plot(signal.times, sigdata[:, i_cell] + y0, label=label,
                    **plot_kwargs)

        if trigger is not None:
            ax.plot(trigger_times, trigger_level + y0, marker='|', linestyle='',
                    snap=True, color='red', markersize=15)

        if multiple_axes:
            ax.grid(True, which='major')
            ax.set_ylabel(label)
            ax.set_yticks(np.arange(-100, 100, dy_minor), minor=True)
            ax.set_yticks(np.arange(-100, 100, dy_major), minor=False)
            ax.set_ylim(ylim)
            ax.set_xlim((times[0].magnitude, times[-1].magnitude))
        else:
            ax.text(int(times[0])-15, y0+ylim[0]+0.5*ydiff,
                    'id {}'.format(gid), ha='right', rotation=90)
            ticks = y0 + np.arange(ylim[0], ylim[1], dy_major)
            trace_yticks_major.extend(ticks)
            trace_yticks_labels.extend(ticks - i_ax * (ylim[1] - ylim[0]))
            # trace_yticks_minor.append(y0 + np.arange(ylim[0], ylim[1], 10))

    if not multiple_axes:
        ax.set_yticks(trace_yticks_major, minor=False)
        ax.set_yticklabels(['{:.0f}'.format(y) for y in trace_yticks_labels])
        ax.set_xlim((times[0].magnitude, times[-1].magnitude))
        ax.grid(True, which='major')

    #ax.set_ylabel("voltage ({})".format(signal.units))
    ax.set_xlabel('time (ms)')
    # fig.text(0.06, 0.5, "voltage ({})".format(signal.units), va='center', rotation='vertical')
    fig.subplots_adjust(bottom=0.15) # prevent clipping of xlabel

    # Save figure
    if _data.export_figs and export:
        pop_label = signal.annotations['source_population']
        fname = 'Vm-sigs_' + pop_label + '_gids-' + '-'.join([str(gid) for gid in cell_gids])
        save_figure(fname, fig=fig, bbox_inches='tight')


def plot_signal_interval(ax, signal, interval, channels, **kwargs):
    """
    Plot neo.AnalogSignal in interval on existing axis.
    """
    rec_dt = signal.sampling_period.magnitude
    tstart = signal.t_start.magnitude
    i_range = [int((t-tstart)/rec_dt) for t in interval]
    i_slice = slice(*i_range)
    ax.plot(signal.times[i_slice], signal[i_slice, channels], **kwargs)


################################################################################
# Spike and PSTH analysis
################################################################################

# @numba.jit("UniTuple(f8[:], 2)(f8[:],f8[:],f8,f8)", nopython=True)
def find_stimlocked_spikes(pulse_times, spike_times, window_lo, window_hi):
    """
    Find pulse-locked spikes in time window after each pulse.

    @param  window_lo : float
            Minimum time after pulse (ms) to look for phase-locked spike

    @param  window_hi : float
            Maximum time after pulse (ms) to look for phase-locked spike

    @return (indices, spike_times) : tuple(numpy.array[int], numpy.array[float])
            Indices of all pulses that have stimulus-spikes, and the spike
            times that are stimulus-locked
    """
    locked_indices = []
    locked_spike_times = []
    for i_pulse, t_pulse in enumerate(pulse_times):
        mask_following = (
            (spike_times > (t_pulse + window_lo)) &
            (spike_times <= (t_pulse + window_hi))
        )
        spikes_following = spike_times[mask_following]
        if spikes_following.size > 0:
            locked_indices.append(i_pulse)
            locked_spike_times.append(spikes_following)

    indices_array = np.array(locked_indices)
    times_array = np.array(locked_spike_times)
    return indices_array, times_array


# @numba.jit("f8[:](f8[:],f8[:],f8,f8)", nopython=True)
def find_stimlocked_indices(pulse_times, spike_times, window_lo, window_hi):
    """
    Find pulse-locked spikes in time window after each pulse.

    Set window_lo = 0.0 and window_hi = inter-pulse-interval to find all pulses
    that have stimulus-locked spikes.

    @param  window_lo : float
            Minimum time after pulse (ms) to look for phase-locked spike

    @param  window_hi : float
            Maximum time after pulse (ms) to look for phase-locked spike
    """
    locked_indices = []
    for i_pulse, t_pulse in enumerate(pulse_times):
        mask_following = (
            (spike_times > (t_pulse + window_lo)) &
            (spike_times <= (t_pulse + window_hi))
        )
        spikes_following = spike_times[mask_following]
        if spikes_following.size > 0:
            locked_indices.append(i_pulse)

    return np.array(locked_indices)


def calc_train_stimlock_fractions(spike_trains, onset_times, interval=None):
    """
    Calculate the percentage of pulses that have stimulus-locked spikes
    for each spike train in the population

    @param    spike_trains : list[Neo.SpikeTrain]
              Population spke trains
    """
    if interval is None:
        interval = _data.ROI_INTERVAL

    inter_pulse_interval = onset_times[1] - onset_times[0]
    mask = (onset_times >= interval[0]) & (onset_times <= interval[1])
    onset_times = onset_times[mask]
    trains_locked_fractions = []

    for st in spike_trains:
        pulse_indices = find_stimlocked_indices(onset_times, st.as_array(),
                                                0.0, inter_pulse_interval)
        locked_fraction = float(len(pulse_indices)) / len(onset_times)
        trains_locked_fractions.append(locked_fraction)

    return trains_locked_fractions


def plot_spiketrain(spike_trains, cell_indices, t_range, pop_label=None,
                    sharex=None, sharey=None, figsize=None, plot_compact=False,
                    export=False, order_by='cell_pop_idx', ticks_dx=1e3,
                    grid=True, y_labels='pop_index', sort_by='cell_gid'):
    """
    Plot spiketrains for one population.

    @param   cell_ids : list(int)
             Cell indices in population that will be visible (y-axis constrainment).

    @param   y_labels : str
            'pop_index' or 'cell_gid'

    @param   sort_by : str
             'pop_index' or 'cell_gid' : method for sorting spike trains
             before indexing them using <cell_indices>
    """
    if figsize is None:
        figsize = (_data.page_width, _data.ax_height)

    if sort_by == 'pop_index':
        spike_trains = sorted(spike_trains, key=lambda st: st.annotations['source_index'])
    elif sort_by == 'cell_gid':
        spike_trains = sorted(spike_trains, key=lambda st: st.annotations['source_id'])
    elif sort_by is not None:
        raise ValueError(sort_by)

    sim_dur = spike_trains[0].t_stop.magnitude

    # Don't plot all rastergrams in same figure
    fig_spikes = plt.figure(figsize=figsize)
    ax = plt.subplot(1,1,1, sharex=sharex, sharey=sharey)
    # fig_spikes, ax = plt.subplot(1, 1, figsize=(page_width,ax_height), sharex=sharex)
    if not plot_compact:
        fig_spikes.suptitle('{} spiketrains'.format(pop_label))

    # Plot all spiketrains but constrain y-axis later (so you can pan & zoom)
    # Only plot selected spike trains
    y_vals = np.empty(len(cell_indices))
    cell_gids = []
    for j, cell_index in enumerate(cell_indices):
        spiketrain = spike_trains[cell_index]

        if order_by == 'cell_pop_idx':
            y_vals[j] = cell_index
        elif order_by == 'caller':
            y_vals[j] = j

        if 'source_id' in spiketrain.annotations:
            cell_gids.append(spiketrain.annotations['source_id'])

        y_vec = np.ones_like(spiketrain) * y_vals[j]
        ax.plot(spiketrain, y_vec,
                marker='|', linestyle='',
                snap=True, color=get_pop_color(pop_label))


    # ax.set_xticks(np.arange(0, sim_dur+5000, 5000), minor=False) # uncomment for long time range
    if plot_compact:
        ax.set_yticks([]) # (np.arange(min(y_vals), max(y_vals)+1, 1), minor=False)
        ax.set_xticks(np.arange(0, sim_dur+1000, ticks_dx), minor=False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    else:
        ax.set_yticks(np.arange(min(y_vals), max(y_vals)+5, 5), minor=False)
        ax.set_xticks(np.arange(0, sim_dur+1000, ticks_dx), minor=False)
        ax.set_ylabel(y_labels)

    # Label each spike train
    if y_labels == 'cell_gid':
        if len(cell_gids) == 0 or len(cell_gids) != len(y_vals):
            print('WARNING: cell GID data not found. '
                  'Using population index for y-axis labels.')
        else:
            # ax.tick_params(axis='y', labelsize=14, labelcolor=grid_color)
            ax.set_yticks(y_vals, minor=False)
            ax.set_yticklabels([str(i) for i in cell_gids])
    elif y_labels != 'pop_index':
        raise ValueError(y_labels)

    ax.set_xlim(t_range)
    ax.set_ylim((min(y_vals)-0.5, max(y_vals)+0.5))
    ax.grid(grid, axis='x', which='major')

    if _data.export_figs and export:
        fname = 'rastergram_{}_cell-{}-{}'.format(pop_label, cell_indices[0], cell_indices[-1])
        save_figure(fname, fig=fig_spikes, bbox_inches='tight')

    return fig_spikes, ax

################################################################################
# Frequency analysis
################################################################################


def calc_psd(pop_label, interval=None, vm_sig=None, save=True):
    """
    Calculate PSD from membrane voltages of population.

    It computes the PSD for the individual Vm traces and then
    averages the resulting PSDs.
    """
    if vm_sig is None:
        segment = _data.pops_segments[pop_label]
        vm_sig = next((sig for sig in segment.analogsignals if sig.name == 'Vm'))

    # Ts = vm_sig.sampling_period.magnitude # ms
    fs = vm_sig.sampling_rate.rescale('Hz').magnitude

    islice = neoutil.make_slice(vm_sig, interval)
    vm_segment = vm_sig[islice, :]

    dF = max(0.5, fs / vm_segment.shape[0]) # dF for max window size is finest possible
    if dF != 0.5:
        print("Adjusted frequency resolution to data size: dF = {}".format(dF))

    # Compute PSD of all 100 Vm signals at the same time
    freqs, psd = elephant.spectral.welch_psd(vm_segment, freq_res=dF)

    # Average individual trace PSDs
    psd_avg = psd.sum(axis=0) / psd.shape[0]
    # psd_rel = psd_avg[0:int(250/dF)] # relevant region of psd

    # Save PSD
    if save:
        sig_label = pop_label + '_' + vm_sig.name
        _data.exported_data['PSD'][sig_label] = (freqs, psd_avg)

    # Save units for other plotting functions
    _data.psd_units = psd.units

    return freqs, psd_avg


def calc_spectrogram(pop_label, signal=None, max_avg=50, freq_res=1.0,
                     t_res=20.0, interval=None, save=True):
    """
    Calculate mean spectrogram (STFT) of membrane voltage traces.

    The resulting time axis is missing 'nperseg' values on each side since
    each PSD sample is calculated on the interval [-nperseg/2, +nperseg/2]
    around it.

    @param     max_avg : int
               Number of spectrograms to compute for averaging.
               WARNING: can lead to very high memory consumption.

    @return    freqs, t, Sxx
               Sxx has t-dimension along axis 0 and f-dimension along axis 1.
    """
    if signal is None:
        segment = _data.pops_segments[pop_label]
        signal = next((sig for sig in segment.analogsignals if sig.name == 'Vm'))

    if interval is None:
        interval = _data.ROI_INTERVAL

    sig_label = pop_label + '_' + signal.name

    # Plot spectrogram using STFT
    dt = signal.sampling_period.magnitude
    fs = 1/dt*1e3
    nperseg = int(fs/freq_res) # determines frequency resolution
    t_res = 20.0 # ms
    noverlap = nperseg - int(t_res/dt)
    islice = neoutil.make_slice(signal, interval)

    if signal.ndim == 2:
        # Spectrogram of all signals along axis 0
        vm_sig = signal[islice, :]
        if vm_sig.shape[1] > max_avg:
            sig_data = vm_sig.as_array()[:, 0:max_avg]
        else:
            sig_data = vm_sig.as_array()

        freqs, t, Sxx = scipy.signal.spectrogram(sig_data, 1/dt, axis=0,
                                                 window='hanning',
                                                 nperseg=nperseg,
                                                 noverlap=noverlap,
                                                 scaling='density')

        # Average over cell dimension -> Sxx dims are now (freqs, t)
        assert sig_data.shape[1] == Sxx.shape[1]
        Sxx = Sxx.mean(axis=1)
    else:
        assert signal.ndim == 1
        sig_data = signal.ravel()[islice]
        freqs, t, Sxx = scipy.signal.spectrogram(sig_data, 1/dt, window='hanning',
                                                 nperseg=nperseg, noverlap=noverlap, scaling='density')

    freqs = freqs * 1000
    t = t + signal.t_start.rescale('ms').magnitude

    # Save spectrogram
    if save:
        df = freqs[1]-freqs[0]
        _data.exported_data['spectrogram'][sig_label] = (
            freqs[0:int(50/df)], t, Sxx[:,0:int(50/df)])

    return freqs, t, Sxx


def integrate_subband(freqs, vals, fband, normalize=True):
    """
    Integrate spectral power in frequency band [f0, f1].
    """
    mask1 = (freqs >= fband[0]) & (freqs <= fband[1])
    mask2 = np.isclose(freqs, fband[0], atol=0.1) | np.isclose(freqs, fband[1], atol=0.1)
    fband_idx, = np.where(mask1 | mask2)
    return np.sum(vals[fband_idx]) / len(fband_idx)


################################################################################
# Phase analysis
################################################################################

def calc_mean_phase_vectors(spiketrains, pop_label, intervals=None,
                            analytic_signal=None,
                            reference_signal=None):
    """
    Calculate mean phase vector of spikes with reference to given analytic signal
    (e.g. BP filtered + Hilbert transformed).

    The mean phase vector for a cell (spike train) is obtained by looking up the
    complex value of the normalized analytic signal at each spike time,
    and averaging over all spikes in the spike train.

    The population phase vector is obtained by taking the mean complex value
    of all cell phase vectors in the population.

    Returns
    -------

    @return    mean_phase_vec : numpy.array(dtype=complex)
               Array containing mean phase vector of each spike train as rows.

    @return    pop_phase_vec : complex
               Mean phase vector over all cells in population.
    """
    # Gather spiketimes that fall within intervals of cortical beta bursts
    spikes_during = [] # list of numpy array
    for i, st in enumerate(spiketrains):
        spiketimes = st.magnitude

        mask = np.zeros_like(spiketimes, dtype=bool)
        for ival in intervals:
            mask = mask | ((spiketimes > ival[0]) & (spiketimes <= ival[1]))

        spikes_during.append(spiketimes[mask])


    Ts = reference_signal.sampling_period.magnitude
    t_start = reference_signal.t_start.magnitude
    mean_phase_vecs = []
    spike_phase_vecs = []
    for i, spiketimes in enumerate(spikes_during):
        analsig_indices = np.round((spiketimes-t_start)/ Ts).astype(int) # indices into analytic signal
        analsig_indices = analsig_indices[analsig_indices < analytic_signal.size]
        if analsig_indices.size > 0:
            cell_spike_phase_vecs = analytic_signal[analsig_indices]
            mean_phase_vecs.append(np.mean(cell_spike_phase_vecs)) # mean of normalized phase vectors
            spike_phase_vecs.append(cell_spike_phase_vecs)
        else:
            mean_phase_vecs.append(0.0 + 0.0j)

    # Save phase vectors for export
    _data.exported_data['cell_phase_vecs'][pop_label] = mean_phase_vecs = np.array(mean_phase_vecs)
    _data.exported_data['pop_phase_vecs'][pop_label] = pop_phase_vec = np.mean(mean_phase_vecs)
    _data.spike_phase_vectors[pop_label] = spike_phase_vecs

    return mean_phase_vecs, pop_phase_vec


def plot_phase_vectors(mean_phase_vecs, pop_phase_vec, pop_label, export=False,
                       rmax=None, rticks=None, rticks_pos=None, plot_cell_vecs=True,
                       cell_color='green', pop_color='red', cell_opacity=0.5,
                       mark_cells_rim=False, ref_vec=None, extra_pop_vecs=None,
                       extra_cell_vecs=None, extra_labels=None, extra_colors=None,
                       extra_cell_colors=None, pop_line_width=3, show_legend=True):
    """
    Plot mean phase vectors for individual neurons and whole population
    in Polar coordinates.
    """
    if extra_pop_vecs is None:
        extra_pop_vecs = []

    # Reference angle for all plotted vectors
    if ref_vec is None:
        ref_ang = 0.0
    else:
        ref_ang = np.angle(ref_vec)

    # Complex vectors to polar coordinates
    vec_angs = np.angle(mean_phase_vecs) - ref_ang
    vec_lens = np.abs(mean_phase_vecs)
    pop_ang = np.angle(pop_phase_vec) - ref_ang
    pop_mag = np.abs(pop_phase_vec)
    if rmax is None:
        vmax = pop_mag
        if extra_pop_vecs:
            vmax = max(vmax, *(np.abs(v) for v in extra_pop_vecs))

        rmax = (vmax // 0.1 + 1) * 0.1

    # if rticks is None:
    #     rticks = np.arange(0.1, 1.1, 0.1)

    fig = plt.figure()
    ax = plt.subplot(111, projection='polar')

    # Plot additional vectors if given
    for i, vec in enumerate(extra_pop_vecs):
        if plot_cell_vecs and extra_cell_vecs is not None:
            extra_angs = np.angle(extra_cell_vecs[i]) - ref_ang
            extra_lens = np.abs(extra_cell_vecs[i])
            ax.vlines(extra_angs, 0, extra_lens, color=extra_cell_colors[i],
                      alpha=cell_opacity, linewidth=1, snap=True)

        # Population vector after/over cell vectors
        ax.vlines(np.angle(vec) - ref_ang, 0, np.abs(vec), color=extra_colors[i],
                  linewidth=pop_line_width, label=extra_labels[i])

    # Plot cell vectors
    if plot_cell_vecs:
        ax.vlines(vec_angs, 0, vec_lens,
                  color=cell_color, alpha=cell_opacity, linewidth=1, snap=True)
        if mark_cells_rim:
            # Mark points on rim of polar plot
            ax.plot(vec_angs, np.zeros_like(vec_angs)+rmax, 'o',
                    color=cell_color, markersize=5)

    # Plot population vector as thick line
    # ax.plot(vec_angs, vec_lens, 'ro')
    ax.vlines(pop_ang, 0, pop_mag, label=pop_label,
              color=pop_color, linewidth=pop_line_width)

    # Format axes
    ax.grid(True)
    if rticks is not None:
        ax.set_rticks(rticks) # less radial ticks

    ax.set_rmax(rmax)
    if rticks_pos is not None:
        ax.set_rlabel_position(rticks_pos)  # Move radial labels away from plotted line

    ax.tick_params(axis='y', labelsize=14)  # labelcolor='blue'
    ax.tick_params(axis='x', labelsize=14)
    if show_legend:
        ax.legend(loc=(1, .8))

    # ax.set_title('Mean angle and vector length of {} neurons'.format(pop_label), va='bottom')

    # kw = dict(arrowstyle="->", color='g')
    # for angle, radius in zip(phases, magnitudes):
    #     ax.annotate("", xy=(angle, radius), xytext=(0, 0), arrowprops=kw)

    if _data.export_figs and export:
        fname = 'phase-vectors_{}'.format(pop_label)
        save_figure(fname, fig=fig)


def plot_phase_histogram(pop_label, ref_vec=None, num_bins=20,
                         face_alpha=0.1, bar_color='blue', export=False,
                         rlabel_angle='default', rmax=None, rticks=None,
                         rlabel_start=0.0, rmax_strict=False):
    """
    Plot mean phase vectors for individual neurons and whole population
    in Polar coordinates.
    """
    # Reference angle for all plotted vectors
    if ref_vec is None:
        ref_ang = 0.0
    else:
        ref_ang = np.angle(ref_vec)

    # Complex vectors to polar coordinates
    all_spike_vecs = np.concatenate(_data.spike_phase_vectors[pop_label], axis=0)
    assert (all_spike_vecs.ndim == 1) or (min(all_spike_vecs.shape) == 1)
    vec_angs = np.angle(all_spike_vecs) - ref_ang

    # Histogram of phase angles
    bin_counts, bin_edges = np.histogram(vec_angs, bins=num_bins, range=(-np.pi, np.pi))
    bin_fractions = bin_counts.astype(float) / np.sum(bin_counts)
    bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(num_bins)]

    fig = plt.figure()
    ax = plt.subplot(111, projection='polar')

    # Use custom colors and opacity
    bar_width = 2*np.pi / len(bin_counts)
    bars = ax.bar(bin_centers, bin_fractions, width=bar_width, bottom=0.0)
    r, g, b = matplotlib.colors.to_rgb(bar_color)
    for bar in bars:
        bar.set_facecolor((r, g, b, face_alpha))
        bar.set_edgecolor(bar_color)
        # bar.set_alpha(0.1)

    # Format axes
    if rticks is not None:
        ax.set_rticks(rticks)

    data_rmax = max(bin_fractions)
    if (rmax is not None) and (data_rmax <= rmax or rmax_strict):
        ax.set_rmax(rmax)

    grid_color = 'gray'
    ax.grid(True, color=grid_color)
    # Appearance of radial tick labels
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14, labelcolor=grid_color)  # labelcolor=grid_color

    if rlabel_start:
        ax.set_yticklabels(['{:.2f}'.format(y) if y >= rlabel_start else '' for y in ax.get_yticks()])
    else:
        ax.set_yticklabels(['{:.2f}'.format(y) for y in ax.get_yticks()])

    if rlabel_angle == 'minimum':
        label_pos = np.degrees(bin_centers[bin_counts.argmin()]) + 9
        ax.set_rlabel_position(label_pos)
    elif isinstance(rlabel_angle, (float, int)):
        ax.set_rlabel_position(rlabel_angle)

    if _data.export_figs and export:
        fname = 'phase-histogram_{}'.format(pop_label)
        save_figure(fname, fig=fig)


################################################################################
# Synaptic current analysis
################################################################################

def get_synapse_index(trace_name):
    matches = re.search(r'(?P<index>\d+)$', trace_name)
    syn_index = matches.group('index')
    assert syn_index is not None
    return int(syn_index)


def sorted_signals(segment, trace_name):
    trace_regex = trace_name + '(\d+)'
    return sorted([sig for sig in segment.analogsignals if re.search(trace_regex, sig.name)],
                  key=lambda sig: get_synapse_index(sig.name))


def plot_synapse_traces(pop_label, max_ncell, max_nsyn, interval, interval_only=True,
                        channel_gates=None, channel_currents=None, beta_phase=True,
                        trace_names=None, extra_sigs=None, extra_plot_with=None,
                        extra_max_plot=5, vm_plot_with=None, export=False):
    """
    Plot Synaptic dynamics for specific cells and their recorded synapses.

    Typically both the cells in the population and their synapses have
    been sampled.


    @pre    If arguments 'channel_gates' and 'channel_currents' are given,
            the signals they refer to are recorded from the same cells
            as the synaptic currents and conductances. I.e. the indices
            into channel-related signals are the same as those into
            the synapse-related signals.
    """
    if vm_plot_with is None:
        vm_plot_with = lambda tracename: tracename.startswith('i')

    segment = _data.pops_segments[pop_label]

    # NOTE: signals are akwardly ordered: one signal is the i-th synapse for all recorded cells
    default_tracenames = 'gAMPA', 'gNMDA', 'iGLU', 'gGABAA', 'gGABAB', 'iGABA'
    if trace_names is None:
        trace_names = default_tracenames

    all_synaptic_sigs = {tn: sorted_signals(segment, tn) for tn in trace_names}
    trace_groups = {k:'GABA' if ('GABA' in k) else 'GLU' for k in all_synaptic_sigs}
    signal = next((sigs[0] for sigs in all_synaptic_sigs.values() if len(sigs)>0), None)

    if signal is None:
        print("No synaptic traces for population {}".format(pop_label))
        return

    # index in all_synaptic_sigs of the synaptic traces that exist (are recorded)
    # existing_traces = [i for i,sigs in enumerate(all_synaptic_sigs.values()) if len(sigs)>0]
    existing_traces = [sig_name for sig_name, sigs in all_synaptic_sigs.items() if len(sigs)>0]
    num_ax_per_cell = len(existing_traces)
    num_cell = min(signal.shape[1], max_ncell)

    # Select synapses to plot.
    selected_synapses = {'GLU': [], 'GABA': []} # trace suffixes of selected synapses
    for group in selected_synapses.keys():
        tname = next((n for n in existing_traces if trace_groups[n] == group), None)
        if tname is None:
            continue
        for j_sig, sig in enumerate(all_synaptic_sigs[tname]):
            if j_sig >= max_nsyn:
                break
            selected_synapses[group].append(get_synapse_index(sig.name))

    # Get signal time data
    rec_dt = signal.sampling_period.magnitude
    tstart = signal.t_start.magnitude
    tstop = signal.t_stop.magnitude
    if interval is None:
        interval = (tstart, tstop)

    irange = [int((t-tstart)/rec_dt) for t in interval]
    times = signal.times[irange[0]:irange[1]]

    # Make the figure
    num_axes = num_cell * num_ax_per_cell
    fig, axes = plt.subplots(num_axes, 1,
                             figsize=(0.75*_data.page_width, num_axes*_data.ax_height),
                             sharex=True)
    # fig.suptitle("{} synapse dynamics".format(pop_label))

    # For each cell we try to plot all synapses of one type on one axis
    # TODO: for each plotted signal, check that the cell_gid is the same
    for i_cell in range(num_cell):
        # One axis for all iGLU, one for all iGABA, one for each conductance type.
        # This makes a maximum of 6 axes per cell
        for i_plotted, tracename in enumerate(existing_traces):
            i_ax = (i_cell * num_ax_per_cell) + i_plotted
            try:
                ax = axes[i_ax]
            except TypeError:
                ax = axes

            # Plot all synapses for this axis (same conductance or current)
            ax_l_sigs = [] # signals plotted on left axis
            for j_sig, sig in enumerate(all_synaptic_sigs[tracename]):
                ax_l_sigs.append(sig)
                i_syn = get_synapse_index(sig.name)
                if i_syn in selected_synapses[trace_groups[tracename]]:
                    label = None # if j_sig>0 else ax_tracename
                    if interval_only:
                        ax.plot(times, sig[irange[0]:irange[1], i_cell], label=label)
                    else:
                        ax.plot(signal.times, sig[:, i_cell], label=label)

            # Save axes limits for main traces
            ymin, ymax = ax.get_ylim()

            # Plot Beta trigger signal (zero phase)
            if beta_phase:
                ax.vlines(_data.phase_zero_times, ymin, ymax, label='$\phi$ = 0',
                          colors='black', linestyle='dashed', linewidths=0.5)

            # NOTE: cell index -> see recorder._get_current_segment() -> should save source_ids/channel_ids
            # TODO: plot spike times or Vm of source_indices in same plot
            ax_r = None

            # Plot additional signals on the right axis
            if channel_gates and tracename.startswith('g'):
                # Plot channel gating variables if we are plotting conductance
                ax_r = ax.twinx()
                gating_sigs = [sig for sig in segment.analogsignals if sig.name in channel_gates]
                for k, csig in enumerate(gating_sigs):
                    color, style = analysis.pick_line_options('red', 'broken', k)
                    plot_signal_interval(ax_r, csig, interval, i_cell, label=csig.name,
                                         color=color, linestyle=style)

                ax_r.legend(loc='upper right')
                ax_r.set_ylabel('open')

            # Plot extra traces specified by caller
            elif extra_sigs and extra_plot_with(tracename):
                ax_r = ax.twinx()
                rax_sigs = [sig for sig in segment.analogsignals if sig.name in extra_sigs]
                for k, csig in enumerate(rax_sigs):
                    if k+1 > extra_max_plot:
                        break
                    color, style = analysis.pick_line_options('red', 'broken', k)
                    plot_signal_interval(ax_r, csig, interval, i_cell, label=csig.name,
                                         color=color, linestyle=style)

                ax_r.legend(loc='upper right')

            # Plot channel currents if we are plotting synaptic currents
            elif channel_currents and tracename.startswith('i'):
                ax_r = ax.twinx()
                curr_sigs = [sig for sig in segment.analogsignals if sig.name in channel_currents]
                for k, csig in enumerate(curr_sigs):
                    color, style = analysis.pick_line_options('red', 'broken', k)
                    plot_signal_interval(ax_r, csig, interval, i_cell, label=csig.name,
                                         color=color, linestyle=style)

                ax_r.legend(loc='upper right')
                ax_r.set_ylabel('current ($mA/cm^2$)')

            # Plot membrane voltages recorded from given cell
            elif vm_plot_with(tracename) and 'source_indices' in ax_l_sigs[-1].annotations:
                src_idxs = ax_l_sigs[-1].annotations['source_indices']
                cell_pop_idx = src_idxs if isinstance(src_idxs, int) else src_idxs[i_cell]
                ax_r = ax.twinx()
                # Plot somatic voltage
                vsoma = next((sig for sig in segment.analogsignals if sig.name == 'Vm'))
                plot_signal_interval(ax_r, vsoma, interval, cell_pop_idx, label=vsoma.name,
                                     color='gray', linewidth=0.2)
                # Plot dendritic voltages
                # vm_sigs = [sig for sig in segment.analogsignals if sig.name.lower().startswith('v')]
                # for k, vsig in enumerate(vm_sigs):
                #     color, style = analysis.pick_line_options('red', 'broken', k)
                #     plot_signal_interval(ax_r, vsig, interval, i_cell, label=vsig.name,
                #                          color=color, linestyle=style)
                ax_r.legend(loc='upper right')
                ax_r.set_ylabel('Vm (mV)')

            # Annotation and axes
            ax.grid(True, axis='y')
            if tracename.startswith('i'):
                ax.set_ylabel('current (nA)')
            else:
                ax.set_ylabel('conductance (uS)')

            if i_plotted == 0:
                ax.set_title('{} cell {}'.format(pop_label, i_cell))

            # ax.set_xlabel('time (ms)')
            ax.set_ylim((ymin, ymax))
            ax.set_xlim((times[0].magnitude, times[-1].magnitude))
            ax.legend(loc='upper left')

    if _data.export_figs and export:
        fname = 'gsyn-{}_cells-{}_syns-{}'.format(pop_label, max_ncell,
                                                  max_nsyn)
        save_figure(fname, fig=fig)


def combine_current_signals(post_pop, afferents_currents, cell_gid=None, cell_pop_idx=None):
    """
    Get combined afferent current signal for a given presynaptic population
    and synaptic current types (e.g. AMPA, NMDA currents)

    @param    afferents_currents : dict[str, list(str)]
              Map of afferent population labels to recorded synaptic current names.
              e.g.: {'CTX': ('i_AMPA', 'i_NMDA')}

    @return   tuple(np.array, dict[str, np.array])
              Tuple containing:
              - np.array: sum total current as a function of time
              - dict[str, np.array]: sum total current per population
    """
    pops_itot = {}
    segment = _data.pops_segments[post_pop]

    for afferent_pop in afferents_currents.keys():
        # Current signals for afferent population, by current type (e.g. GABA-A/GABA-B)
        isyn_names = afferents_currents[afferent_pop]
        currents_traces = {curr: sorted_signals(segment, curr) for curr in isyn_names}

        # Sum current signals of same type
        currents_itot = {}
        for current, traces in currents_traces.items():
            # traces is all signals (synapses) recorded from the same synapse type
            currents_itot[current] = sum_cell_signals(traces, cell_gid=cell_gid,
                                                      cell_pop_idx=cell_pop_idx)

        # Sum current signals of different types
        pops_itot[afferent_pop] = sum(currents_itot.values())

    # Sum total current of all afferent populations
    allpops_itot = sum(pops_itot.values())
    return allpops_itot, pops_itot


def sum_cell_signals(traces=None, segment=None, sig_name=None,
                     cell_gid=None, cell_pop_idx=None, average=False):
    """
    Sum all signals recorded from the same cell.

    @param    traces : list(Neo.AnalogSignal)
              Set of signals that will be searched for traces
              recorded from the given cell.
    """
    if traces is None:
        traces = sorted_signals(segment, sig_name)

    itot = None
    for sig in traces:
        if cell_gid is not None:
            trace_idx = list(sig.annotations['source_ids']).index(cell_gid)
        elif cell_pop_idx is not None:
            trace_idx = list(sig.annotations['source_indices']).index(cell_pop_idx)

        isyn = sig.magnitude[:, trace_idx]
        itot = isyn if (itot is None) else (itot + isyn)

    if average:
        itot /= len(traces)

    print("Summed {} signals for cell gid {}".format(len(traces), cell_gid))
    return itot


def plot_oscillatory_traces(pop, exc_currents, inh_currents, ranking_slice=None,
                            interval=None, lpf_current=100, gating_signals=None,
                            trace_group_member=None, export_indices=None):
    """
    Plot factors contributing to network oscillations & synchronization.

    @param    trace_group_member : str
              Name of any recorded signal in the same recorded trace group as
              plotted variabed. Used for finding cell indices.
    """
    if interval is None:
        interval = _data.ROI_INTERVAL

    currents_by_action = {'EXC': exc_currents, 'INH': inh_currents}
    segment = _data.pops_segments[pop]

    # Get cells with recorded synapses and rank by phase locking strength
    cell_pop_idx_ranked = list(_data.exported_data['phaselock_ranking_source_indices'][pop])
    test_signal = next((sig for sig in segment.analogsignals if
                        sig.name.startswith(trace_group_member)), None)
    if test_signal is None:
        print("No synaptic traces for population {}".format(pop))
        return
    recorded_cell_idx = list(test_signal.annotations['source_indices'])
    print("Recorded cell indices: {}".format(recorded_cell_idx))
    print("Phase-locking ranked indices: {}".format(cell_pop_idx_ranked))
    recorded_cell_idx_ranked = sorted(recorded_cell_idx, key=lambda i: cell_pop_idx_ranked.index(i))
    plotted_cell_idx = recorded_cell_idx_ranked[ranking_slice]
    print("Ranked recorded cells: {}".format(recorded_cell_idx_ranked))

    # Get signal time data
    rec_dt = test_signal.sampling_period.magnitude
    tstart = test_signal.t_start.magnitude
    tstop = test_signal.t_stop.magnitude
    t0, t1 = interval
    irange = [int((t-tstart)/rec_dt) for t in interval]
    times = test_signal.times[irange[0]:irange[1]]

    # Phase signal
    phase_ref = _data.sigmean_bpss
    phase_slice = neoutil.make_slice(_data.sigmean_bpss, interval)
    t_zcross = _data.phase_zero_times
    phase_zcross = t_zcross[(t_zcross > interval[0]) & (t_zcross < interval[1])]

    # Voltage signal
    vm_sig = next((sig for sig in _data.pops_segments[pop].analogsignals if sig.name == 'Vm'))
    vm_slice = neoutil.make_slice(vm_sig, interval)

    # Current signal
    isig_ref = sorted_signals(_data.pops_segments[pop], exc_currents.values()[0][0])[0]
    isig_slice = neoutil.make_slice(isig_ref, interval)
    isig_times = isig_ref.times[isig_slice]

    # Gating variables
    if isinstance(gating_signals, (list, tuple)):
        # Ensure format is {<name of gating variable>: <list of signal labels>}
        gating_signals = {sig_name: [sig_name] for sig_name in gating_signals}

    # Filter design
    Fs = isig_ref.sampling_rate.rescale('Hz').magnitude
    Fn = Fs / 2. # Nyquist frequency
    hpfreq, lpfreq = 1.0, lpf_current
    low, high = hpfreq / Fn, lpfreq / Fn
    sos = scipy.signal.butter(4, high, btype='lowpass', analog=False, output='sos')

    # Make the figure
    num_cell = len(plotted_cell_idx)
    num_ax_per_cell = 2 + int(len(gating_signals.keys()) > 0)
    num_axes = num_cell * num_ax_per_cell

    def add_rastergram(ax):
        """
        Add rastergram of afferents on new background axis
        """
        axb = ax.twinx()
        axb.set_ylim((0, 1))
        hide_axis(axb)
        for EI_label, pops_currents in currents_by_action.items():
            base_color = 'green' if EI_label.startswith('E') else 'red'
            y_mid = 0.15 if EI_label.startswith('E') else 0.05
            rastergram_presynaptic_pops_pooled(axb, post_gid, pop, pops_currents.keys(), interval, base_color, y_mid)

    # For each selected cell, plot currents and presynaptic rates
    for i_cell, cell_idx in enumerate(plotted_cell_idx):

        fig, axes = plt.subplots(num_ax_per_cell, 1,
                                 figsize=(0.75*_data.page_width, num_ax_per_cell*_data.ax_height),
                                 sharex=False)
        ax_i_offset = 0 # i_cell * num_ax_per_cell

        # FIXME: comment temp solution
        # post_gid = network_params[pop]['gids'][cell_idx]
        vm_trace_idx = list(vm_sig.annotations['source_indices']).index(cell_idx)
        post_gid = vm_sig.annotations['source_ids'][vm_trace_idx]
        cell_ranking = cell_pop_idx_ranked.index(cell_idx)

        #########################################################
        # FIGURE AXIS 1: membrane voltage
        ax1a = axes[ax_i_offset + 0]
        ax1a.set_title('Cell index {} - gid {} - rank {}'.format(cell_idx, post_gid, cell_ranking))
        ## voltage on left axis
        ax1a.plot(vm_sig.times[vm_slice], vm_sig[vm_slice, vm_trace_idx], label='V_{m}')
        ax1a.set_ylim((-100, 25)) # create empty space above trace
        ax1a.set_ylabel('voltage (mV)')

        ## Rastergram of afferent spikes
        add_rastergram(ax1a)

        ## hilbert phase
        # ax2 = ax1.twinx()
        # ax2.plot(phase_ref.times[phase_slice], analphase_sigmean[phase_slice],
        #          label='$\phi$', color='magenta', alpha=0.5)
        # ax2.set_ylim((-4*np.pi, np.pi)) # above voltage traces
        # ax2.set_ylabel('phase (rad)')
        plot_phase_grid(phase_zcross, ax1a, False, False)

        #########################################################
        # FIGURE AXIS 2: synaptic currents
        ## combined current on left axis
        ax2a = axes[ax_i_offset + 1]
        # another current on third axis
        # ax2c = ax2a.twinx()

        # For each group of afferents (exc & inh pops), plot combined current and spikes
        for EI_label, pops_currents in currents_by_action.items():
            # Sum excitatory/inhibitory current signals
            itot_sig, itot_bypop = combine_current_signals(pop, pops_currents, cell_gid=post_gid)
            itot_filt = scipy.signal.sosfiltfilt(sos, itot_sig) # Filter current signal

            EXCITATORY = EI_label.startswith('E')
            axi = ax2a # if EXCITATORY else ax2c
            base_color = 'green' if EXCITATORY else 'red'
            flip = -1.0 if EXCITATORY else 1.0

            # Plot current
            base_color = 'green' if EI_label.startswith('E') else 'red'
            axi.plot(isig_times, flip*itot_filt[isig_slice], label='$i_{{{0:}}}$'.format(EI_label),
                     color=base_color, alpha=1.0)
            axi.set_ylabel('current (nA)') # color=base_color
            # axi.tick_params(axis='y', colors=base_color, size=4, width=1.5)

            # Make space below
            y0, y1 = axi.get_ylim()
            y0 = 0 # y0-0.2*(y1-y0)
            axi.set_ylim(y0, y1)

        ax2a.legend()

        ## Rastergram of afferent spikes
        add_rastergram(ax2a)

        ## Phase as gridlines
        plot_phase_grid(phase_zcross, ax2a, True, True)

        #########################################################
        # FIGURE AXIS 3: gating variables
        if gating_signals:
            ## gating variables on left axis
            ax3a = axes[ax_i_offset + 2]
            ## rastergram on background
            ax3b = ax2a.twinx()
            ax3b.set_ylim((0, 1))
            hide_axis(ax3b)

            # plot mean of gating variable recorded from distinct compartments
            for gate_name, sig_names in gating_signals.items():
                traces = [sig for sig in segment.analogsignals if any(
                    sig.name.startswith(n) for n in sig_names)]
                if len(traces) == 0:
                    continue
                msig_slice = neoutil.make_slice(traces[0], interval)
                msig_times = traces[0].times[msig_slice]
                m_avg = sum_cell_signals(traces, cell_pop_idx=cell_idx, average=True)
                ax3a.plot(msig_times, m_avg[msig_slice],
                          label="{} (mean)".format(gate_name))

            # plot gating variables recorded from each compartment
            # for sig in gating_sigs:
            #     trace_idx = list(sig.annotations['source_indices']).index(cell_idx)
            #     ax3a.plot(msig_times, sig[msig_slice, trace_idx], label=sig.name)

            # Legend and axes
            ax3a.legend()
            plot_phase_grid(phase_zcross, ax3a, True, True)

        # Cleanup axes
        for ax in fig.axes:
            ax.set_xlim(interval)

        if _data.export_figs and export_indices and cell_idx in export_indices:
            fname = 'phaselocking_cell-idx{}-gid{}-rank{}'.format(cell_idx, post_gid, cell_ranking)
            fpath = save_figure(fname, fig=fig)
            print('Saved figure as ' + fpath)


def get_presynaptic_gids(post_gid, pre_pop, post_pop):
    """
    Get presynaptic cell GIDs for given post-synaptic GID.
    """
    return [
        i for i,j in _data.network_params[pre_pop][post_pop]['conpair_gids']
            if j==post_gid
    ]

def get_presynaptic_pop_indices(post_index, pre_pop, post_pop):
    """
    Get presynaptic cell index in its population based on post-synaptic cell
    index in its own population.
    """
    return [
        i for i,j in _data.network_params[pre_pop][post_pop]['conpair_pop_indices']
            if j==post_index
    ]


def rastergram_presynaptic_pops_pooled(
        ax, post_gid, post_pop, pre_pops, interval, color, y_mid):
    """
    Add rastergram of POOLED presynaptic spikes (one population per line).
    """
    # Presynaptic spike trains from gid
    pre_spikes = []
    for pre_pop in pre_pops:
        pre_gids = get_presynaptic_gids(post_gid, pre_pop, post_pop)
        pre_spikes.extend(_data.spiketrains_by_gid(pre_gids, pre_pop.split('.')[0]))

    # Combined rastergram for presynaptic cells
    t0, t1 = interval
    all_spiketimes = [t.magnitude for t in pre_spikes]
    tot_spiketrain = np.concatenate([t[(t>t0) & (t<t1)] for t in all_spiketimes])
    y_vec = np.ones_like(tot_spiketrain) * y_mid
    ax.plot(tot_spiketrain, y_vec, marker='|', linestyle='', snap=True, color=color)


def plot_phase_grid(zero_crossings, ax, set_xticks=False, label_xticks=False):
    """
    Plot phase zero crossings as vertical grid lines on axis.
    """
    y0, y1 = ax.get_ylim()
    ax.vlines(zero_crossings, y0, y1, label='$\phi$=0',
              colors='black', linestyle='dashed', linewidths=0.5)
    if set_xticks:
        ax.set_xticks(zero_crossings)
    if set_xticks and label_xticks:
        ax.set_xticklabels(['{}$\pi$'.format(2*i) for i in range(len(zero_crossings))])


def sum_total_current(currents_traces, cell_idx, interval=None):
    """
    Sum synaptic currents for all recorded synapses of a given cell.

    @param    currents_traces : dict[str, list(Neo.AnalogSignal)]
              Map of trace names (without index) to signals


    """
    itot = 0.0
    itot_tracetype = {k: 0.0 for k in currents_traces.keys()}
    for current, traces in currents_traces.items():
        # traces is all signals (synapses) recorded from the same synapse type
        for sig in traces:
            if interval is not None:
                islice = neoutil.make_slice(sig, interval)
            else:
                islice = slice(None) # same as np.s_[:]

            isyn = sig.magnitude[islice, cell_idx]
            itot_tracetype[current] += isyn.sum()

        itot += itot_tracetype[current]

    return itot, itot_tracetype


def calc_exc_inh_ratio(pop_label, exc_currents, inh_currents, rec_ids,
                       interval=None, print_report=True, conductance=False):
    """
    Get total synaptic current for all afferent population and the ratio
    of excitatory to inhibitory currents.

    @param    rec_ids : list(int)
              which cell to use out of all recorded cells

    @param    conductance : bool
              signals are conductances instead of currents
    """
    segment = _data.pops_segments[pop_label]

    # Signal units
    test_current = exc_currents.values()[0][0]
    test_signals = sorted_signals(segment, test_current)
    rec_units = test_signals[0].sampling_period
    rec_dt = rec_units.magnitude # ms
    assert rec_units.dimensionality.string == 'ms'
    if interval is None:
        interval = [test_signals[0].t_start.magnitude, test_signals[0].t_stop.magnitude]

    delta_t = (interval[1] - interval[0])

    # Report intermediate calculations
    reports = []

    def afferents_total_current(afferents_currents, cell_idx, reports=reports):
        """
        Get total current for all given afferents onto cell.

        @param    afferents_currents : dict[str, list(str)]
                  Map of afferent population labels to recorded synaptic current names
        """
        itot = 0.0 # total current for all given afferents
        pops_currents_itot = {}

        for afferent_pop in afferents_currents.keys():
            # Afferent currents, by current type
            isyn_names = afferents_currents[afferent_pop]
            aff_traces = {curr: sorted_signals(segment, curr) for curr in isyn_names}
            num_rec_aff = len(aff_traces.values()[0]) # number of synapses for afferent population
            itot_sum, itot_bytrace = sum_total_current(aff_traces, cell_idx,
                                                       interval=interval)
            pops_currents_itot[afferent_pop] = itot_bytrace

            # Sum number of afferents (matrix column), should be same for each cell
            num_syn_aff = sum(
                _data.network_params[afferent_pop][pop_label]['conn_matrix'][:, 0] > 0)

            reports.append("\t- {} afferents ({}): {}/{} synapses recorded".format(
                afferent_pop, ",".join(isyn_names), num_rec_aff, num_syn_aff))

            # Multiply total currents by ratio of recorded to afferent synapses
            # FIXME: if separate connection matrix for surrogate/real cells but same synaptic trace name
            #        this is not a good estimate
            itot += itot_sum * num_syn_aff / num_rec_aff
        return itot, pops_currents_itot

    cells_i_info = []
    for rec_id in rec_ids:
        reports.append("{} cell {}:".format(pop_label, rec_id))
        itot_exc, iexc_bytype = afferents_total_current(exc_currents, rec_id)
        itot_inh, iinh_bytype = afferents_total_current(inh_currents, rec_id)

        # excitatory currents are negative by convention
        if not conductance:
            itot_exc *= -1.0

        if itot_inh == 0:
            ratio = np.inf
        else:
            ratio = itot_exc / itot_inh

        # Save all current info per recorded cell
        # units for synaptic currents: nA = nC/s = pC/ms
        iaff_info = {}
        iaff_info.update(iexc_bytype)
        iaff_info.update(iinh_bytype)
        iaff_info['EXC:INH'] = ratio
        iaff_info['EXC'] = itot_exc
        iaff_info['INH'] = itot_inh
        if conductance:
            iaff_info['EXC_integral_aoc'] = itot_exc * rec_dt # pC/ms * ms = pC
            iaff_info['INH_integral_aoc'] = itot_inh * rec_dt
            iaff_info['EXC_integral_avg'] = itot_exc * rec_dt / delta_t
            iaff_info['INH_integral_avg'] = itot_inh * rec_dt / delta_t
        else:
            iaff_info['EXC_charge_pC'] = itot_exc * rec_dt # pC/ms * ms = pC
            iaff_info['INH_charge_pC'] = itot_inh * rec_dt
            iaff_info['EXC_avg_nA'] = itot_exc * rec_dt / delta_t
            iaff_info['INH_avg_nA'] = itot_inh * rec_dt / delta_t
        cells_i_info.append(iaff_info)

    # Population average ratio EXC/INH
    cells_i_ratio = [info['EXC:INH'] for info in cells_i_info]
    ratio = sum(cells_i_ratio) / len(cells_i_ratio)

    # Save data
    if conductance:
        key_pop = "G_exc_inh_ratio"
        key_cells = "G_afferents"
    else:
        key_pop = "I_exc_inh_ratio"
        key_cells = "I_afferents"

    _data.exported_data[key_pop][pop_label] = ratio
    _data.exported_data[key_cells][pop_label] = cells_i_info

    if conductance:
        reports.append("""
{}: TOTAL conductance estimate
    G_EXC_avg = {} uS
    G_INH_avg = {} uS
""".format(pop_label,
           [i['EXC_integral_avg'] for i in cells_i_info],
           [i['INH_integral_avg'] for i in cells_i_info]))

        reports.append("""
{}: Ratio of integrated conductance EXC / INH:
    => cell ratios = {}
    => pop average ratio = {}""".format(pop_label, cells_i_ratio, ratio))

    else:
        reports.append("""
{}: TOTAL currents estimate
    I_EXC_avg = {} nA
    I_INH_avg = {} nA
""".format(pop_label, [i['EXC_avg_nA'] for i in cells_i_info],
           [i['INH_avg_nA'] for i in cells_i_info]))

        reports.append("""
{}: Ratio of integrated current EXC / INH:
    => cell ratios = {}
    => pop average ratio = {}""".format(pop_label, cells_i_ratio, ratio))

    report = "\n".join(reports)
    if print_report:
        print(report)

    return ratio, cells_i_info, report


################################################################################
# Sweep analysis
################################################################################

# Compare computed metrics over multiple simulations

def gather_metrics(all_cell_metrics, metric):
    """
    Gather burst metric from all cells into a list
    """
    if isinstance(all_cell_metrics[0][metric], (float, int)):
        map_func = lambda x,y: x + [y]
    elif isinstance(all_cell_metrics[0][metric], list):
        map_func = lambda x,y: x + y
    else:
        map_func = lambda x,y: x + list(y)

    all_cell_vals = reduce(map_func,
                           (cell_metrics[metric] for cell_metrics in all_cell_metrics), [])
    return all_cell_vals


def plot_metric_sweep(pop_labels, metric_name=None, metric_func=None,
                      independent='sweep', pop_colors=None, regress='linear',
                      regression_report=False, metric_kwargs={}, plot_kwargs={},
                      export=False, **kwargs):
    """
    Linear regression of metric

    @param    independent : str
              Independent variable for plot x-axis and linear regression

    @param    kwargs
              title, color, x_label, y_label, x_ticks, y_lim, ...
    """
    if regression_report:
        from statsmodels.formula.api import ols
        import pandas

    fig, ax = plt.subplots(figsize=(_data.fig_width, _data.fig_height))

    # Function for converting saved metric to y-value
    if metric_func is None:
        metric_func = lambda x: x

    if isinstance(pop_labels, str):
        pop_labels = [pop_labels]

    if pop_colors is None:
        pop_cmap = get_pop_color
    elif isinstance(pop_colors, (list, tuple)):
        pop_cmap = lambda pop: pop_colors[pop_labels.index(pop)]
    else:
        # assume function
        pop_cmap = pop_colors

    # Set independent variable for regression
    sweep_vals = np.array(sorted(_data.analysis_results.keys()))
    if independent == 'sweep':
        x = sweep_vals
        x_label = _data.sweep_var_legend
        x_ticks = sweep_vals
    elif independent == 'custom':
        x = kwargs['x']
        x_label = kwargs['x_label']
        x_ticks = kwargs.get('x_ticks', None)
    else:
        raise ValueError("Argument <independent> must be either 'sweep', 'exc_inh_ratio', or 'custom'.")


    ax.set_title(kwargs.get('title', metric_name))

    # Plot for each population
    for i_pop, pop_label in enumerate(pop_labels):
        sig_label = pop_label

        # Function for calculating metric
        y_vals = []
        for i, sweep_value in enumerate(sweep_vals):
            metric_vals = _data.analysis_results[sweep_value][metric_name][sig_label]
            pop_metric_kwargs = {
                k: metric_kwargs[k][i_pop] for k in metric_kwargs.keys()
            }
            y_vals.append(metric_func(metric_vals, **pop_metric_kwargs))

        # Linear regression
        y = np.array(y_vals)
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)


        # Scatter plot + linear fit
        color = pop_cmap(pop_label)
        plt_kwargs = {
            'color': color,
            'linestyle': '--', 'linewidth': 1,
            'marker': 'o', 'markersize': 6
        }
        plt_kwargs.update(plot_kwargs)
        ax.plot(x, y, label=pop_label, **plt_kwargs)

        if regress == 'linear':
            rcolor = 'k' # color
            ax.plot(x, intercept + slope*x, '-', color=rcolor, lw=1, label='linear fit')

            # Print regression statistics
            ax.text(.98, .05 + i_pop * .20,
                    '$R^2$ = {:.2f}\n$p^r$ = {:f}'.format(r_value**2, p_value),
                    color=rcolor, transform=ax.transAxes, ha='right')

            if regression_report:
                data = pandas.DataFrame({'x': x, metric_name: y})
                model = ols("{} ~ {}".format(metric_name, 'x'), data).fit()
                print(model.summary()) # summary2()

    if x_ticks is not None:
        ax.set_xticks(x_ticks)
    if 'y_lim' in kwargs:
        ax.set_ylim(kwargs['y_lim'])

    ax.set_ylabel(kwargs.get('y_label', metric_name))
    ax.set_xlabel(x_label)
    ax.grid(True)
    ax.legend()

    fig.subplots_adjust(bottom=0.15) # prevent clipped xlabel

    # set_axes_size(ax_width, ax_height, ax)
    if _data.export_figs and export:
        fname = ax.get_title().replace(' ', '_') + '_pop-' + '-'.join(pop_labels)
        save_figure(fname, fig=fig)


def boxplots_burst_metrics(pop_label, metric_names=None, export_metrics=None):
    """
    Compare burst metrics for same population across sweep values
    """
    sweep_vals = np.array(sorted(_data.analysis_results.keys()))
    if metric_names == None:
        metric_names = _data.analysis_results[sweep_vals[0]]['burst_metrics'][pop_label][0].keys()

    sweep_metrics = {}
    for metric_name in metric_names:
        for i, sweep_value in enumerate(sweep_vals):
            cell_metrics = _data.analysis_results[sweep_value]['burst_metrics'][pop_label]
            all_cell_vals = gather_metrics(cell_metrics, metric_name)
            sweep_metrics.setdefault(metric_name, []).append(all_cell_vals)

        # Plot boxplots
        fig, ax = plt.subplots(figsize=(_data.ax_width, _data.ax_height))
        ax.set_title('{} {}'.format(pop_label, metric_name))
        bp = ax.boxplot(sweep_metrics[metric_name], 0, 'g+')
        ax.set_xticklabels(sweep_vals)
        ax.set_ylim((0, ax.get_ylim()[1]))
        ax.set_ylabel(metric_name)
        ax.set_xlabel(_data.sweep_var_legend)
        ax.grid(True, which='major', axis='y')
        fig.subplots_adjust(bottom=0.15) # prevent clipped xlabel

        # Save figure
        if (_data.export_figs and (export_metrics is not None)
            and (metric_name in export_metrics)):
            fname = '{}_{}'.format(metric_name, pop_label)
            save_figure(fname, fig=fig)


def regression_burst_metrics(pop_label, metric_names=None, detailed=True, export_metrics=None):
    """
    Linear regression of burst metrics vs sweep variable.
    """
    from statsmodels.formula.api import ols
    import pandas

    sweep_vals = np.array(sorted(_data.analysis_results.keys()))
    if metric_names == None:
        metric_names = _data.analysis_results[sweep_vals[0]]['burst_metrics'][pop_label][0].keys()

    sweep_metrics = {}
    for metric_name in metric_names:
        for i, sweep_value in enumerate(sweep_vals):
            cell_metrics = _data.analysis_results[sweep_value]['burst_metrics'][pop_label]
            all_cell_vals = gather_metrics(cell_metrics, metric_name)
            sweep_metrics.setdefault(metric_name, []).append(all_cell_vals)

        # Linear regression
        x = sweep_vals
        y = [np.median(mvals) for mvals in sweep_metrics[metric_name]]
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

        # Detailed linear regression
        if detailed:
            data = pandas.DataFrame({'x': x, metric_name: y})
            model = ols("{} ~ {}".format(metric_name, 'x'), data).fit()
            print(model.summary()) # summary2()

        fig, ax = plt.subplots(figsize=(_data.ax_width, _data.ax_height))
        ax.set_title('{} {}'.format(pop_label, metric_name))
        ax.set_xticks(sweep_vals)
        ax.plot(x, y, 'o', color='g', label='original data')
        ax.plot(x, intercept + slope*x, 'k--', label='linear fit')
        ax.set_ylabel(metric_name)
        ax.set_xlabel(_data.sweep_var_legend)
        ax.grid(True)
        ax.text(.98, .05, '$R^2$ = {:.2f}\n$p^r$ = {:f}'.format(r_value**2, p_value),
                transform=ax.transAxes, ha='right')
        fig.subplots_adjust(bottom=0.15) # prevent clipped xlabel

        # Save figure
        if (_data.export_figs and (export_metrics is not None)
            and (metric_name in export_metrics)):
            fname = '{}_{}'.format(metric_name, pop_label)
            save_figure(fname, fig=fig)


def summarize_burst_metrics(pop_label, metric_names, axsize=None, export=False):
    """
    Plot multiple burst metrics (median across cells) in single figure with different y-axes.
    """
    sweep_vals = _data.sweep_vals
    if axsize is None:
        axsize = (_data.ax_width, _data.ax_height)


    fig, ax1 = plt.subplots()
    axes = [ax1] + [ax1.twinx() for m in metric_names[1:]]
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    sweep_metrics = {}
    for i_metric, metric_name in enumerate(metric_names):
        for i, sweep_value in enumerate(sweep_vals):
            cell_metrics = _data.analysis_results[sweep_value]['burst_metrics'][pop_label]
            all_cell_vals = gather_metrics(cell_metrics, metric_name)
            sweep_metrics.setdefault(metric_name, []).append(all_cell_vals)

        # Linear regression
        x = sweep_vals
        y = [np.median(mvals) for mvals in sweep_metrics[metric_name]]
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

        # Plot metric
        ax = axes[i_metric]
        act, = ax.plot(x, y, '--', marker='.', linewidth=1, markersize=6,
                       color=color_cycle[i_metric], label='original data')
        # ax.plot(x, intercept + slope*x, 'k--', label='linear fit')

        ax.set_ylabel(metric_name, color=act.get_color())
        ax.tick_params(axis='y', colors=act.get_color(), size=4, width=1.5)
        # ax.text(.98, .05, '$R^2$ = {:.2f}\n$p^r$ = {:f}'.format(r_value**2, p_value),
        #         transform=ax.transAxes, ha='right')

        if i_metric == 0:
            ax.set_xticks(sweep_vals)
            ax.set_xlabel(_data.sweep_var_legend)
            ax.grid(True)
        elif i_metric == 2:
            offset_show_twin_yax(ax, offset=1.12)

    fig.subplots_adjust(bottom=0.15) # prevent clipped xlabel
    set_axes_size(axsize[0], axsize[1], ax1)

    # Save figure
    if _data.export_figs and export:
        fname = 'burst-metrics_{}_{}'.format(_data.sweep_var_name, pop_label)
        save_figure(fname, fig=fig)


def compare_covariance_complexity(sig_label, f_band=None, detailed=False, ymax=None):
    """
    Compare distribution of Morgera-index values for time segments
    in simulation.
    """
    from statsmodels.formula.api import ols
    import pandas

    sweep_vals = np.array(sorted(_data.analysis_results.keys()))
    M_datasets = []
    M_intervals = []
    for i, sweep_value in enumerate(sweep_vals):
        fbands_M = _data.analysis_results[sweep_value]['Morgera_index'][sig_label]
        if isinstance(fbands_M, dict):
            # Morgera index was computed for several frequency bands
            if f_band:
                t, M = fbands_M[f_band]
            else:
                t, M = fbands_M[(5, 200)]
        else:
            # Morgera index was only computed in default band (5, 200) Hz
            if f_band:
                raise ValueError('Morgera index was only computed in f = (5, 200)')

            t, M = fbands_M

        M_datasets.append(M)
        M_intervals.append(t)

    M_array = np.array(M_datasets)
    mean_M = np.mean(M_array, axis=1)
    std_M = np.std(M_array, axis=1)
    med_M = np.median(M_array, axis=1)
    ival_width = M_intervals[0][0][1] - M_intervals[0][0][0]

    # Plot boxplots
    fig, ax = plt.subplots(figsize=(_data.ax_width, _data.ax_height))
    ax.set_title('Covariance complexity over {} ms intervals'.format(ival_width))
    bp = ax.boxplot(M_datasets, 0, 'g+')
    ax.set_xticklabels(sweep_vals)
    # ax.set_ylim((0, 1))
    ax.set_ylabel('M (0-1)')
    ax.set_xlabel(_data.sweep_var_legend)
    ax.grid(True, which='major', axis='y')
    fig.subplots_adjust(bottom=0.15) # prevent clipped xlabel

    # Linear regression
    metric_name = 'M'
    x = sweep_vals
    y = med_M
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

    # Detailed linear regression
    if detailed:
        data = pandas.DataFrame({'x': x, metric_name: y})
        model = ols("{} ~ {}".format(metric_name, 'x'), data).fit()
        print(model.summary()) # summary2()

    # Plot linear regression
    fig, ax = plt.subplots(figsize=(_data.ax_width, _data.ax_height))
    ax.set_title('Median covariance complexity over {} ms intervals'.format(ival_width))

    ax.plot(x, y, 'o', color='g', label='original data')
    ax.plot(x, intercept + slope*x, 'k--', label='linear fit')
    if ymax:
        ax.set_ylim((0, ymax))

    ax.set_xticks(sweep_vals)
    ax.set_ylabel('M (0-1)')
    ax.set_xlabel(_data.sweep_var_legend)
    ax.grid(True)
    ax.text(.98, .05, '$R^2$ = {:.2f}\n$p^r$ = {:f}'.format(r_value**2, p_value),
            transform=ax.transAxes, ha='right')
    fig.subplots_adjust(bottom=0.15) # prevent clipped xlabel
