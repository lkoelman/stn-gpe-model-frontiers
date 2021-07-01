"""
Functions for setting up neuronal stimulation.

@author     Lucas Koelman

@date       26-04-2019
"""
from __future__ import division

from scipy import signal
import numpy as np


def make_pulse_train(
        frequency=250.0,
        pulse_width_ms=0.300,
        amp0=1.0,
        amp1=0.0, 
        duration=10e3,
        dt=0.01,
        phase_deg=0.0,
        coincident_discontinuities=False,
        off_intervals=None):
    """
    Make monphasic or biphasic pulse train.

    @param  amp0 : float
            Amplitude of first phase (positive pulse)

    @param  amp1 : float
            Amplitude of second phase (negative pulse).
            Set to zero for monophasic pulse train.

    @param  duration : float
            Duration (ms)
    
    @param  dt : float
            Sampling period (ms)

    @param  coincident_discontinuities
            If true, time vector will have adjacent elements with the same time
            value at points of d
    """
    
    freq_ms = frequency * 1e-3
    duty_cycle = pulse_width_ms * freq_ms # duty cycle = PW / period
    time_axis = np.arange(0, duration + dt, dt)
    omega = 2 * np.pi * freq_ms
    phase_rad = phase_deg * np.pi / 180.0

    if pulse_width_ms  >= 1 / freq_ms:
        raise ValueError('Pulse width must be smaller than period')
    
    # square wave between [+1, -1] (2 peak-to-peak, baseline 0)
    dbs_vec = signal.square(omega * time_axis + phase_rad, duty_cycle)
    pos_phase = (dbs_vec == 1)
    neg_phase = (dbs_vec == -1)
    dbs_vec[pos_phase] = amp0
    dbs_vec[neg_phase] = amp1

    # turn DBS off in silent intervals
    if off_intervals is not None:
        silent_mask = np.zeros_like(time_axis, dtype=bool)
        for (t_off, t_on) in off_intervals:
            silent_mask = silent_mask | ((time_axis >= t_off) & (time_axis <= t_on))
        dbs_vec[silent_mask] = 0.0

    if coincident_discontinuities:
        edge_indices = np.where(np.diff(dbs_vec[:-1]))[0]
        time_axis[edge_indices + 1] = time_axis[edge_indices]
    
    return dbs_vec, time_axis