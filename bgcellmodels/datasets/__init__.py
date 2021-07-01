"""
Datasets for use in basal ganglia simulations.
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import submodule as submodule_alias
# from .submodule import myclass


import os.path
import numpy as np

here = os.path.abspath(os.path.dirname(__file__))

################################################################################
# Little and Brown, 2012
################################################################################

fs_brown_beta_envelope = 10.0

def load_brown_beta_envelope(remove_outliers=True):
    """
    Load Brown beta envelope waveform.

    Experimental data
    -----------------

    Experimental LFP data that was recorded from a single patient with bilateral 
    STN DBS at the Nuffield Department of Clinical Neurosciences in the 
    University of Oxford (Little and Brown, 2012). 

    The time-ranging magnitude of beta-band activity in the experimental 
    LFP signal with DBS off was calculated by band-pass filtering the LFP data 
    between 10 and 35 Hz, full wave rectifying, and then lowpass filtering the 
    rectified data with a cutoff frequency of 5 Hz. This provided an estimate of
    temporal modulation of the beta frequency band activity of a typical LFP 
    signal in the beta frequency range with DBS turned off. The LFP modulation 
    was centered at zero by subtracting the mean and downsampled to 10 Hz.

    Returns
    -------

    beta_envelope : np.array
        Beta envelope centered at zero and sampled at 10 Hz
    """
    beta_filename = os.path.join(here, 'Brown_Beta_Modulation_Values.txt')
    beta = np.loadtxt(beta_filename, dtype=float)
    thresh_outliers = 0.1
    beta_upper = max(beta[beta<thresh_outliers])
    if remove_outliers:
        beta[beta>=thresh_outliers] = beta_upper
    return beta


def make_brown_modulated_sinusoid(f_sin=20.0, Ts=5.0, tstart=0.0, tstop=1000.0):
    """
    Make sinusoidal waveform modulated by the Brown beta envelope.

    @param  fsin : float
            Base frequency of sinusoid

    @param  Ts : float
            Sampling period (ms) of resulting signal.

    @param  tstart : float
            Starting time within Beta envelope signal

    @param  tstop : float
            Stopping time within Beta envelope signal
    """
    import scipy.signal

    dur = tstop - tstart
    time_ms = np.arange(0.0, dur+Ts, Ts)
    base_sin = np.sin(2.0*np.pi*f_sin*1e-3*time_ms)
    envelope_full = load_brown_beta_envelope(remove_outliers=True)
    Tenv = 1.0/fs_brown_beta_envelope * 1e3
    envelope_seg = envelope_full[int(tstart//Tenv):int(tstop//Tenv)]
    envelope_upsample = scipy.signal.resample(envelope_seg, base_sin.size)
    return base_sin * envelope_upsample