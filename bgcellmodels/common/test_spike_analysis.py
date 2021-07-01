"""
Test suite for signal.py and spike_analysis.jl
"""

import numpy as np
from bgcellmodels.common import spikelib, signal

spike_gen = spikelib.make_oscillatory_bursts(50.0, 30.0, 180.0, 20.0, 5e3)
spike_times = np.fromiter(spike_gen, float)
isi_values = np.diff(spike_times)


def test_burst_metrics_surprise_python():
    """
    Test Python implementation of burst_metrics_surprise()
    """
    metrics = signal.burst_metrics_surprise(np.diff(spike_times))
    print(metrics)


def test_burst_metrics_surprise_julia():
    """
    Test Julia implementation of burst_metrics_surprise()
    """
    print("Loading Julia interpreter ...")
    from julia import Main as jl
    jl.include("spike_analysis.jl")

    print("Executing algorithm ...")
    metrics = jl.burst_metrics_surprise(isi_values)
    print(metrics)


if __name__ == '__main__':
    test_burst_metrics_surprise_python()
    test_burst_metrics_surprise_julia()