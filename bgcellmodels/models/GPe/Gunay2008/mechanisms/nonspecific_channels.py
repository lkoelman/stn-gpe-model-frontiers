import solvers.mechanism as mechanism
from solvers.mechanism import (
    State, Current, Parameter,
    state_derivative, state_steadystate, state_timeconst,
    mV, ms, S, cm2, nodim,
    v, exp, log,
)

nrn_mechs = [ # suffixes for Na-channels
    'HCN',
    'HCN2',
]

def plot_IV_curves():
    """
    Plot I-V curves for K channel mechanisms
    """
    from bgcellmodels.common.channel_analysis import ivcurve

    from matplotlib import pyplot as plt
    import neuron
    h = neuron.h

    # Load NEURON libraries, mechanisms
    # import os
    # nrn_dll_path = os.path.dirname(__file__)
    # neuron.load_mechanisms(nrn_dll_path)

    # h.CVode().active(1)

    for mech_name in nrn_mechs:
        print('Generating I-V curve for {} ...'.format(mech_name))
        ik, v = ivcurve(mech_name, 'ih_'+mech_name)

        plt.figure()
        plt.plot(v, ik, label='ih_'+mech_name)

        plt.suptitle(mech_name)
        plt.xlabel('v (mV)')
        plt.ylabel('current (mA/cm^2)')
        plt.legend()
        plt.show(block=False)

if __name__ == '__main__':
    # Plot response of mod file channel
    plot_IV_curves()
