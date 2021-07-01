import solvers.mechanism as mechanism
from solvers.mechanism import (
    State, Current, Parameter,
    state_derivative, state_steadystate, state_timeconst,
    mV, ms, S, cm2, nodim,
    v, exp, log,
)

nrn_mechs = [ # suffixes for K-channels
    'Kv2',
    'Kv3',
    'Kv4f',
    'Kv4s',
    'KCNQ',
    'SK',
]

class Kv3_channel(mechanism.MechanismBase):
    """
    fast activated potassium Kv3 (Kv3.1/3.4) channel for GPe neuron
    """

    # PARAMETER block
    # gmax = Parameter('gmax', 0.001, 'S/cm2')
    # ek = Parameter('ek', -90.0, 'mV')
    # theta_m0 = Parameter('theta_m0', -13.0, 'mV')
    # k_m = Parameter('k_m', 7.8, 'mV')

    @classmethod
    def PARAMETERS(mech):
        """
        Define mechanism parameters
        """

        gmax = 0.001 * S/cm2
        ek = -90.0 * mV

        theta_m0 = -13.0 * mV
        k_m = 7.8 * mV
        tau_m0 = 0.1 * (ms)
        tau_m1 = 14.0 * (ms)
        sigma_m0 = -12.0 * (mV)
        sigma_m1 = -13.0 * (mV)

        ## h-gate
        h0 = 0.6 * nodim
        theta_h = -20.0 * (mV)
        k_h = 10.0 * (mV)
        tau_h0 = 7.0 * (ms)
        tau_h1 = 33.0 * (ms)
        phi_h = 0.0 * (mV)
        sigma_h0 = 10.0 * (mV)

        super(mech, mech).PARAMETERS() # NOTE: cannot use explicit MechanismBase because we want base method with arg equal to subclass

    @classmethod
    def DYNAMICS(mech):
        """
        Define mechanism dynamics.
        """
        # self.inject_parameters() # equivalent to locals().update(self._IMECH_PARAMS)

        # STATE block
        # TODO: make sympy symbol + register power
        m = State('m', power=4)
        h = State('h', power=1)

        # RHS_EXPRESSIONS
        ## m gate
        theta_m = mech.theta_m0 + (mech.k_m * log((1.0 / pow(0.5, 1.0/4.0)) - 1.0))

        minf = state_steadystate(
                    1.0 / (1.0 + exp((theta_m - v)/mech.k_m)),
                    state='m')
        
        tau_m = state_timeconst(
                    mech.tau_m0 + (mech.tau_m1 - mech.tau_m0) / (
                        exp((theta_m - v)/mech.sigma_m0) + exp(-(theta_m - v)/mech.sigma_m1)),
                    state='m')

        ## h gate
        hinf = state_steadystate(
                    mech.h0 + (1.0 - mech.h0) / (1.0 + exp((v - mech.theta_h)/mech.k_h)),
                    state='h')
        
        tau_h = state_timeconst(
                    mech.tau_h0 + (mech.tau_h1 - mech.tau_h0)/(1 + exp((v - mech.phi_h)/mech.sigma_h0)),
                    state='h')

        # DERIVATIVE block
        dm = state_derivative((minf - m) / tau_m, state='m')
        dh = state_derivative((hinf - h) / tau_h, state='h')

        # BREAKPOINT block
        ik = Current(mech.gmax * m**m.power * h**h.power * (v-mech.ek), ion='k')

        super(mech, mech).DYNAMICS()

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
        ik, v = ivcurve(mech_name, 'ik')

        plt.figure()
        plt.plot(v, ik, label='ik_' + mech_name)

        plt.suptitle(mech_name)
        plt.xlabel('v (mV)')
        plt.ylabel('current (mA/cm^2)')
        plt.legend()
        plt.show(block=False)

if __name__ == '__main__':
    # Create channel
    # chan = Kv3_channel()
    # chan.plot_steadystate_gating()

    # Plot response of mod file channel
    plot_IV_curves()