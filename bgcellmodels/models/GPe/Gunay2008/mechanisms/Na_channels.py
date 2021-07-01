import solvers.mechanism as mechanism
from solvers.mechanism import (
    State, Current, Parameter,
    state_derivative, state_steadystate, state_timeconst,
    mV, ms, S, cm2, nodim,
    v, exp, log,
)

nrn_mechs = [ # suffixes for Na-channels
    'NaF',
    'NaP',
]


class NaP_channel(mechanism.MechanismBase):
    """
    Persistent Na channel

    Based on Magistretti & Alonso (1999), JGP 114:491-509
    and Magistretti & Alonso (2002), JGP 120: 855-873.
    """

    @classmethod
    def PARAMETERS(mech):
        """
        Define mechanism parameters
        """

        gmax  = 0.001 * S/cm2
        ena = 50.0 * mV

        Q10 = 3.0 * nodim

        # Activation & Deactivation (m-gate)
        theta_m0 = -50.0 * mV
        k_m = 5.7 * mV
        phi_m = -41.6 * mV
        sigma_m = 14.4 * mV
        alp0 = 2.130 * 1/ms
        bet0 = 2.460 * 1/ms

        # Slow Inactivation (s-gate)
        h0 = 0.154 * nodim
        theta_h = -57.0 * mV
        k_h = -4.0 * mV
        tau_h0 = 10.0 * ms
        tau_h1 = 17.0 * ms
        phi_h = -34.0 * mV
        sigma_h0 = -26.0 * mV
        sigma_h1 = -31.9 * mV

        # Slow Inactivation (s-gate)
        theta_s = -10.0 * mV
        k_s = -4.9 * mV
        Aa_s = -2.88e-6 * 1/ms/mV
        Ba_s = -4.9e-5 * 1/ms
        Ka_s = 4.63 * mV
        Ab_s = 6.94e-6 * 1/ms/mV
        Bb_s = 4.47e-4 * 1/ms
        Kb_s = -2.63 * mV

        super(mech, mech).PARAMETERS() # NOTE: cannot use explicit MechanismBase because we want base method with arg equal to subclass

    @classmethod
    def DYNAMICS(mech):
        """
        Define mechanism dynamics.
        """
        # self.inject_parameters() # equivalent to locals().update(self._IMECH_PARAMS)

        # STATE block
        # TODO: make sympy symbol + register power
        m = State('m', power=3)
        h = State('h', power=1)
        s = State('s', power=1)

        # RHS_EXPRESSIONS

        T_Q10 = mech.Q10 ** ((35-22)/10)

        ## m gate
        # TODO: fix rate, maybe use alpha/beta formulation
        # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=1075
        # https://www.neuron.yale.edu/phpBB/viewtopic.php?f=16&t=110
        theta_m = mech.theta_m0 + (mech.k_m * log((1.0 / pow(0.5, 1.0/3.0)) - 1.0))

        minf = state_steadystate(
                    1.0 / (1.0 + exp((theta_m - v)/mech.k_m)),
                    state='m')

        alpha_m = mech.alp0 * exp((v - mech.phi_m) / mech.sigma_m) * T_Q10
        beta_m = mech.bet0 * exp((mech.phi_m - v) / mech.sigma_m) * T_Q10
        
        tau_m = state_timeconst(
                    1.0 / (alpha_m + beta_m),
                    state='m')

        ## h gate
        hinf = state_steadystate(
                    mech.h0 + (1.0 - mech.h0)/ (1.0 + exp((mech.theta_h - v)/mech.k_h)),
                    state='h')
        
        tau_h = state_timeconst(
                    mech.tau_h0 + (mech.tau_h1 - mech.tau_h0) / 
                        (exp((v - mech.phi_h)/mech.sigma_h0) + exp((mech.phi_h - v)/mech.sigma_h1)),
                    state='h')
        

        ## s gate
        sinf = state_steadystate(
                    1.0 / (1.0 + exp((mech.theta_s - v)/mech.k_s)),
                    state='s')

        alphas = (mech.Aa_s * v + mech.Ba_s)/(1.0 - exp((v + mech.Ba_s/mech.Aa_s)/mech.Ka_s)) * T_Q10
        betas = (mech.Ab_s * v + mech.Bb_s)/(1.0 - exp((v + mech.Bb_s/mech.Ab_s)/mech.Kb_s)) * T_Q10
        
        tau_s = state_timeconst(
                    1.0 / (alphas + betas),
                    state='s')

        # DERIVATIVE block
        dm = state_derivative((minf - m) / tau_m, state='m')
        dh = state_derivative((hinf - h) / tau_h, state='h')
        ds = state_derivative((sinf - s) / tau_s, state='s')

        # BREAKPOINT block
        ina = Current(mech.gmax * m**m.power * h**h.power * s**s.power * (v-mech.ena), ion='na')

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
        ik, v = ivcurve(mech_name, 'ina')

        plt.figure()
        plt.plot(v, ik, label='ina_' + mech_name)

        plt.suptitle(mech_name)
        plt.xlabel('v (mV)')
        plt.ylabel('current (mA/cm^2)')
        plt.legend()
        plt.show(block=False)

if __name__ == '__main__':
    # Create channel
    # chan = NaP_channel()
    # chan.plot_steadystate_gating()

    # Plot response of mod file channel
    plot_IV_curves()
