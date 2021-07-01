import solvers.mechanism as mechanism
from solvers.mechanism import (
    State, Current, Parameter,
    state_derivative, state_steadystate, state_timeconst,
    mV, ms, S, cm2, m2, nodim,
    v, exp, log,
)

import math
PI = math.pi

nrn_mechs = [ # suffixes for Ca-channels
    'CaHVA',
    'CaHVA2',
    'Calcium', # not a channel; calcium buffering mechanism
]


class CaHVA_channel(mechanism.MechanismBase):
    """
    Persistent Na channel

    Based on Magistretti & Alonso (1999), JGP 114:491-509
    and Magistretti & Alonso (2002), JGP 120: 855-873.
    """

    @classmethod
    def PARAMETERS(mech):
        """
        Define mechanism parameters

        NOTES from GENESIS manual
        -------------------------

        See documentation for 'setuptau': 
        http://genesis-sim.org/GENESIS/Hyperdoc/Manual-25.html#ss25.160

        And documentation for 'tabchannel':
        http://genesis-sim.org/GENESIS/Hyperdoc/Manual-26.html#ss26.62
        """

        G_Ca_HVA_dend = 0.15 * S/m2
        G_Ca_HVA_soma = 2.0 * S/m2
        
        gmax  = G_Ca_HVA_dend.to('S/cm2')
        eca = 130.0 * mV

        tau =  0.2 * ms   # estimated from traces Q10=2.5 adjusted
        K = -7.0 * mV # -1*{Kn_CaHVA} (V)
        V0 = -20.0 * mV # {Vhalfn_CaHVA} (V)

        # setuptau Ca_HVA_GP X  {tau} {tau*1e-6} 0 0 1e6 1 0 1 {-1.0*V0} {K} -range {xmin} {xmax}
        # unit scaling: D and F need to be scaled from V to mV

        Atau = tau
        Btau = tau * 1e-6
        Ctau = 0.0 * nodim
        Dtau = 0.0 * mV
        Ftau = 1e9 * mV
        # def tauX(v,A,B,C,D,F):
        #     # A,B,C,D,F = tau, tau*1e-6, 0., 0., 1e9
        #     return (A + B * v) / (C + exp((v + D) / F))

        Ainf = 1.0 * nodim
        Binf = 0.0 * 1/mV
        Cinf = 1.0 * nodim
        Dinf = -V0 # already converted
        Finf = K   # already converted
        # def Xinf(v,A,B,C,D,F):
        #     # A,B,C,D,F = 1., 0., 1., -1.0*V0, K
        #     return (A + B * v) / (C + exp((v + D) / F))


        super(mech, mech).PARAMETERS() # NOTE: cannot use explicit MechanismBase because we want base method with arg equal to subclass

    @classmethod
    def DYNAMICS(mech):
        """
        Define mechanism dynamics.
        """

        # STATE block
        m = State('m', power=1)

        # RHS expressions
        ## m gate
        minf = state_steadystate(
                    (mech.Ainf + mech.Binf * v) / (mech.Cinf + exp((v + mech.Dinf) / mech.Finf)),
                    state='m')

        
        tau_m = state_timeconst(
                    (mech.Atau + mech.Btau * v) / (mech.Ctau + exp((v + mech.Dtau) / mech.Ftau)),
                    state='m')

        # DERIVATIVE block
        dm = state_derivative((minf - m) / tau_m, state='m')

        # BREAKPOINT block
        ina = Current(mech.gmax * m**m.power * (v-mech.eca), ion='ca')

        super(mech, mech).DYNAMICS()


def plot_Ca_buffering(export_locals=False):
    """

    NOTES from GENESIS manual
    -------------------------

    http://genesis-sim.org/GENESIS/Hyperdoc/Manual-26.html#ss26.1

    Single shell model for Ca concentration.
    Solves  dC/dt = B*Ik - C/tau.
    Ca = Ca_base + C.

    In SI units, where concentration is moles/m^3
    (milli-moles/liter) and current is in amperes, theory gives B
    = 5.2e-6/(shell volume).  In practice, B is a parameter to be
    fitted or estimated from experiment, as buffering, non-uniform
    distribution of Ca, etc., will modify this value.  If thick =
    0, the readcell routine calculates B by dividing the "density"
    parameter in the cell parameter file by the volume of the
    compartment.  Otherwise, it scales as a true shell, with the
    volume of a shell having thickness thick.  A negative value of
    the "density" parameter may be used to indicate that it should
    be taken as an absolute value of B, without scaling.  
    
    """

    # Hendrickson (2011) - GPchans.g
    print(plot_Ca_buffering.__doc__)
    
    ## Parameters in simdefaults.g
    B_Ca_GP_conc = 4.0/3.0*5.2e-12
    shell_thick  = 20e-9        #  meters (NOTE: seems to be global param: same for all compartments)
    tau_CaClearance = 0.001     #  time constant for Ca2+ clearance (sec)
    
    ## Settings in GPchans.g (creation of mechanism)
    class Ca_GP_conc:
        """ Ca buffering mechanism """
        
        tau = tau_CaClearance   # sec
        B = 5.2e-6              # Moles per coulomb, later scaled to conc
        Ca_base = 5e-05         # Units in mM, so = 50 nM.
        
        def __init__(self):
            self.tau = Ca_GP_conc.tau
            self.B = Ca_GP_conc.B
            self.Ca_base = Ca_GP_conc.Ca_base
    
    ## Settings in GPcomps.g
    ### Soma sections
    dia = 1. 
    rad = dia/2.
    rad_core = rad - shell_thick # shell_thick is in (m) so rad too
    surf = 4*PI*rad*rad
    vol = 4.0/3.0*PI*rad*rad*rad
    core_vol = 4.0/3.0*PI*rad_core*rad_core*rad_core
    shell_vol = vol - core_vol # m^3
    
    soma_Ca_conc = Ca_GP_conc()
    soma_Ca_conc.B = B_Ca_GP_conc / shell_vol

    ### Axon hillock sections
    L = 1.
    dia = 1.
    rad = dia / 2.
    surf = 2*PI*rad*L
    vol = PI*rad*rad*L
    if dia > (shell_thick*2):
        rad_core = rad - shell_thick
        core_vol = PI*rad_core*rad_core*L
        shell_vol = vol - core_vol
    else:
        shell_vol = vol

    hillock_Ca_conc = Ca_GP_conc()
    hillock_Ca_conc.B = B_Ca_GP_conc / shell_vol

    ### Dendritic sections
    L = 1.
    dia = 1.
    rad = dia / 2.
    surf = 2*PI*rad*L
    vol = PI*rad*rad*L
    if dia > (shell_thick*2):
        rad_core = rad - shell_thick
        core_vol = PI*rad_core*rad_core*L
        shell_vol = vol - core_vol
    else:
        shell_vol = vol

    dend_Ca_conc = Ca_GP_conc()
    dend_Ca_conc.B = B_Ca_GP_conc / shell_vol


    if export_locals:
        globals().update(locals())

    return soma_Ca_conc, hillock_Ca_conc, dend_Ca_conc


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
        
        if mech_name=='Calcium':
            quantity = 'cai'
            descr = 'concentration (mM)'
        else:
            quantity = 'ica'
            descr = 'current (mA/cm^2)'
            
        
        ik, v = ivcurve(mech_name, quantity)

        plt.figure()
        plt.plot(v, ik, label=quantity+'_'+mech_name)

        plt.suptitle(mech_name)
        plt.xlabel('v (mV)')
        plt.ylabel(descr)
        plt.legend()
        plt.show(block=False)

if __name__ == '__main__':
    # Create channel
    chan = CaHVA_channel()
    chan.plot_steadystate_gating()

    # Plot response of mod file channel
    # plot_IV_curves()
