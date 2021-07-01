"""
Test cases for synapses and connections

@author Lucas Koelman
@date   17/10/2018
"""

# PyNN library
import pyNN
import pyNN.neuron as sim
from pyNN.utility import init_logging

# Custom PyNN extensions
import bgcellmodels.extensions.pynn.synapses as custom_synapses

# Custom cell models
import bgcellmodels.models.STN.GilliesWillshaw.gillies_pynn_model as gillies
from bgcellmodels.mechanisms import noise # loads MOD files



def test_MultiMechanismConnection(export_locals=True):
    """
    Test case for MultiMechanismConnection and NativeMultiSynapse.
    """

    init_logging(logfile=None, debug=True)
    sim.setup(timestep=0.025, min_delay=0.1, max_delay=10.0, use_cvode=False)

    # Global vars used by our custom cell models
    sim.state.duration = 500.0 
    sim.state.rec_dt = 0.05
    sim.state.mcellran4_rng_indices = {} 
    sim.state.shared_rng_seed = shared_seed = 888
    sim.state.rank_rng_seed = rank_seed = sim.state.native_rng_baseseed + sim.state.mpi_rank
    sim.state.shared_rng = sim.NumpyRNG(seed=shared_seed)
    sim.state.rank_rng = sim.NumpyRNG(seed=rank_seed)

    # --------------------------------------------------------------------------
    # CTX Population
    pop_ctx = sim.Population(
                    20,
                    sim.SpikeSourcePoisson(rate=30),
                    label='CTX')

    # --------------------------------------------------------------------------
    # STN Population
    stn_type = gillies.StnCellType(
                        calculate_lfp=False,
                        membrane_noise_std=0.0075)

    initial_values = {
        'v': pyNN.random.RandomDistribution('uniform', (60, 70))
    }

    pop_stn = sim.Population(5, 
                         cellclass=stn_type, 
                         label='STN',
                         initial_values=initial_values)

    # --------------------------------------------------------------------------
    # CTX -> STN Connection
    ctx_stn_conn = sim.FixedNumberPreConnector(10)

    ctx_stn_syn = custom_synapses.NativeMultiSynapse(**{
        'mechanisms_receptors':{
            'AMPAsynTM': 'distal.AMPA',
            'NMDAsynTM': 'proximal.NMDA'
        },
        'weight':       1.0,
        'delay':        1.0, # [ms] delay from literature
        # AMPA receptor
        'AMPAsynTM_U1':           0.1, # baseline release probability
        'AMPAsynTM_tau_rec':      200.0, # [ms] recovery from STD
        'AMPAsynTM_tau_facil':    800.0, # [ms] recovery from facilitation
        'AMPAsynTM_gmax_AMPA':    0.25e-3, # [uS], adjusted for U1
        'AMPAsynTM_tau_r_AMPA':   1.0, # [ms] rise time
        'AMPAsynTM_tau_d_AMPA':   4.0, # [ms] decay time
        # NMDA receptor
        'NMDAsynTM_U1':           0.1, # baseline release probability
        'NMDAsynTM_tau_rec':      200.0, # [ms] recovery from STD
        'NMDAsynTM_tau_facil':    800.0, # [ms] recovery from facilitation
        'NMDAsynTM_gmax_NMDA':    0.251e-3, # [uS], adjusted for U1
        'NMDAsynTM_tau_r_NMDA':   3.7,    # [ms] rise time
        'NMDAsynTM_tau_d_NMDA':   80.0,   # [ms] decay time
    })
    

    ctx_stn_proj = sim.Projection(pop_ctx, pop_stn, 
                                 connector=ctx_stn_conn,
                                 synapse_type=ctx_stn_syn,
                                 receptor_type='distal.AMPA+NMDA')

    # --------------------------------------------------------------------------
    # Simulate
    sim.run(sim.state.duration)
    print("Simulation finished.")

    if export_locals:
        globals().update(locals())

if __name__ == '__main__':
    test_MultiMechanismConnection()