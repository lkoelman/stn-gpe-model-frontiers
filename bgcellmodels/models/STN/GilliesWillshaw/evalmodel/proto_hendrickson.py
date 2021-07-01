"""
Synaptic stimulation protocols like in Hendrickson, Edgerton, Jaeger (2011)
"""

# Python stdlib
from collections import namedtuple

# NEURON
import neuron
h = neuron.h

# Physiological parameters
from bgcellmodels.cellpopdata import (
        PhysioState,
        Populations as Pop,
        NTReceptors as NTR,
        ParameterSource as Cit,
        CellConnector
)
from proto_common import pick_random_segments

from bgcellmodels.common import logutils
logger = logutils.getBasicLogger(name='stn_protos')


################################################################################
# Interface functions
################################################################################

def init_sim_impl(**kwargs):
    """
    Initialize simulator (implementation adhering to interface in proto common).
    """
    stim_data   = kwargs['stim_data']
    h           = kwargs['nrnsim']
    tstop       = kwargs.get('tstop', 2000.0)
    dt          = kwargs.get('dt', 0.025)

    # Initialize phsyiology
    h.celsius = 35
    h.v_init = -60
    h.set_aCSF(4)
    h.init()

    h.dt = dt
    h.tstop = tstop

    # Reset each instance of Hoc.Random
    for RNG_data in stim_data['RNG_data']:

        # Get RNG and sequence number (highindex)
        random = RNG_data['HocRandomObj']
        start_seq = RNG_data['seq']

        # Reset counter
        old_seq = random.seq()
        random.seq(start_seq)
        logger.anal("Changing RNG seq from {} to {}".format(old_seq, start_seq))


def make_all_synapses(**kwargs):
    """
    Make inputs (implementation adhering to interface in proto common).

    PARAMETER Arguments
    -------------------

    @param  num_syn_<ampa/nmda/gabaa> : int
            Number of AMPA/NMDA/GABAA synapses

    @param  gbar_<ampa/nmda/gabaa> : int
            Maximum synaptic conductance of AMPA/NMDA/GABAA synapses
    """
    # Optional arguments
    evaluator = kwargs.get('evaluator', None)
    connector = kwargs.get('connector', None)
    
    if evaluator is None:
        # Assume user provided parameters
        physio_state    = PhysioState.from_descr(kwargs['physio_state'])
        rng             = kwargs['rng']
        base_seed       = kwargs['base_seed']
    else:
        # Get params from Evaluator
        physio_state = evaluator.physio_state
        rng = evaluator.rng
        base_seed = evaluator.base_seed
    
    if connector is None:
        connector = CellConnector(physio_state, rng)

    # Required arguments
    icell           = kwargs['icell']
    stim_data       = kwargs['stim_data']
    gid             = kwargs['gid']
    rng_info        = kwargs['rng_info']
    
    # Random Number Generation
    rng_state = rng.get_state()
    rng_pos = rng_state[2]
    logger.debug('Using NumPy Random object with position = {}'.format(rng_pos))
    highest_indices = rng_info['stream_indices']

    # NOTE: in exp2syn.mod : NetCon.weight is interpreted in units [uS]
    input_types = [
        {
            'neurotransmitter': NTR.AMPA,
            'mech':             'Exp2Syn',
            'stim_rate':        5.0,
            'tstart':           300.0,
            'num_syn':          kwargs['num_syn_ampa'],
            'con_par': {
                'gbar':         kwargs.get('gbar_ampa', 0.25e-3),
                'tau_rise_g':   1.0,
                'tau_decay_g':  3.0,
                'Erev':         0.0,
            },
            'target_func':      lambda seg: True,
            'pre_pop':          Pop.STR, # only used for labeling
        },
        {
            'neurotransmitter': NTR.NMDA,
            'mech':             'Exp2Syn',
            'stim_rate':        5.0,
            'tstart':           300.0,
            'num_syn':          kwargs['num_syn_nmda'],
            'con_par': {
                'gbar':         kwargs.get('gbar_nmda', 0.25e-3),
                'tau_rise_g':   10.0,
                'tau_decay_g':  30.0,
                'Erev':         0.0,
            },
            'target_func':      lambda seg: True,
            'pre_pop':          Pop.CTX, # only used for labeling
        },
        {
            'neurotransmitter': NTR.GABAA,
            'mech':             'Exp2Syn',
            'stim_rate':        5.0,
            'tstart':           300.0,
            'num_syn':          kwargs['num_syn_gabaa'],
            'con_par': {
                'gbar':         kwargs.get('gbar_gabaa', 0.25e-3),
                'tau_rise_g':   1.0,
                'tau_decay_g':  12.0,
                'Erev':         -80.0,
            },
            'target_func':      lambda seg: True,
            'pre_pop':          Pop.STN, # only used for labeling
        },
    ]


    for i_input, input_data in enumerate(input_types):

        # Get target segments: distribute synapses over dendritic trees
        target_segs = pick_random_segments(
                        icell.dendritic,
                        input_data['num_syn'],
                        input_data['target_func'],
                        rng=rng)

        post_pop = Pop.GPE

        for i_seg, target_seg in enumerate(target_segs):

            # MCellRan4: each stream should be statistically independent as long as the highindex values differ by more than the eventual length of the stream. See http://www.neuron.yale.edu/neuron/static/py_doc/programming/math/random.html?highlight=MCellRan4

            dur_max_ms = 10000.0
            stim_interval = input_data['stim_rate']**-1*1e3

            # RNG settings
            num_indep_repicks = dur_max_ms / stim_interval + 1000
            low_index = gid+250+base_seed
            highest_index = highest_indices.get(low_index, 0)
            high_index = int(highest_index + num_indep_repicks)
            highest_indices[low_index] = high_index # update highest index

            # Make RNG for spikes
            stimrand = h.Random() # see CNS2014 Dura-Bernal example or EPFL cell synapses.hoc file
            stimrand.MCellRan4(high_index, low_index) # high_index can also be set using .seq()
            stimrand.negexp(1) # if num arrivals is poisson distributed, ISIs are negexp-distributed

            # make NetStim spike generator
            stimsource = h.NetStimExt()
            stimsource.interval = stim_interval # Interval between spikes
            stimsource.start = input_data['tstart']
            stimsource.number = 1e9 # inexhaustible for our simulation
            stimsource.noise = 1.0
            stimsource.noiseFromRandom(stimrand) # Set it to use this random number generator

            # Make synapse and NetCon
            syn, nc, wvecs = connector.make_synapse(
                                (input_data['pre_pop'], post_pop),
                                (stimsource, target_seg), 
                                input_data['mech'],
                                (input_data['neurotransmitter'],),
                                con_par_data=input_data['con_par'])

            # Save inputs
            stim_data.setdefault('NetStims', []).append(stimsource)
            stim_data.setdefault('RNG_data', []).append({
                'HocRandomObj': stimrand, 
                'seq': stimrand.seq(),
            })
            stim_data.setdefault('syn_NetCons', []).append(nc)
            stim_data.setdefault('synapses', []).append(syn)
            stim_data.setdefault('stimweightvec', []).append(wvecs)


if __name__ == '__main__':
    print("""
    Running example script for proto_hendrickson...
    """)

    import re, functools
    import proto_common
    from stn_model_evaluation import StnModel, StnModelEvaluator

    # Make cell model and evaluator
    full_model = StnModel.Gillies2005
    evaluator = StnModelEvaluator(full_model, PhysioState.NORMAL)

    # Register protocol functions
    init_funcs = [init_sim_impl]
    setup_funcs = [make_all_synapses]
    rec_funcs = [proto_common.rec_spikes]

    plot_funcs = [
        functools.partial(proto_common.plot_all_sikes,
                          trace_filter=lambda trace: re.search('AP_'+Pop.CTX.name, trace)),
        functools.partial(proto_common.plot_all_sikes,
                          trace_filter=lambda trace: re.search('AP_'+Pop.STR.name, trace)),
    ]

    setup_kwargs = {
        'num_syn_ampa': 100,
        'num_syn_nmda': 100,
        'num_syn_gabaa': 0,
        'rec_pre_pop_spikes': ['str', 'ctx', 'stn'],
    }

    evaluator.register_keyword_funcs(
                proto_init_funcs = init_funcs,
                proto_setup_funcs_pre = setup_funcs,
                proto_rec_funcs = rec_funcs,
                proto_plot_funcs = plot_funcs,
                proto_setup_kwargs_const = setup_kwargs)

    evaluator.setup_keyword_protocol(pre_model=full_model)