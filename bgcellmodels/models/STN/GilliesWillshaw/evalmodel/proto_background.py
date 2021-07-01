"""
Functions to set up STN experimental protocol with low and high
background activity


Physiological & Anatomical parameters: see notes/BG_observations_experiments.md

"""

# Standard library
import re
from collections import namedtuple

# NEURON
import neuron
h = neuron.h

# Physiological parameters
import bgcellmodels.cellpopdata as cpd
from bgcellmodels.cellpopdata import (
        PhysioState,
        Populations as Pop,
        NTReceptors as NTR,
        ParameterSource as Cit,
        MSRCorrection as MSRC
)

# Stimulation protocols
from bgcellmodels.common import stimprotocols
from bgcellmodels.common.stimprotocols import StimProtocol, EvaluationStep, register_step

from bgcellmodels.common import logutils
logger = logutils.getBasicLogger(name='stn_protos')
# from IPython.config import Application
# logger = Application.instance().log 

# PROBLEM:
#   - we have estimates for the total number of synapses
#   - we have measurements of the post-synaptic response to axonal stimulation
#   - However, this response is most likely due to multiple synapses (MULTI SYNAPSE RULE: axons make multi-synaptic contacts)
#   - so when we calibrate a synapse mechanism to match this response, our synapse mechanism emulates the effect of multiple real synapses

# Global parameters
n_syn_stn_tot = 300     # [SETPARAM] 300 total synapses on STN
frac_ctx_syn = 3.0/4.0  # [SETPARAM] fraction of CTX/GPE of STN afferent synapses

gsyn_single = 0.8e-3    # [SETPARAM] average unitary conductance [uS])

MSR_NUM_SYN = {
    Pop.CTX: 5.0,       # [SETPARAM] average number of contacts per multi-synaptic contact (CTX)
    Pop.GPE: 5.0,       # [SETPARAM] average number of contacts per multi-synaptic contact (GPE)
}

FRAC_SYN = {
    Pop.CTX: frac_ctx_syn,
    Pop.GPE: 1.0 - frac_ctx_syn,
}

MSR_METHOD = MSRC.SCALE_NUMSYN_MSR # How to take into account multi-synapse rule
# NOTE: if method is SCALE_NUMSYN_MSR => n_syn = n_syn_stn_tot * FRAC_SYN / MSR_NUM_SYN

USE_BURSTSTIM = True

################################################################################
# Interface functions
################################################################################

@register_step(EvaluationStep.INIT_SIMULATION, StimProtocol.SYN_BACKGROUND_HIGH)
def init_sim_wrapper(self, protocol):
    """
    Initialize simulator to simulate background protocol
    """

    # Only adjust duration
    self._init_sim(dur=10000)

    # Reset RNGs
    for pre_pop in (Pop.CTX, Pop.GPE):
        input_dict = self.get_inputs(pre_pop, cpd.StnModel.Gillies2005)

        # Reset each instance of Hoc.Random
        for RNG_data in input_dict['RNG_data']:

            # Get RNG and sequence number (highindex)
            random = RNG_data['HocRandomObj']
            start_seq = RNG_data['seq']

            # Reset counter
            old_seq = random.seq()
            random.seq(start_seq)
            logger.anal("Changing RNG seq from {} to {}".format(old_seq, start_seq))


def init_sim_impl(**kwargs):
    """
    Initialize simulator (implementation adhering to interface in proto common).
    """
    stim_data   = kwargs['stim_data']
    h           = kwargs['nrnsim']

    # Initialize phsyiology
    h.celsius = 35
    h.v_init = -60
    h.set_aCSF(4)
    h.init()

    # Reset each instance of Hoc.Random
    for RNG_data in stim_data['RNG_data']:

        # Get RNG and sequence number (highindex)
        random = RNG_data['HocRandomObj']
        start_seq = RNG_data['seq']

        # Reset counter
        old_seq = random.seq()
        random.seq(start_seq)
        logger.anal("Changing RNG seq from {} to {}".format(old_seq, start_seq))


@register_step(EvaluationStep.MAKE_INPUTS, StimProtocol.SYN_BACKGROUND_HIGH)
def make_inputs_wrapper(self, connector=None):
    """
    Make a realistic number of CTX and GPe synapses that fire
    a background firing pattern onto STN.
    """
    model = self.target_model

    # Prepare inputs
    ICell = namedtuple("ICell", ['dendritic'])
    icell = ICell(dendritic=[sec for sec in self.model_data[model]['dend_refs']])

    common_args = {
        'nrnsim':       neuron.h,
        'connector':    connector,
        'evaluator':    self,
        'gid':          self.model_data[model]['gid'],
        'icell':        icell,
        'rng_info':     self.rng_info,
    }

    # Make CTX inputs
    ctx_stim_data = {}
    make_inputs_ctx_impl(
            stim_data = ctx_stim_data,
            **common_args)

    self.add_inputs(Pop.CTX.name.lower(), model, **ctx_stim_data)

    # Make GPE inputs
    gpe_stim_data = {}
    make_inputs_gpe_impl(
            stim_data = gpe_stim_data,
            **common_args)
    
    self.add_inputs(Pop.GPE.name.lower(), model, **gpe_stim_data)


def make_inputs_ctx_impl(**kwargs):
    """
    Make inputs (implementation adhering to interface in proto common).

    Interface keyword arguments:
    ----------------------------

    @param  nrnsim
            NEURON top-level Hoc interpreter


    @param  icell : object
            
            instantianed cell object with attributes:
            seclist_names:  names of available SectionList attributes
            seclist_x:      one SectionList attribute for each
                            name in seclist_names


    @param  stim_data : dict
            
            a dictionary where synapses, NetCons, stimulators etc. can be stored


    @param  rng_info
            
            dictionary with info about random streams, with at least
            the following entries: 'stream_indices': dict(int -> int)



    Additional required keyword arguments:
    --------------------------------------

    @param  gid :int
            global cell identifier

    @param  base_seed : int
            base seed for the RNGs


    Additonal optional keyword arguments:
    -------------------------------------

    @param  evaluator : StnModelEvaluator
            The evaluator object

    @param  connector : CellConnector
            Cell connector for looking up connection parameters
    """
    # Optional arguments
    evaluator = kwargs.get('evaluator', None)
    connector = kwargs.get('connector', None)

    if evaluator is None:
        physio_state    = PhysioState.from_descr(kwargs['physio_state'])
        rng             = kwargs['rng']
        base_seed       = kwargs['base_seed']
    else:
        physio_state = evaluator.physio_state
        rng = evaluator.rng
        base_seed = evaluator.base_seed

    ###########################################################################
    # CTX inputs

    # Filter to select distal, smaller-diam dendritic segments
    is_ctx_target = lambda seg: seg.diam <= 1.0 # see Gillies cell diagram

    # Parameters for making connection
    syn_mech_NTRs = ('GLUsyn', [NTR.AMPA, NTR.NMDA])
    refs_con = [Cit.Chu2015]
    refs_fire = [Cit.Custom, Cit.Bergman2015RetiCh3]

    # Get connection & firing parameters
    POP_PRE = Pop.CTX
    nsyn_per_ax = MSR_NUM_SYN[POP_PRE] # if MSR_METHOD==MSRC.SCALE_GSYN_MSR else 0
    con_par = connector.getPhysioConParams(POP_PRE, Pop.STN, refs_con, 
                            adjust_gsyn_msr=nsyn_per_ax)

    fire_par = connector.getFireParams(POP_PRE, physio_state, refs_fire, 
                            custom_params={'rate_mean': 20.0})

    passed_args = ['stim_data', 'icell', 'gid', 'rng_info']
    passed_kwargs = {k:v for k,v in kwargs.items() if (k in passed_args)}

    # Let user override number of placed synapses
    if 'num_syn_ctx' in kwargs:
        passed_kwargs['num_syn_pp'] = kwargs['num_syn_ctx']

    # Make CTX GLUtamergic inputs
    make_background_inputs(
                    pre_pop         =POP_PRE,
                    syn_seg_filter  =is_ctx_target, 
                    syn_mechs_NTRs  =syn_mech_NTRs,
                    stim_params     =fire_par,
                    con_params      =con_par,
                    connector       =connector,
                    rng             =rng,
                    base_seed       =base_seed,
                    **passed_kwargs)


def make_inputs_gpe_impl(**kwargs):
    """
    Make inputs (implementation adhering to interface in proto common).

    Expected keyword arguments: see make_inputs_ctx_impl

    @see    make_inputs_ctx_impl
    """

    # Optional arguments
    evaluator = kwargs.get('evaluator', None)
    connector = kwargs.get('connector', None)

    if evaluator is None:
        physio_state    = PhysioState.from_descr(kwargs['physio_state'])
        rng             = kwargs['rng']
        base_seed       = kwargs['base_seed']
    else:
        physio_state = evaluator.physio_state
        rng = evaluator.rng
        base_seed = evaluator.base_seed

    ###########################################################################
    # GPe inputs

    # Filter to select proximal, larger-diam dendritic segments
    is_gpe_target = lambda seg: seg.diam > 1.0 # see Gillies cell diagram

    # Parameters for making connection
    syn_mech_NTRs = ('GABAsyn', [NTR.GABAA, NTR.GABAB])
    refs_con = [Cit.Chu2015, Cit.Fan2012, Cit.Atherton2013]
    refs_fire = [Cit.Bergman2015RetiCh3]

    # Get connection & firing parameters
    # SETPARAM: method for scaling gbar of synapse point process
    nsyn_per_ax = MSR_NUM_SYN[Pop.GPE] # if MSR_METHOD==MSRC.SCALE_GSYN_MSR else 0
    con_par = connector.getPhysioConParams(Pop.GPE, Pop.STN, refs_con,
                        adjust_gsyn_msr=nsyn_per_ax)

    fire_par = connector.getFireParams(Pop.GPE, physio_state, refs_fire)

    # Make GPe GABAergic inputs
    passed_args = ['stim_data', 'icell','gid', 'rng_info']
    passed_kwargs = {k:v for k,v in kwargs.items() if (k in passed_args)}

    # Let user override number of placed synapses
    if 'num_syn_gpe' in kwargs:
        passed_kwargs['num_syn_pp'] = kwargs['num_syn_gpe']

    make_background_inputs(
                    pre_pop         =Pop.GPE,
                    syn_seg_filter  =is_gpe_target, 
                    syn_mechs_NTRs  =syn_mech_NTRs,
                    stim_params     =fire_par,
                    con_params      =con_par,
                    connector       =connector,
                    rng             =rng,
                    base_seed       =base_seed,
                    **passed_kwargs)



@register_step(EvaluationStep.RECORD_TRACES, StimProtocol.SYN_BACKGROUND_HIGH)
def rec_traces(self, protocol, traceSpecs):
    """
    Record all traces for this protocol.
    """
    model = self.target_model
    # stim_data_gpe = self.model_data[model]['inputs']['gpe']
    # stim_data_ctx = self.model_data[model]['inputs']['ctx']
    rec_kwargs = {
        'stim_data': self.merged_inputs(['gpe','ctx'], model),
        'trace_specs': traceSpecs,
        'rec_hoc_objects': self.model_data[model]['rec_segs'][protocol],
        'rec_pre_pop_spikes': [Pop.GPE.name, Pop.CTX.name],
    }

    # record synaptic traces
    stimprotocols.rec_GABA_traces(**rec_kwargs)
    stimprotocols.rec_GLU_traces(**rec_kwargs)

    # record membrane voltages
    stimprotocols.rec_Vm(**rec_kwargs)

    # Record input spikes
    stimprotocols.rec_spikes(**rec_kwargs)


@register_step(EvaluationStep.PLOT_TRACES, StimProtocol.SYN_BACKGROUND_HIGH)
def plot_traces(self, model, protocol):
    """
    Plot all traces for this protocol
    """

    # Plot rastergrams
    gpe_filter = lambda trace: re.search('AP_'+Pop.GPE.name, trace)
    self._plot_all_spikes(model, protocol, trace_filter=gpe_filter, color='r')

    ctx_filter = lambda trace: re.search('AP_'+Pop.CTX.name, trace)
    self._plot_all_spikes(model, protocol, trace_filter=ctx_filter, color='g')

    # Plot Vm in select number of segments
    self._plot_all_Vm(model, protocol, fig_per='cell')


################################################################################
# Building block functions
################################################################################


def make_background_inputs(**kwargs):
    """
    Make synapses with background spiking from given population.

    @param refs_fire        References for firing parameters

    @param refs_con         References for connectivity parameters
    """

    POP_PRE         = kwargs['pre_pop']
    cc              = kwargs['connector']
    con_par         = kwargs['con_params']
    fire_par        = kwargs['stim_params']
    syn_mech_NTRs   = kwargs['syn_mechs_NTRs']
    rng             = kwargs['rng']
    icell           = kwargs['icell']
    stim_data       = kwargs['stim_data']
    gid             = kwargs['gid']
    base_seed       = kwargs['base_seed']
    rng_info        = kwargs['rng_info']
    is_target_seg   = kwargs['syn_seg_filter']

    rng_state = rng.get_state()
    rng_pos = rng_state[2]
    logger.debug('Using NumPy Random object with position = {}'.format(rng_pos))

    # RNG info
    highest_indices = rng_info['stream_indices']

    # Get max synaptic conductance
    syn_mech, syn_NTRs = syn_mech_NTRs
    gsyn_multi = max((con_par[NTR].get('gbar', 0) for NTR in syn_NTRs))

    # Calculate number of synapses
    n_syn_single = int(FRAC_SYN[POP_PRE] * n_syn_stn_tot) # number of unitary synapses (actual synaptic contacts) for this population

    if 'num_syn_pp' in kwargs:
        # Let user override number of actual synapses
        n_syn_multi = kwargs['num_syn_pp']
        logger.debug('Override number of synapses for {} : n={}'.format(POP_PRE, n_syn_multi))
    
    elif MSR_METHOD == MSRC.SCALE_NUMSYN_GSYN:
        gsyn_tot = n_syn_single * gsyn_single # total parallel condutance desired [uS]
        n_syn_multi = int(gsyn_tot / gsyn_multi)
    
    elif MSR_METHOD == MSRC.SCALE_NUMSYN_MSR:
        # One synapse represents multiple contacts: divide number of observed synapses by number of contacts per synapse
        n_syn_multi = int(n_syn_single / MSR_NUM_SYN[POP_PRE])

    elif MSR_METHOD == MSRC.SCALE_GSYN_MSR:
        # gbar is adjusted, so don't adjust number of synapses
        n_syn_multi = int(n_syn_single)
    
    n_syn = n_syn_multi
    logger.debug("Number of {}->STN MSR synapses = {}".format(POP_PRE.name, n_syn))

    # Get target segments: distribute synapses over dendritic trees
    dendritic_secs = icell.dendritic
    target_segs = stimprotocols.pick_random_segments(dendritic_secs, n_syn, is_target_seg, rng=rng)

    # Data for configuring inputs
    tstart = 300
    stim_rate = fire_par['rate_mean']
    pause_rate_hz = fire_par.get('pause_rate_mean', 0)
    pause_dur = fire_par.get('pause_dur_mean', 0)
    discharge_dur = fire_par.get('discharge_dur_mean', 0)
    # TODO: set intra-burst rate (higher than mean rate), set burst dur, calculate 'number' from these two, then set controlling stim rate to pause rate

    # Make synapses
    for i_seg, target_seg in enumerate(target_segs):

        # MCellRan4: each stream should be statistically independent as long as the highindex values differ by more than the eventual length of the stream. See http://www.neuron.yale.edu/neuron/static/py_doc/programming/math/random.html?highlight=MCellRan4

        dur_max_ms = 10000.0
        stim_interval = stim_rate**-1*1e3

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

        if pause_rate_hz > 0 and USE_BURSTSTIM: # make bursting spike generator

            stimsource = h.BurstStim()
            stimsource.fast_invl = stim_interval
            stimsource.slow_invl = pause_dur*1e3
            stimsource.burst_len = discharge_dur*stim_rate # ((1.0/pause_rate_hz) - pause_dur) * stim_rate
            stimsource.start = tstart
            stimsource.noise = 1.0
            stimsource.noiseFromRandom(stimrand) # Set it to use this random number generator

        else: # make poisson spike generator
            
            # make NetStim spike generator
            stimsource = h.NetStimExt()
            stimsource.interval = stim_interval # Interval between spikes
            stimsource.start = tstart
            stimsource.number = 1e9 # inexhaustible for our simulation
            stimsource.noise = 1.0
            stimsource.noiseFromRandom(stimrand) # Set it to use this random number generator

        # Make synapse and NetCon
        syn, nc, wvecs = cc.make_synapse((POP_PRE, Pop.STN), (stimsource, target_seg), 
                            syn_mech, syn_NTRs, con_par_data=con_par)

        # Save inputs
        stim_data.setdefault('NetStims', []).append(stimsource)
        stim_data.setdefault('RNG_data', []).append({
            'HocRandomObj': stimrand, 
            'seq': stimrand.seq(),
        })
        stim_data.setdefault('syn_NetCons', []).append(nc)
        stim_data.setdefault('synapses', []).append(syn)
        stim_data.setdefault('stimweightvec', []).append(wvecs)
        

        if pause_rate_hz > 0 and not USE_BURSTSTIM:

            # Sanity check timing params
            assert (discharge_dur < 1/pause_rate_hz), "Discharge duration must be smaller than inter-pause interval"

            # Make spike generator exhaustible
            stimsource.number = discharge_dur*stim_rate # expected number of spikes in discharge duration
            stimsource.ispike = int(discharge_dur*stim_rate*rng.rand())

            pause_interval = pause_dur*1e3
            num_indep_repicks = dur_max_ms / pause_interval + 1000
            low_index = gid+100+base_seed
            highest_index = highest_indices.get(low_index, 0)
            high_index = int(highest_index + num_indep_repicks)
            highest_indices[low_index] = high_index # update highest index

            # Make RNG for spike control
            ctlrand = h.Random()
            ctlrand.MCellRan4(high_index, low_index) # high_index can also be set using .seq()
            ctlrand.negexp(1) # num arrivals is poisson distributed, ISIs are negexp-distributed

            # control spike generator spiking pattern
            stimctl = h.NetStim()
            # stimctl.interval = pause_rate_hz**-1*1e3
            stimctl.interval = pause_interval # replenish only works when spikes exhaused
            stimctl.number = 1e9
            stimctl.noise = 1.0
            stimctl.noiseFromRandom(ctlrand)
            stimctl.start = tstart

            # Connect to spike generator
            # off_nc = h.NetCon(stimctl, stimsource)
            # off_nc.weight[0] = -1 # turn off spikegen
            # off_nc.delay = 0
            
            ctl_nc = h.NetCon(stimctl, stimsource)
            ctl_nc.weight[0] = 1 # turn on spikegen (resets available spikes)
            # ctl_nc.delay = pause_dur*1e3

            # extend_dictitem(new_inputs, 'com_NetCons', off_nc)
            stim_data.setdefault('com_NetCons', []).append(ctl_nc)
            stim_data.setdefault('NetStims', []).append(stimctl)
            stim_data.setdefault('RNG_data', []).append({
                'HocRandomObj': ctlrand, 
                'seq': ctlrand.seq(),
            })



            