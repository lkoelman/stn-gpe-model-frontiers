"""
Evaluation of STN cell models under different physiological and stimulus conditions.
"""

# Standard library modules
import collections
import re
import pickle
from copy import deepcopy

import logging
logging.basicConfig(format='%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logger = logging.getLogger('stn_protos') # create logger for this module

# Third party modules
import numpy as np
from scipy import signal

# NEURON modules
import neuron
from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# Gillies-Willshaw STN model
import gillies_model as gillies

# Plotting & recording
from bgcellmodels.common import analysis

# Physiological parameters
import bgcellmodels.cellpopdata as cpd
from bgcellmodels.cellpopdata import (
    StnModel,
    PhysioState,
    Populations,
    NTReceptors as NTR,
    ParameterSource as Cit
)
Pop = Populations

# Experimental protocols
import proto_common
proto_common.logger = logger
from proto_common import (
    StimProtocol, SynapticProtocols, EvaluationStep,
    extend_dictitem, pick_random_segments
)

# Importing protocols will run module and register functions
import proto_gillies_article
import proto_simple_syn as proto_simple
import proto_background
import proto_passive_syn

# Adjust verbosity of loggers
logging.getLogger('redops').setLevel(logging.WARNING)
logging.getLogger('folding').setLevel(logging.WARNING)
logging.getLogger('marasco').setLevel(logging.WARNING)


class StnModelEvaluator(object):
    """
    Evaluate STN models

    Inspired by:
    - optimization.py/Simulation
    - optimization.py/StnCellController
    - bgmodel/models/kang/model.py/BGSim

    Improvements:
    - make evaluator subclasses for the different protocols
    """

    # SIGNATURE: make_inputs(evaluator, connector)
    _MAKE_INPUT_FUNCS = {}

    # SIGNATURE: rec_traces(evaluator, protocol, traceSpecs)
    _REC_TRACE_FUNCS = {}

    # SIGNATURE: plot_traces(evaluator, model, protocol)
    _PLOT_TRACE_FUNCS = {}

    # SIGNATURE: init_sim(evaluator, protocol)
    _INIT_SIM_FUNCS = {}

    # Make accessible by step
    _EVALUATION_STEP_FUNCS = {
        EvaluationStep.INIT_SIMULATION :    _INIT_SIM_FUNCS,
        EvaluationStep.MAKE_INPUTS :        _MAKE_INPUT_FUNCS,
        EvaluationStep.RECORD_TRACES :      _REC_TRACE_FUNCS,
        EvaluationStep.PLOT_TRACES :        _PLOT_TRACE_FUNCS,
    }

    # Fill dictionaries with registered functions (see proto_common.register_step)
    for step in list(EvaluationStep):
        for proto in list(StimProtocol):
            
            # If no function set for this protocol and step
            if (proto not in _EVALUATION_STEP_FUNCS[step].keys()):

                # Check if a function was registered for this step and protocol
                step_func = proto_common.EVALUATION_FUNCS[proto].get(step, None)
                if step_func is not None:
                    _EVALUATION_STEP_FUNCS[step][proto] = step_func


    def __init__(self, target_model, physio_state=PhysioState.NORMAL):
        """
        Initialize new evaluator in given physiological state.
        """
        self._physio_state = physio_state

        self.model_data = dict(((model, {}) for model in list(StnModel)))
        self.target_model = target_model

        # Initialize containers for model data
        for model in list(StnModel):
            self.model_data[model]['rec_data'] = {proto: {} for proto in StimProtocol}
            self.model_data[model]['rec_segs'] = {proto: {} for proto in StimProtocol}
            self.model_data[model]['inputs'] = {} # key for each pre-synaptic population
            self.model_data[model]['user_params'] = {} # parameters for model building
            self.model_data[model]['proto_vars'] = {proto: {} for proto in StimProtocol}

        self.sim_dur = 1000.
        self.sim_dt = 0.025
        self.base_seed = 25031989 # used: 25031989
        self.rng = np.random.RandomState(self.base_seed)
        self.rng_info = {
            'base_seed': self.base_seed,
            'stream_indices': {},
        }

        # Functions to set up stimulation protocols, with signature taking keyword arguments
        # (for compatibility with Ephys SelfContainedProtocol)
        self._iproto_data = {} # instantiated protocol
        
        # init_physio(sim, icell)
        self._proto_init_funcs = [] # initialization of simulation, physiology
        self._proto_rec_funcs = []
        self._proto_plot_funcs = []

        # Protocol setup functions applied to full model
        self._proto_setup_funcs_pre = []
        # Protocol setup functions applied to reduced model
        self._proto_setup_funcs_post = []

        # Common keyword arguments for all functions
        self._proto_setup_kwargs_const = {}
        self._proto_setup_kwargs_getters = collections.OrderedDict()
        self._proto_setup_kwargs_setters = []

        # Data used by protocol setup functions
        self.recorded_trace_vectors = {}


    @property
    def physio_state(self):
        """
        Get cell physiological state.
        """
        return self._physio_state


    @physio_state.setter
    def set_physio_state(self, state):
        """
        Set cell physiological state.
        """
        # ! All changes should be enacted in build_cell() and make_inputs()
        # Set the state flag
        self._physio_state = state


    def build_cell(self, model, state=None):
        """
        Build cell model using current physiological state

        @param state        PhysioState enum member
                            (if none given, use self.physio_state)
        """
        if state is None:
            state = self.physio_state

        if self.model_data[model].get('built', False):
            logger.warning("Attempting to build model {} which has already been built.".format(
                            model))

        if model == StnModel.Gillies2005:

            # Build gillies STN cell model
            soma_refs, dend_refs = self._build_gillies(state)

            self.model_data[model]['sec_refs'] = {
                'soma': soma_refs[0],
                'dendrites': list(gillies.get_each_dend_refs(dend_refs))
            }

        elif model in cpd.ReducedModels:

            # Reduce Gillies model
            soma_refs, dend_refs = self._reduce_map_gillies(model)

            dendL = [ref for ref in dend_refs if 'dend0' in ref.sec.name()]
            dendR = [ref for ref in dend_refs if 'dend1' in ref.sec.name()]
            
            self.model_data[model]['sec_refs'] = {
                'soma': soma_refs[0],
                'dendrites': dendL + dendR,
            }

        else:
            raise Exception("Model '{}' not supported".format(
                    model))

        # Set Gillies mechanism
        self.model_data[model]['mechs_gbars_dict'] = gillies.gillies_gdict
        self.model_data[model]['gleak_name'] = gillies.gleak_name
        self.model_data[model]['active_gbar_names'] = gillies.active_gbar_names

        # Save Sections
        self.model_data[model]['soma_refs'] = soma_refs
        self.model_data[model]['dend_refs'] = dend_refs

        # Indicate that given model has been built
        self.model_data[model]['built'] = True
        self.model_data[model]['gid'] = 1 # we only have one cell


    def _build_gillies(self, state):
        """
        Build the Gillies2005 model.
        """

        # Make Gillies STN cell
        soma, dends, stims = gillies.stn_cell_gillies()
        somaref, dendL_refs, dendR_refs = gillies.get_stn_refs()
        
        soma_refs = [somaref]
        dend_refs = dendL_refs + dendR_refs
        all_refs = soma_refs + dend_refs

        # State-dependent changes to model specification
        if state == PhysioState.NORMAL:
            pass # this case corresponds to default model parameters

        elif (state == PhysioState.PARKINSONIAN or state == PhysioState.PARK_DBS):

            # TODO: decide parameter modifications for DA-depleted state from literature

            # 1. Reduce sKCA channel conductance by 90%, from sources:
            #   - Gillies & Willshaw 2005 (see refs)
            for secref in all_refs:
                for seg in secref.sec:
                    seg.gk_sKCa = 0.1 * seg.gk_sKCa

            # 2. Modifications to GPE GABA IPSPs
            #   - DONE: see changes in cellpopdata.py
            #   - Changes:
            #       - Increased strength of GABA IPSPs
            #       - longer decay kinetics, 
            #       - increase in number of functional synapses (1 afferent axon has more activated synaptic contacts)
            #   - References:
            #       - Fan (2012), "Proliferation of External Globus Pallidus-Subthalamic Nucleus Synapses following Degeneration of Midbrain Dopamine Neurons"

            # 3. Modifications to GPE AMPA EPSCs (see hLTP)
            #   - DONE: see changes in cellpopdata.py
            #   - NMDA is involved in this hLTP mechanism

            # 4. Changes to regularity/variability of spontaneous firing (summarize literature)

            if state == PhysioState.PARK_DBS:
                # 5. Neurochemical effects of DBS?
                raise NotImplementedError()

        else:
            raise NotImplementedError('Unrecognized state %s' % repr(state))

        return soma_refs, dend_refs


    def _reduce_map_gillies(self, model):
        """
        Reduce Gillies model and map existing synapses to reduced model 
        """

        full_model = StnModel.Gillies2005

        # Make sure gillies model is built
        if not self.model_data[full_model].get('built', False):
            logger.info("Building Gillies original model first...")
            self.build_cell(full_model)
        else:
            logger.warn("Gillies STN model will be modified")

        ########################################################################
        # Gather full model info

        # Restore conductances
        gillies.reset_channel_gbar()

        # Get synapses and NetCon we want to map
        inputs = self.model_data[full_model]['inputs']
        pre_pops = inputs.keys()
        syns_tomap = self.get_synapse_data(full_model)

        ########################################################################
        # Reduce
        logger.debug("Starting model reduction...")
        
        # Create reduced cell
        # TODO: replace these by CERSEI reduction code, with option to use legacy code
        if model in (StnModel.Gillies_FoldMarasco_Legacy,
                        StnModel.Gillies_BranchZip_Legacy):
            from reducemodel import reduce_cell
            reduction = reduce_cell.gillies_marasco_reduction()
            reduction.set_syns_tomap(syns_tomap)
            reduction.reduce_model(num_passes=7)

        elif model == StnModel.Gillies_FoldStratford_Legacy:
            from reducemodel import reduce_cell
            reduction = reduce_cell.gillies_stratford_reduction()
            reduction.set_syns_tomap(syns_tomap)
            reduction.reduce_model(num_passes=1)

        elif model in (StnModel.Gillies_FoldBush_Tapered,
                        StnModel.Gillies_FoldMarasco_Tapered):
            import cersei_reduce
            model_to_method = {
                StnModel.Gillies_FoldBush_Tapered: 'BushSejnowski',
                StnModel.Gillies_FoldMarasco_Tapered: 'Marasco',
            }
            reduction = cersei_reduce.make_reduction(
                                        model_to_method[model],
                                        tweak=False)
            reduction.set_syns_tomap(syns_tomap)
            reduction.reduce_model(num_passes=1, map_synapses=True)

        else:
            return ValueError("Model not supported: {}".format(model))
            
        soma_refs, dend_refs = reduction._soma_refs, reduction._dend_refs

        ########################################################################
        # Save mapped inputs

        self._init_con_dict(model, pre_pops)
        new_inputs = self.model_data[model]['inputs']

        # Save mapped inputs
        for syn in reduction.map_syn_info:
            pre = syn.pre_pop.name.lower()
            new_inputs[pre]['synapses'].append(syn.mapped_syn)
            new_inputs[pre]['syn_NetCons'].extend(syn.afferent_netcons)

        return soma_refs, dend_refs


    def all_sec_refs(self, model):
        """
        Get SectionRef's to all sections in model.
        """
        return self.model_data[model]['soma_refs'] + self.model_data[model]['dend_refs']


    def reset_handlers(self, pre_pop):
        """
        Reset Initialize Handlers for given pre-synaptic input.
        """
        py_handlers = self.model_data[self.target_model]['inputs'][pre_pop].setdefault('PyInitHandlers', [])
        hoc_handlers = []

        for pyh in py_handlers:
            fih = h.FInitializeHandler(pyh)
            hoc_handlers.append(fih)

        self.model_data[self.target_model]['inputs'][pre_pop]['HocInitHandlers'] = hoc_handlers


    def print_synapse_info(self, model):
        """
        Print all synapse and NetCon parameter values
        """

        mech_param_names = {}
        
        for pre_pop, input_dict in self.model_data[model]['inputs'].items():
            
            for nc in input_dict['syn_NetCons']:
                
                syn = nc.syn()
                mech_name = cpd.get_mod_name(syn)

                # Load parameter names (once per mechanism)
                if mech_name not in mech_param_names:
                    mech_param_names[mech_name] = cpd.getSynMechParamNames(mech_name)

                print '\nInfo for synapse {} @ {}'.format(syn, syn.get_segment())
                print 'NetCon.weight[0] : {}'.format(nc.weight[0])

                for param_name in mech_param_names[mech_name]:
                    print '{} : {}'.format(param_name, getattr(syn, param_name, 'not found'))


    def get_synapse_data(self, model, pre_pops=None):
        """
        Get list of SynInfo properties containing a reference
        to each synapse, its NetCon and pre-synaptic population.

        @param  pre_pops    pre-synaptic populations for synapses. If none
                            are given, return synapses for all populations.

        @return             list(SynInfo)
        """

        inputs = self.model_data[model]['inputs']
        if pre_pops is None:
            pre_pops = inputs.keys()

        cc = cpd.CellConnector(self.physio_state, self.rng)

        all_syns = set((syn for pop in pre_pops for syn in inputs[pop].get('synapses', [])))
        all_ncs = set((nc for pop in pre_pops for nc in inputs[pop].get('syn_NetCons', [])))

        return cpd.get_synapse_data(cc, all_syns, all_ncs)
        
        # syn_list = []
        # for pop in pre_pops:

        #     # Get NetCon and unique synapses
        #     ncs = inputs[pop].get('syn_NetCons', [])
        #     syns = set(inputs[pop].get('synapses', [])) # unique Synapses

        #     # Get connection parameters
        #     # TODO: this should be saved when actually making the connections
        #     if not isinstance(pop, Populations):
        #         pre = Populations.from_descr(pop)
        #     else:
        #         pre = pop
        #     con_par = cc.getPhysioConParams(pre, Pop.STN, [Cit.Default]) 

        #     # Save properties
        #     for syn in syns:
        #         syn_info = mapsyn.SynInfo()

        #         # HocObjects
        #         syn_info.orig_syn = syn
        #         syn_info.afferent_netcons = [nc for nc in ncs if nc.syn().same(syn)]
                
        #         # meta-information
        #         syn_info.pre_pop = pop

        #         # For PSP frequency: need to find the receptor types that this synaptic mechanism implements
        #         modname = cpd.get_mod_name(syn)
        #         syn_info.mod_name = modname

        #         syn_receptors = cc.getSynMechReceptors(modname)
        #         freqs = [con_par[receptor]['f_med_PSP_burst'] for receptor in syn_receptors]
        #         syn_info.PSP_median_frequency = max(freqs)

        #         syn_list.append(syn_info)

        # return syn_list


    def add_inputs(self, pre_pop, model, **kwargs):
        """
        Add inputs to dict of existing inputs.
        """
        pop_inputs = self.model_data[model]['inputs'].get(pre_pop, {})

        for input_type, input_objs in kwargs.items():
            if input_type in pop_inputs:
                # Add to list (don't overwrite)
                pop_inputs[input_type].extend(input_objs)
            else:
                pop_inputs[input_type] = input_objs

        self.model_data[model]['inputs'][pre_pop] = pop_inputs


    def get_inputs(self, pre_pop, model):
        """
        Get inputs from given pre-synaptic population.

        @return     dictionary containing (optional) keys:
                        'synapses'
                        'syn_NetCons'
                        'NetStims'
                        'RNG_data'
                        ...
        """
        if isinstance(pre_pop, cpd.Populations):
            pre_pop = pre_pop.name.lower()
        
        return self.model_data[model]['inputs'].get(pre_pop, None)


    def get_all_inputs(self, model):
        """
        Get inputs from all pre-synaptic cells.
        """
        pre_pops = self.model_data[model]['inputs'].keys()
        return self.merged_inputs(pre_pops, model)


    def merged_inputs(self, pops, model):
        """
        Return input dict with all data from multiple pre-synaptic populations
        merged into one dict.
        """
        merged = {}
        for pop in pops:
            if not isinstance(pop, str):
                pop = pop.name.lower()
            if not pop in self.model_data[model]['inputs']:
                logger.warning('No inputs found for pre-synaptic population {}'.format(pop))
                continue
            for input_type, inputs in self.model_data[model]['inputs'][pop].items():
                extend_dictitem(merged, input_type, inputs, append=False)
        return merged


    def save_input_data():
        pass


    def load_input_data():
        """
        Restore all NetStim, NetCon, Synapse, Random objects from pickle file.
        """
        pass


    def _init_con_dict(self, model, pre_pops):
        """
        Initialize dict with connection data for given model
        and presynaptic populations.
        """
        for pop in pre_pops:
            self.model_data[model]['inputs'][pop] = {
                'stimweightvec': [],
                'synapses': [],
                'syn_NetCons': [], # NetCons targetting synapses
                'com_NetCons': [], # NetCons for control events
                'NetStims': [],
                'HocInitHandlers': [],
                'PyInitHandlers': [],
                'RNG_data': [],
            }


    def get_pre_pop(self, syn, model):
        """
        Get pre-synaptic population for given synapse.

        @return     str: abbreviation for pre-synaptic population
        """
        inputs = self.model_data[model]['inputs']

        if 'NetCon' in syn.hname(): # dirty hack: cast to NetCon if named 'NetCon[xyz]'
            gen_pop = (pop for (pop, conn_data) in inputs.items() if syn in conn_data['syn_NetCons'])
        else:
            gen_pop = (pop for (pop, conn_data) in inputs.items() if syn in conn_data['synapses'])
        
        return next(gen_pop, None)


    def get_num_syns(self, model):
        """
        Get total number of synapses that have been creates on given model.
        """
        num_syn = 0
        for pop in list(Populations):
            pops_inputs = self.model_data[model]['inputs']
            if pop in pops_inputs.keys():
                num_syn += len(pops_inputs[pop]['synapses'])
        
        return num_syn


    def _add_recorded_obj(self, tag, obj, protocol, model=None):
        """
        Add given object to list of recorded objects for given model & protocol.
        """
        if model is None:
            model = self.target_model

        rec_objs = self.model_data[model]['rec_segs'][protocol]
        
        # Check if tag already exists
        prev_obj = rec_objs.get(tag, None)
        if prev_obj is not None:
            logger.warning("A recorded object with name {} already exists. Overwriting.".format(tag))
        
        rec_objs[tag] = obj


    ############################################################################
    # SETUP functions
    ############################################################################

    def make_inputs(self, stim_protocol):
        """
        Make the inputs for given stimulation protocol.
        """
        if ((self.target_model in cpd.ReducedModels) and 
            (stim_protocol in SynapticProtocols)):
            
            logger.info("Skip making inputs for model. Inputs for reduced model "
                        "should have been mapped from full model.")
            return

        cc = cpd.CellConnector(self.physio_state, self.rng)

        if stim_protocol == StimProtocol.SPONTANEOUS:
            # Spontaneous firing has no inputs
            pass

        elif stim_protocol == StimProtocol.SYN_PARK_PATTERNED:
            
            # Method 1:
            #   - add netstim and play input signal into weight (see Dura-Berndal arm example)
            #   - advantage: low weight in off-period will provide background noise
            #   - if no input desired in off-period: periodically turn it on/off using events or other NetStim
            
            # Method 2:
            #   - use nsloc.mod (=netstim with variable rate)

            ####################################################################
            # CTX inputs
            ####################################################################

            # TODO: calibrate cortical inputs
            #   - test that EPSP/IPSP have desired size
            #   - in DA-depleted state: should trigger bursts

            n_ctx_syn = 10
            new_inputs = {}


            # Make weight signal representing oscillatory pattern
            ctx_timevec = np.arange(0, self.sim_dur, 0.05) # update every 0.05 ms
            ctx_pattern_freq = 8.0 # frequency [Hz]
            ctx_pattern_phase = 0.0 # anti-phase from GPE
            ctx_radvec = 2*np.pi*ctx_pattern_freq*1e-3*ctx_timevec + ctx_pattern_phase
            
            # Set ON and OFF phase
            duty_ms = 50.0 # max is 1/freq * 1e3
            duty = duty_ms / (1./ctx_pattern_freq*1e3)
            ctx_pattern = signal.square(ctx_radvec, duty)
            ctx_pattern[ctx_pattern<0.0] = 0.05 # fractional noise amplitude

            stimweightvec = h.Vector(ctx_pattern)
            stimtimevec = h.Vector(ctx_timevec)

            # Distribute synapses over dendritic trees
            is_ctx_target = lambda seg: seg.diam <= 1.0         
            dend_secrefs = self.model_data[self.target_model]['dend_refs']
            ctx_target_segs = pick_random_segments(dend_secrefs, n_ctx_syn, is_ctx_target, rng=self.rng)

            # Make synapses
            for target_seg in ctx_target_segs:

                # Make poisson spike generator
                stim_rate = 80.0
                stimsource = h.NetStim() # Create a NetStim
                stimsource.interval = stim_rate**-1*1e3 # Interval between spikes
                stimsource.number = 1e9 # max number of spikes
                stimsource.noise = 0.25 # Fractional noise in timing
                # stimsource.noiseFromRandom(stimrand) # Set it to use this random number generator

                # TODO SETPARAM: (CTX) set poisson noise & rate parameters dependent on PhysioState
                #   use references for this (reported in vivo firing rates, traces etc)

                # Make synapse and NetCon
                syn, nc, wvecs = cc.make_synapse((Pop.CTX, Pop.STN), (stimsource, target_seg), 
                                    'GLUsyn', (NTR.AMPA, NTR.NMDA), use_sources=(Cit.Chu2015,), 
                                    weight_scales=[stimweightvec], weight_times=[stimtimevec])

                # Save inputs
                extend_dictitem(new_inputs, 'synapses', syn)
                extend_dictitem(new_inputs, 'syn_NetCons', nc)
                extend_dictitem(new_inputs, 'NetStims', stimsource)
                extend_dictitem(new_inputs, 'stimweightvec', wvecs)

            # Save inputs to model
            extend_dictitem(new_inputs, 'stimweightvec', stimweightvec)
            extend_dictitem(new_inputs, 'stimtimevec', stimtimevec)
            self.add_inputs('ctx', self.target_model, **new_inputs)

            ####################################################################
            # GPe inputs
            ####################################################################

            # Add GPe inputs using Tsodyks-Markram synapses
            n_gpe_syn = 10 # NOTE: one synapse represents a multi-synaptic contact from one GPe axon
            new_inputs = {}

            # Make weight signal representing oscillatory pattern
            gpe_timevec = np.arange(0, self.sim_dur, 0.05) # update every 0.05 ms
            gpe_pattern_freq = 8.0 # frequency [Hz]
            gpe_pattern_phase = np.pi # anti-phase from CTX
            gpe_radvec = 2*np.pi*gpe_pattern_freq*1e-3*gpe_timevec + gpe_pattern_phase
            
            # Set ON and OFF phase
            duty_ms = 80.0 # max is 1/freq * 1e3
            duty = duty_ms / (1./gpe_pattern_freq*1e3)
            gpe_pattern = signal.square(gpe_radvec, duty)
            gpe_pattern[gpe_pattern<0.0] = 0.05 # fractional noise amplitude

            stimweightvec = h.Vector(gpe_pattern)
            stimtimevec = h.Vector(gpe_timevec)

            # Pick random segments in dendrites for placing synapses
            is_gpe_target = lambda seg: seg.diam > 1.0 # select proximal dendrites
            dend_secrefs = self.model_data[self.target_model]['dend_refs']
            gpe_target_segs = pick_random_segments(dend_secrefs, n_gpe_syn, is_gpe_target, rng=self.rng)

            # Make synapses
            for target_seg in gpe_target_segs:

                # Make poisson spike generator
                stim_rate = 100.0 # hz
                stimsource = h.NetStim() # Create a NetStim
                stimsource.interval = stim_rate**-1*1e3 # Interval between spikes
                stimsource.number = 1e9 # max number of spikes
                stimsource.noise = 0.25 # Fractional noise in timing
                # stimsource.noiseFromRandom(stimrand) # Set it to use this random number generator

                # TODO SETPARAM: (GPe) set poisson noise & rate parameters, dependent on PhysioState
                # HallworthBevan2005_DynamicallyRegulate: 
                #       Fig 6: 10 pulses @ 100 Hz and 20 pulses @ 100 Hz triiger rebound bursts
                # Bevan2006_CellularPrinciples:
                #       Fig 2: 10 stimuli @ 100 Hz trigger pause + rebound burst

                # Make synapse and NetCon
                syn, nc, wvecs = cc.make_synapse((Pop.GPE, Pop.STN), (stimsource, target_seg), 
                                    'GABAsyn', (NTR.GABAA, NTR.GABAB), 
                                    use_sources=(Cit.Chu2015, Cit.Fan2012, Cit.Atherton2013), 
                                    weight_scales=[stimweightvec], weight_times=[stimtimevec])

                # Save inputs
                extend_dictitem(new_inputs, 'synapses', syn)
                extend_dictitem(new_inputs, 'syn_NetCons', nc)
                extend_dictitem(new_inputs, 'NetStims', stimsource)
                extend_dictitem(new_inputs, 'stimweightvec', wvecs)

            # Save inputs to model
            extend_dictitem(new_inputs, 'stimtimevec', stimtimevec)
            extend_dictitem(new_inputs, 'stimweightvec', stimweightvec)
            self.add_inputs('gpe', self.target_model, **new_inputs)

        else: # standard action: look up in dict

            try:
                make_inputs_func = self._MAKE_INPUT_FUNCS[stim_protocol]
            except KeyError:
                raise NotImplementedError("Make inputs function for protocol {} not implemented".format(stim_protocol))

            make_inputs_func(self, connector=cc)


    def map_inputs(self, cand_model):
        """
        Map inputs from target model to candidate model.
        """
        raise NotImplementedError()

    def rec_GABA_traces(self, protocol, traceSpecs, n_syn=1):
        """
        Set up recording Vectors

        @param n_syn        number of synaptic traces to record
        """

        rec_segs = self.model_data[self.target_model]['rec_segs'][protocol]
        model = self.target_model
        
        # Add synapse and segment containing it
        nc_list = self.model_data[model]['inputs']['gpe']['syn_NetCons']
        for i_syn, nc in enumerate(nc_list):
            if i_syn > n_syn-1:
                break

            syn_tag = 'GABAsyn%i' % i_syn
            seg_tag = 'GABAseg%i' % i_syn
            rec_segs[syn_tag] = nc.syn()
            rec_segs[seg_tag] = nc.syn().get_segment()

            # Record synaptic variables
            traceSpecs['gA_GABAsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'g_GABAA'}
            traceSpecs['gB_GABAsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'g_GABAB'}
            traceSpecs['Rrp_GABAsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'Rrp'}
            traceSpecs['Use_GABAsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'Use'}
            traceSpecs['Hill_GABAsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'G'}


    def rec_GLU_traces(self, protocol, traceSpecs, n_syn=1):
        """
        Set up recording Vectors
        """
        rec_segs = self.model_data[self.target_model]['rec_segs'][protocol]
        model = self.target_model
        
        # Add synapse and segment containing it
        nc_list = self.model_data[model]['inputs']['ctx']['syn_NetCons']
        for i_syn, nc in enumerate(nc_list):
            if i_syn > n_syn-1:
                break

            syn_tag = 'GLUsyn%i' % i_syn
            seg_tag = 'GLUseg%i' % i_syn
            rec_segs[syn_tag] = nc.syn()
            rec_segs[seg_tag] = nc.syn().get_segment()

            # Record synaptic variables
            traceSpecs['gA_GLUsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'g_AMPA'}
            traceSpecs['gN_GLUsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'g_NMDA'}
            traceSpecs['Rrp_GLUsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'R'}
            traceSpecs['Use_GLUsyn%i' % i_syn] = {'pointp':syn_tag, 'var':'Use'}


    def rec_Vm(self, protocol, traceSpecs):
        """
        Record membrane voltages in all recorded segments
        """
        rec_segs = self.model_data[self.target_model]['rec_segs'][protocol]
        
        for seclabel, seg in rec_segs.items():
            if isinstance(seg, neuron.nrn.Segment):
                traceSpecs['V_'+seclabel] = {'sec':seclabel, 'loc':seg.x, 'var':'v'}


    def rec_traces(self, stim_protocol, recordStep=0.025):
        """
        Set up recording Vectors to record from relevant pointers
        """
        # Initialize data
        model = self.target_model
        self.model_data[model]['rec_data'][stim_protocol] = {}

        # Specify sections to record from
        if model == StnModel.Gillies2005:

            somasec = h.SThcell[0].soma
            dendsec = h.SThcell[0].dend1[7]

            # Assign label to each recorded section
            rec_segs = {
                'soma': somasec(0.5), # middle of soma
                'dist_dend': dendsec(0.8), # approximate location along dendrite in fig. 5C
            }

        elif model in cpd.ReducedModels:

            somasec = self.model_data[model]['soma_refs'][0].sec
            dendrefs = self.model_data[model]['dend_refs']
            dendsec = next(ref.sec for ref in dendrefs if not any(ref.sec.children()))

            # Default recorded segments
            rec_segs = {
                'soma': somasec(0.5), # middle of soma
                'dist_dend': dendsec(0.9), # approximate location along dendrite in fig. 5C
            }

        else:
            raise NotImplementedError("""Recording from other models 
                    besides {} not yet implemented""".format(StnModel.Gillies2005))

        # Save recorded segments list
        self.model_data[model]['rec_segs'][stim_protocol] = rec_segs

        # Start trace specification
        traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)
        traceSpecs['t_global'] = {'var':'t'}
        self.rec_dt = recordStep

        # PROTOCOL-SPECIFIC TRACES
        if stim_protocol == StimProtocol.SPONTANEOUS: # spontaneous firing (no inputs)
            
            # Trace specs for membrane voltages
            for seclabel, seg in rec_segs.items():
                traceSpecs['V_'+seclabel] = {'sec':seclabel, 'loc':seg.x, 'var':'v'}

            # Trace specs for recording ionic currents, channel states
            analysis.rec_currents_activations(traceSpecs, 'soma', 0.5)
            

        elif stim_protocol == StimProtocol.SYN_PARK_PATTERNED: # pathological input, strong patterned cortical input with strong GPi input in antiphase
            ####################################################################
            # Record SYN_PARK_PATTERNED
            ####################################################################

            # See diagram in marasco_reduction.pptx
            # dist_dend0_ids = [8,9,7,10,12,13,18,19,17,20,22,23]
            # prox_dend0_ids = [2,4,5,3,14,15]
            # dist_dend1_ids = [6,7,5,8,10,11]
            # prox_dend1_ids = [1,2,3]

            # Pick some proximal and distal dendritic sections for recording
            dist_secs = [(0,9), (0,10), (0,17), (0,23), (1,6), (1,8)]
            prox_secs = [(0,2), (0,3), (1,1)]

            # Add to recorded sections
            for tree_id, sec_id in dist_secs:
                tree_name = 'dend' + str(tree_id)
                sec = getattr(h.SThcell[0], tree_name)[sec_id-1]
                rec_segs['dist_' + repr(sec)] = sec(0.9)

            for tree_id, sec_id in prox_secs:
                tree_name = 'dend' + str(tree_id)
                sec = getattr(h.SThcell[0], tree_name)[sec_id-1]
                rec_segs['prox_' + repr(sec)] = sec(0.9)

            # Pick some segments that are targeted by synapse
            gpe_ncs = self.model_data[self.target_model]['inputs']['gpe']['syn_NetCons']
            ctx_ncs = self.model_data[self.target_model]['inputs']['ctx']['syn_NetCons']
            gpe_picks = [gpe_ncs[i] for i in self.rng.choice(len(gpe_ncs), 3, replace=False)]
            ctx_picks = [ctx_ncs[i] for i in self.rng.choice(len(ctx_ncs), 3, replace=False)]
            for nc in gpe_picks + ctx_picks:
                seg = nc.syn().get_segment()
                rec_segs['postsyn_' + repr(seg.sec)] = seg


            # Specify which traces you want in these sections
            for seclabel, seg in rec_segs.items():
                # Membrane voltages
                traceSpecs['V_'+seclabel] = {'sec':seclabel, 'loc':seg.x, 'var':'v'}

                # # K currents (dendrite)
                # traceSpecs['I_KCa_'+seclabel] = {'sec':'dist_dend','loc':seg.x,'mech':'sKCa','var':'isKCa'}
                
                # # Ca currents (dendrite)
                # traceSpecs['I_CaL_'+seclabel] = {'sec':'dist_dend','loc':seg.x,'mech':'HVA','var':'iLCa'}
                # traceSpecs['I_CaN_'+seclabel] = {'sec':'dist_dend','loc':seg.x,'mech':'HVA','var':'iNCa'}
                # traceSpecs['I_CaT_'+seclabel] = {'sec':'dist_dend','loc':seg.x,'mech':'CaT','var':'iCaT'}

        else: # standard action: look up in dict

            try:
                rec_trace_func = self._REC_TRACE_FUNCS[stim_protocol]
            except KeyError:
                raise NotImplementedError("Recording function for protocol {} not implemented".format(stim_protocol))

            rec_trace_func(self, stim_protocol, traceSpecs)



        # Prepare dictionary (label -> Section)
        rec_secs = {}
        for seclabel, hobj in rec_segs.items():
            if isinstance(hobj, neuron.nrn.Segment):
                rec_secs[seclabel] = hobj.sec
            else:
                rec_secs[seclabel] = hobj # point process

        # Use trace specs to make Hoc Vectors
        recData, markers = analysis.recordTraces(rec_secs, traceSpecs, recordStep)

        # Save trace specs and recording Vectors
        self.model_data[model]['rec_data'][stim_protocol].update({
            'trace_specs': traceSpecs,
            'trace_data': recData,
            'rec_dt': recordStep,
            'rec_markers': markers,
        })

    ############################################################################
    # PLOT functions
    ############################################################################

    def save_proto_traces(self, model, protocol, filename):
        """
        Save recorded traces for given model and stimulation protocol.
        """
        rec_data = dict(self.model_data[model]['rec_data'][protocol]) # copy it
        
        # remove unnecessary data
        rec_data.pop('rec_markers') # can't pickle/deepcopy HocObject
        rec_data = deepcopy(rec_data) # don't modify original trace data

        # convert h.Vector to numpy array
        trace_data = rec_data['trace_data']
        for k,v in trace_data.items():
            trace_data[k] = v.as_numpy()

        # Write arrays to npz file
        with open(filename, 'w') as recfile:
            pickle.dump(rec_data, recfile)


    def _plot_all_Vm(self, model, protocol, fig_per='trace'):
        """
        Plot all membrane voltages.
        """
        # Get data
        rec_dict = self.model_data[model]['rec_data'][protocol]
        recData, recordStep = (rec_dict[k] for k in ('trace_data', 'rec_dt'))

        return proto_common.plot_all_Vm(recData, 
                                recordStep=recordStep,
                                fig_per=fig_per)
    

    def _plot_all_spikes(self, model, protocol, **kwargs):
        """
        Plot all recorded spikes.

        @pre        all traces containing spike times have been tagged
                    with prefix 'AP_'. If not, provide a custom filter
                    function in param 'trace_filter'

        @param trace_filter     filter function for matching spike traces.

        @param kwargs           can be used to pass any arguments of analysis.plotRaster()
        """
        # Get data
        rec_dict = self.model_data[model]['rec_data'][protocol]
        recData, rec_dt = (rec_dict[k] for k in ('trace_data', 'rec_dt'))

        return proto_common.plot_all_spikes(recData, **kwargs)


    def _plot_GABA_traces(self, model, protocol, fig_per='cell'):
        """
        Plot GABA synapse traces.
        """
        # Get data
        rec_dict = self.model_data[model]['rec_data'][protocol]
        recData, recordStep = (rec_dict[k] for k in ('trace_data', 'rec_dt'))

        return proto_common.plot_GABA_traces(
                                recData, 
                                recordStep=recordStep, 
                                fig_per=fig_per)


    def _plot_GLU_traces(self, model, protocol, fig_per='cell'):
        """
        Plot GABA synapse traces.
        """
        # Get data
        rec_dict = self.model_data[model]['rec_data'][protocol]
        recData, recordStep = (rec_dict[k] for k in ('trace_data', 'rec_dt'))

        return proto_common.plot_GLU_traces(
                                recData, 
                                recordStep=recordStep, 
                                fig_per=fig_per)


    ############################################################################
    # SIMULATION functions
    ############################################################################

    def plot_traces(self, protocol, model=None):
        """
        Plot relevant recorded traces for given protocol
        """

        # Get recorded data
        if model is None:
            model = self.target_model

        recData = self.model_data[model]['rec_data'][protocol]['trace_data']
        recordStep = self.model_data[model]['rec_data'][protocol]['rec_dt']

        # Plot membrane voltages
        

        # Extra plots depending on simulated protocol
        if protocol == StimProtocol.SPONTANEOUS:

            self._plot_all_Vm(model, protocol)

            # Plot ionic currents, (in)activation variables
            figs, cursors = analysis.plot_currents_activations(self.recData, recordStep)

        elif protocol == StimProtocol.SYN_PARK_PATTERNED:

            V_prox = analysis.match_traces(recData, lambda t: t.startswith('V_prox'))
            V_dist = analysis.match_traces(recData, lambda t: t.startswith('V_dist'))
            V_postsyn = analysis.match_traces(recData, lambda t: t.startswith('V_postsyn'))

            analysis.plotTraces(V_prox, recordStep, yRange=(-80,40), traceSharex=True, 
                                title='Proximal sections')
            analysis.plotTraces(V_dist, recordStep, yRange=(-80,40), traceSharex=True,
                                title='Distal sections')
            analysis.plotTraces(V_postsyn, recordStep, yRange=(-80,40), traceSharex=True,
                                title='Post-synaptic segments')
        
        
        else: # standard action: look up in dict

            try:
                plot_trace_func = self._PLOT_TRACE_FUNCS[protocol]  
            except KeyError:
                raise NotImplementedError("Plotting function for protocol {} not implemented".format(protocol))

            plot_trace_func(self, model, protocol)


    def run_sim(self, nthread=1):
        """
        Run NEURON simulator for `dur` or `self.sim_dur` milliseconds
        with precise measurement of runtime
        """
        # enable multithreaded execution
        if nthread > 1:
            h.cvode_active(0)
            h.load_file('parcom.hoc')
            pct = h.ParallelComputeTool[0]
            pct.nthread(nthread)
            pct.multisplit(1)
            pct.busywait(1)

        # Simulate
        logger.debug("Simulating...")
        t0 = h.startsw()
        h.run()
        t1 = h.startsw()
        h.stopsw() # or t1=h.startsw(); runtime = t1-t0
        logger.debug("Simulated for {:.6f} seconds".format(t1-t0))


    def init_sim(self, stim_protocol, dur=2000., dt=0.025, celsius=35., v_init=-60):
        """
        Initialize simulation.
        """
        try:
            init_sim_func = self._INIT_SIM_FUNCS[stim_protocol]
        
        except KeyError:
            logger.warning("Simulator initializaton function for protocol {} not implemented. Falling back to default initialization function.".format(stim_protocol))

            self._init_sim()

        else:
            init_sim_func(self, stim_protocol)

    
    def _init_sim(self, dur=2000.0, dt=0.025, celsius=35., v_init=-60):
        """
        Default initialization function.

        Sets simulation duration, time step, temperature, 
        initial membrane voltage and ion concentrations.
        """

        self.sim_dur = dur
        h.tstop = dur
        
        self.sim_dt = dt
        h.dt = dt

        h.celsius = celsius # different temp from paper
        h.v_init = v_init # paper simulations sue default v_init
        gillies.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

        # Initialize NEURON simulator
        h.init() # calls finitialize()


    def _setup_proto(self, proto):
        """
        Standard setup function.

        Dispatches protocol-specific setup functions.
        """
        # Make inputs
        self.make_inputs(proto)

        # Set up recording
        self.rec_traces(proto, recordStep=0.05)


    def _run_proto(self, proto, stdinit=False):
        """
        Standard simulation function.

        Dispatches protocol-specific init & run functions.
        """

        # Initialize
        if stdinit:
            h.stdinit()
        else:
            self.init_sim(proto)

        # Simulate
        self.run_sim()


    def register_keyword_funcs(
            self,
            proto_init_funcs=None,
            proto_setup_funcs_pre=None,
            proto_setup_funcs_post=None,
            proto_setup_kwargs_const=None,
            proto_setup_kwargs_getters=None,
            proto_setup_kwargs_setters=None,
            proto_rec_funcs=None,
            proto_plot_funcs=None,
            reset=False
        ):
        """
        Register functions that setup a stimulation protocol using functions
        that accept only keyword arguments (for compatibility with Ephys
        SelfContainedProtocol).
        """
        # Copy each container to avoid later modification of pass-by-reference containers
        proto_init_funcs = [] if proto_init_funcs is None else list(proto_init_funcs)
        proto_rec_funcs = [] if proto_rec_funcs is None else list(proto_rec_funcs)
        proto_plot_funcs = [] if proto_plot_funcs is None else list(proto_plot_funcs)
        
        proto_setup_funcs_pre = [] if proto_setup_funcs_pre is None else list(proto_setup_funcs_pre)
        proto_setup_funcs_post = [] if proto_setup_funcs_post is None else list(proto_setup_funcs_post)
        
        proto_setup_kwargs_const = {} if proto_setup_kwargs_const is None else dict(proto_setup_kwargs_const)
        proto_setup_kwargs_getters = collections.OrderedDict() if proto_setup_kwargs_getters is None else collections.OrderedDict(proto_setup_kwargs_getters)
        proto_setup_kwargs_setters = [] if proto_setup_kwargs_setters is None else list(proto_setup_kwargs_setters)

        if reset:
            self._proto_init_funcs = proto_init_funcs
            self._proto_setup_funcs_pre = proto_setup_funcs_pre
            self._proto_setup_funcs_post = proto_setup_funcs_post
            self._proto_setup_kwargs_const = proto_setup_kwargs_const
            self._proto_setup_kwargs_getters = proto_setup_kwargs_getters
            self._proto_setup_kwargs_setters = proto_setup_kwargs_setters
            self._proto_rec_funcs = proto_rec_funcs
            self._proto_plot_funcs = proto_plot_funcs
        else:
            self._proto_init_funcs.extend(proto_init_funcs)
            self._proto_setup_funcs_pre.extend(proto_setup_funcs_pre)
            self._proto_setup_funcs_post.extend(proto_setup_funcs_post)
            self._proto_setup_kwargs_const.update(proto_setup_kwargs_const)
            self._proto_setup_kwargs_getters.update(proto_setup_kwargs_getters)
            self._proto_setup_kwargs_setters.extend(proto_setup_kwargs_setters)
            self._proto_rec_funcs.extend(proto_rec_funcs)
            self._proto_plot_funcs.extend(proto_plot_funcs)


    def setup_keyword_protocol(self, pre_model=None, post_model=None):
        """
        Set up stimulation protocols using all keyword functions registered
        until now.

        @note   call after self.build_cell()

        @note   call self.run_sim() after this function

        @note   if you want to run full model first, only provide pre_model.
                Then after running you can call register_keyword_funcs(..., reset=True)
                and call this function again with only post_model

        @post   The dictionary with all shared data for protocol is stored
                as self._iproto_data
        """
        if pre_model is not None:

            # Prepare inputs
            self.build_cell(pre_model)
            model = pre_model
            ICell = collections.namedtuple("ICell", ['dendritic', 'somatic'])
            icell = ICell(
                    dendritic=[ref.sec for ref in self.model_data[model]['dend_refs']],
                    somatic=[ref.sec for ref in self.model_data[model]['soma_refs']])

            # Create keyword arguments for protocol setup functions
            kwargs_default = {
                'nrnsim': h,
                'stim_data':    {}, # synapses and netcons
                'evaluator':    self,
                'gid':          self.model_data[model]['gid'],
                'icell':        icell,
                'rng_info':     self.rng_info,
            }

            # Assemble all keyword arguments
            logger.debug("Instantiating protocol setup kwargs...")
            self._iproto_data = kwargs_default
            self._iproto_data.update(self._proto_setup_kwargs_const)

            for kwarg_name, kwarg_getter in self._proto_setup_kwargs_getters.items():
                self._iproto_data[kwarg_name] = kwarg_getter(self._iproto_data)
            for kwarg_setter in self._proto_setup_kwargs_setters:
                kwarg_setter(self._iproto_data)

            # Set up stimulators, synaptic inputs
            for setup_func in self._proto_setup_funcs_pre:
                setup_func(**self._iproto_data)
        
        # POST MODEL (reduced) setup funcs
        if post_model is not None:
            self.build_cell(post_model)
            model = post_model
            post_icell = ICell(dendritic=[sec for sec in self.model_data[model]['dend_refs']])
            self._iproto_data['icell'] = post_icell

            for setup_func in self._proto_setup_funcs_post:
                setup_func(**self._iproto_data)

        # Prepare recording data containers
        iproto_rec_data = {
            'icell': icell,
            'trace_specs': collections.OrderedDict(),
            'trace_vectors': {},
            'rec_hoc_markers': {},
        }
        self._iproto_data.update(iproto_rec_data)
        # params that may have been set by kwarg getter
        self._iproto_data.setdefault('rec_hoc_objects', {})
        self._iproto_data.setdefault('record_step', 0.05)

        # Custom recording functions
        for rec_func in self._proto_rec_funcs:
            rec_func(**self._iproto_data) # updates trace_specs
            logger.debug("Executed recording function {}".format(rec_func))

        # Create recording vectors
        if not 't_global' in self._iproto_data['trace_specs']:
            self._iproto_data['trace_specs']['t_global'] = {'var':'t'}
        
        recData, markers = analysis.recordTraces(
                                self._iproto_data['rec_hoc_objects'], 
                                self._iproto_data['trace_specs'],
                                recordStep=self._iproto_data['record_step'])

        self._iproto_data['trace_vectors'] = recData
        self._iproto_data['rec_hoc_markers'] = markers

        # Initialize physiological conditions
        for init_func in self._proto_init_funcs:
            init_func(**self._iproto_data)


    def run_keyword_protocol(self, plot_traces=True):
        """
        Run stimulation protocols that was set up using keyword functions.
        """
        self.run_sim()
        if plot_traces:
            remaining_kwargs = dict(self._iproto_data)
            trace_vecs = remaining_kwargs.pop('trace_vectors')
            for plot_func in self._proto_plot_funcs:
                plot_func(trace_vecs, **remaining_kwargs)


    def setup_run_protocol(self, protocol, model=None):
        """
        Simulate cell in physiological state, under given stimulation protocol.

        @param model    Model to simulate protocol with. If no model is given,
                        the target model is used.

        @pre        cell must be built using build_cell

        @pre        physiological state must be set

        @effect     makes inputs (afferents, stimulation) based on physiological 
                    state and stimulation protocol (see make_inputs)

        @effect     runs a simulation
        """
        
        if protocol == StimProtocol.SPONTANEOUS: # spontaneous firing (no inputs)
            
            # Set simulation parameters
            self.sim_dur = 1000
            h.dt = 0.025
            self.sim_dt = h.dt

            h.celsius = 35 # different temp from paper (fig 3B: 25degC, fig. 3C: 35degC)
            h.v_init = -60 # paper simulations use default v_init
            gillies.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

            # Set up recording
            self.rec_traces(protocol, recordStep=0.05)

            # Simulate
            h.init()
            self.run_sim()

        else:
            # Standard simulation function
            self._setup_proto(protocol) # MAKE_INPUTS + RECORD_TRACES
            self._run_proto(protocol) # INIT_SIMULATION + RUN_SIMULATION


################################################################################
# EXPERIMENTS
################################################################################

def index_or_name(trace_data):
    """
    Ordering function: gets trace name or index.
    """
    name = trace_data[0]
    m = re.search(r'\[([\d]+)\]', name) # get index [i]
    key = int(m.groups()[0]) if m else name
    return key


def run_protocol(proto, model, export_locals=True):
    """
    Run given stimulation protocol.
    """

    # Make cell model and evaluator
    evaluator = StnModelEvaluator(model, PhysioState.NORMAL)
    evaluator.build_cell(model)
    
    # Run protocol
    evaluator.setup_run_protocol(proto)
    evaluator.plot_traces(proto)

    if export_locals:
        globals().update(locals())

    return evaluator


def map_protocol(proto, model, export_locals=True, pause=False):
    """
    Run given stimulation protocol.

    @param  model       the reduced model you want to test
    """
    full_model = StnModel.Gillies2005
    red_model = model

    # Make cell model and evaluator
    evaluator = StnModelEvaluator(full_model, PhysioState.NORMAL)
    evaluator.build_cell(full_model)
    
    # Run protocol
    evaluator.setup_run_protocol(proto)
    evaluator.plot_traces(proto)

    # Inspect model before continuing
    if pause:
        from neuron import gui
        raw_input("Use ModelView to inspect model. Press enter to continue.")
    ##################################################

    # Model reduction
    evaluator.build_cell(red_model)
    evaluator.target_model = red_model

    # Run Protocol
    evaluator.setup_run_protocol(proto)
    evaluator.plot_traces(proto)

    if export_locals:
        globals().update(locals())


if __name__ == '__main__':
    # map_protocol_MIN_SYN_BURST()
    # map_protocol_SYN_BACKGROUND_HIGH()
    # map_protocol(StimProtocol.CLAMP_REBOUND, StnModel.Gillies_FoldMarasco)
    # map_protocol(StimProtocol.CLAMP_REBOUND, StnModel.Gillies_FoldStratford)
    map_protocol(StimProtocol.SYN_BACKGROUND_HIGH, StnModel.Gillies_FoldMarasco)
    # map_protocol(StimProtocol.PASSIVE_SYN, StnModel.Gillies_FoldMarasco, pause=True)
