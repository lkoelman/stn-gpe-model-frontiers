"""
Evaluation of cell models: setting up experimental protocols, recordings,
and plotting.
"""

# Standard library
import collections

# Third party
import numpy as np

from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# Plotting & recording
from . import analysis

# Physiological parameters
import bgcellmodels.cellpopdata as cpd
from bgcellmodels.cellpopdata import PhysioState, Populations
Pop = Populations

# Logging
import logging
logging.basicConfig(format='%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logger = logging.getLogger('protocols') # create logger for this module

# Experimental protocols
from .stimprotocols import StimProtocol, extend_dictitem


class CellEvaluator(object):
    """
    Evaluate Cell models
    """


    def __init__(self, model_id, physio_state=PhysioState.NORMAL):
        """
        Initialize new evaluator in given physiological state.
        """
        self._physio_state = physio_state

        self.model_data = {
            model_id: {
                'rec_data' : {proto: {} for proto in StimProtocol},
                'rec_segs' : {proto: {} for proto in StimProtocol},
                'inputs' : {}, # key for each pre-synaptic population,
                'user_params' : {}, # parameters for model building
                'proto_vars' : {proto: {} for proto in StimProtocol},
            }
        }
        self.target_model = model_id

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

        # Save Sections
        # TODO: pass in _init_ or let uses make custom build func
        self.model_data[model]['soma_refs'] = TODO
        self.model_data[model]['dend_refs'] = TODO

        # Indicate that given model has been built
        self.model_data[model]['built'] = True
        self.model_data[model]['gid'] = 1 # we only have one cell


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
        
        for pre_pop, input_dict in self.model_data[model]['inputs'].iteritems():
            
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


    def add_inputs(self, pre_pop, model, **kwargs):
        """
        Add inputs to dict of existing inputs.
        """
        pop_inputs = self.model_data[model]['inputs'].get(pre_pop, {})

        for input_type, input_objs in kwargs.iteritems():
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
            for input_type, inputs in self.model_data[model]['inputs'][pop].iteritems():
                extend_dictitem(merged, input_type, inputs, append=False)
        return merged


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
            gen_pop = (pop for (pop, conn_data) in inputs.iteritems() if syn in conn_data['syn_NetCons'])
        else:
            gen_pop = (pop for (pop, conn_data) in inputs.iteritems() if syn in conn_data['synapses'])
        
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
    # SIMULATION functions
    ############################################################################


    def run_sim(self, nthread=1):
        """
        Run NEURON simulator for `dur` or `self.sim_dur` milliseconds
        with precise measurement of runtime

        @note   this method does NOT call h.init()/h.stdinit()/
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

        # Initialize NEURON simulator
        h.init() # calls finitialize()


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

            for kwarg_name, kwarg_getter in self._proto_setup_kwargs_getters.iteritems():
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
