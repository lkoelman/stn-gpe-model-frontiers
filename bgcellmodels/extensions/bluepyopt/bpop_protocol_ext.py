"""
Extension of BluePyOpt protocol that allow instantiating stimuli, synapses, etc.
without specifying them in the form of ephys mechanisms & parameters (using bpop.ephys module).

@author Lucas Koelman

@date   7/10/2017

"""

import bluepyopt.ephys as ephys

import logging
logger = logging.getLogger('bpop_ext')

import collections
from bgcellmodels.common import analysis


PROTOCOL_WRAPPERS = {} # StimProtocol : enum -> BpopProtocolWrapper


class BpopProtocolWrapper(object):
	"""
	Common attributes:

		ephys_protocol:		PhysioProtocol instance
		
		proto_vars:			dict with protocol variables
		
		response_interval:	expected time interval of response
	"""

	# SETPARAM: spike threshold for cell or specific protocol
	spike_threshold = -10.0
	
	@classmethod
	def make(cls, stim_proto, *args, **kwargs):
		"""
		Instantiate procol by name/enum value.
		"""
		if len(args) > 0:
			kwargs['stn_model_type'] = args[0]

		wrapper_class = PROTOCOL_WRAPPERS[stim_proto]
		return wrapper_class(**kwargs)


	def get_mechs_params(self):
		"""
		Get all ephys.mechanisms and ephys.parameters used by the protocols.

		These need to be assigned to the cell model to run the protocol.

		@return		tuple(mechs, params) containing a list of ephys.mechanisms 
					and ephys.parameters respectively
		"""

		proto_mechs = self.proto_vars.get('pp_mechs', []) + \
		              self.proto_vars.get('range_mechs', [])

		proto_params = self.proto_vars.get('pp_mech_params', [])

		return proto_mechs, proto_params


	@classmethod
	def all_mechs_params(cls, proto_wrappers):
		"""
		Concatenate all mechanisms and all parametes for given protocol wrappers.

		This is useful for assigning to a cell model that must be able to run multiple protocols.
		"""
		all_mechs, all_params = [], []
		for proto in proto_wrappers:
			mechs, params = proto.get_mechs_params()
			all_mechs.extend(mechs)
			all_params.extend(params)

		return all_mechs, all_params


# Helper functions for SelfContainedProtocol

def rng_getter(setup_kwargs):
	"""
	Function to get Numpy.Random object for stimulation protocol setup functions.
	"""
	import numpy as np
	base_seed = setup_kwargs['base_seed']
	return np.random.RandomState(base_seed)


def connector_getter(setup_kwargs):
	"""
	Function to get CellConnector for stimulation protocol setup functions.
	"""
	import bgcellmodels.cellpopdata as cpd
	physio_state = setup_kwargs['physio_state']
	rng = setup_kwargs['rng']
	return cpd.CellConnector(physio_state, rng)


class PhysioProtocol(ephys.protocols.SweepProtocol):
    """
    Wrapper for SweepProtocol that adds an initialization function.
    """

    def __init__(
            self,
            name=None,
            stimuli=None,
            recordings=None,
            cvode_active=None,
            init_func=None):
        """
        Constructor
        
        Args:
            init_func:  function(sim, model) that takes Simulator and instantiated 
                        CellModel (icell) as arguments in that order
        """

        self._init_func = init_func

        super(PhysioProtocol, self).__init__(
            name,
            stimuli=stimuli,
            recordings=recordings,
            cvode_active=cvode_active)


    def instantiate(self, sim=None, icell=None):
        """
        Instantiate

        NOTE: called in self._run_func() after model.instantiate()
        """
        # First apply physiological conditions
        self._init_func(sim, icell)

        # Then instantiate stimuli and recordings
        super(PhysioProtocol, self).instantiate(
            sim=sim,
            icell=icell)

    def destroy(self, sim=None):
        """
        Destroy protocol
        """

        # Make sure stimuli are not active in next protocol if cell model reused
        # NOTE: should better be done in Stimulus objects themselves for encapsulation, but BluePyOpt built-in Stimuli don't do this
        for stim in self.stimuli:
            if hasattr(stim, 'iclamp'):
                stim.iclamp.amp = 0
                stim.iclamp.dur = 0
            elif hasattr(stim, 'seclamp'):
                for i in range(3):
                    setattr(stim.seclamp, 'amp%d' % (i+1), 0)
                    setattr(stim.seclamp, 'dur%d' % (i+1), 0)

        # Calls destroy() on each stimulus
        super(PhysioProtocol, self).destroy(sim=sim)


class SelfContainedProtocol(ephys.protocols.SweepProtocol):
	"""
	Extension of BluePyOpt protocol that allow instantiating stimuli, synapses, etc. 
	without specifying them in the form of ephys mechanisms & parameters.
	"""

	def __init__(
			self,
			name=None,
			stimuli=None,
			recordings=None,
			cvode_active=None,
			total_duration=None,
			proto_init_funcs=None,
			proto_setup_funcs_pre=None,
			proto_setup_funcs_post=None,
			proto_setup_kwargs_const=None,
			proto_setup_kwargs_getters=None,
			rec_traces_funcs=None,
			plot_traces_funcs=None):
		"""
		Constructor
		
		@param init_func	function(sim, model) that takes Simulator and instantiated 
							CellModel (icell) as arguments in that order

		@param proto_setup_funcs_pre	protocol setup functions that need to be called
										before cell model is instantiated

		@param proto_setup_funcs_post	protocol setup functions that need to be called
										after cell model is instantiated

		@param proto_setup_funcs_kwargs_const	Constant keyword arguments to protocol
												setup functions

		@param proto_setup_funcs_kwargs_getters	Getter functions for keyword arguments to 
												protocol setup functions. These functions
												may create new entries in the passed dict.
		"""

		# init_physio(sim, icell)
		self._proto_init_funcs = [] if proto_init_funcs is None else list(proto_init_funcs)
		# rec_traces(icell, stim_data_dict, trace_spec_data, recorded_hoc_objects)
		self._funcs_rec_traces = [] if rec_traces_funcs is None else list(rec_traces_funcs)
		# plot_traces(trace_rec_data)
		self._funcs_plot_traces = [] if plot_traces_funcs is None else list(plot_traces_funcs)

		# Protocol setup functions applied to full model
		self._proto_setup_funcs_pre = [] if proto_setup_funcs_pre is None else list(proto_setup_funcs_pre)
		self._proto_setup_funcs_post = [] if proto_setup_funcs_post is None else list(proto_setup_funcs_post)

		# Common keyword arguments for all functions
		self._proto_setup_kwargs_const = proto_setup_kwargs_const if proto_setup_kwargs_const is not None else {}
		self._proto_setup_kwargs_getters = proto_setup_kwargs_getters if proto_setup_kwargs_getters is not None else {}

		# Data used by protocol setup functions
		self.recorded_trace_vectors = {}
		self.record_contained_traces = False
		self.autoplot_contained_traces = False

		########################################################################
		# SweepProtocol parameters

		if stimuli is None:
			stimuli = []

		self._total_duration = total_duration

		super(SelfContainedProtocol, self).__init__(
			name,
			stimuli=stimuli,
			recordings=recordings,
			cvode_active=cvode_active)


	@property
	def total_duration(self):
		"""
		Total duration of protocol (hides SweepProtocol.total_duration)
		"""
		return self._total_duration


	def run(self, cell_model, param_values, sim=None, isolate=None):
		"""
		Wrapper for _run_func() for execution over multiple threads.

		@post	if there are self-contained traces, they will available on
				attribute self.recorded_trace_vectors
		"""

		if isolate is None:
			isolate = True

		if isolate:
			def _reduce_method(meth):
				"""Overwrite reduce"""
				return (getattr, (meth.__self__, meth.__func__.__name__))

			import copyreg
			import types
			copyreg.pickle(types.MethodType, _reduce_method)

			import multiprocessing

			pool = multiprocessing.Pool(1, maxtasksperchild=1)

			# NOTE: the only difference with SweepProtocol is the second return
			# value, which is needed to marshal data between processes
			responses, traces = pool.apply(
				self._run_func,
				kwds={
					'cell_model': cell_model,
					'param_values': param_values,
					'sim': sim})

			pool.terminate()
			pool.join()
			del pool
		else:
			responses, traces = self._run_func(
				cell_model=cell_model,
				param_values=param_values,
				sim=sim)

		# 'responses' are protocol responses used by BluePyOpt, 'traces' are self-contained traces not used in optimization process
		self.recorded_trace_vectors = traces
		return responses


	def _run_func(self, cell_model, param_values, sim=None):
		"""
		Run protocols.

		Overrides SweepProtocol._run_func()

		NOTE: call graph for SweepProtocol is as follows:

			protocol.run(cell_model, param_values, sim) -> _run_func(...) :
				
				model.instantiate()
				
				protocol.instantiate() : <=== THIS FUNCTION
					stimulus.instantiate() :
						location.instantiate()
					recording.instantiate()
				
				sim.run()
		"""

		# try:

		# Fixes each param.value to individual's 'genes'
		cell_model.freeze(param_values)

		# Pass functions and parameters to cell_model before instantiation
		self.pre_model_instantiate(cell_model=cell_model, sim=sim)
		
		# Make final cell model
		cell_model.instantiate(sim=sim)
		
		# Instatiate things that need final cell model (original instantiate())
		self.post_model_instantiate(cell_model=cell_model, sim=sim)

		try:
			sim.run(self.total_duration, cvode_active=self.cvode_active)
		
		except (RuntimeError, ephys.simulators.NrnSimulatorException):
			
			logger.debug(
				'SelfContainedProtocol: Running of parameter set {%s} generated '
				'an exception, returning None in responses',
				str(param_values))
			
			responses = {recording.name: None 
							for recording in self.recordings}
		else:
			responses = {recording.name: recording.response
							for recording in self.recordings}

		# Cleanup functions before model & protocol destruction
		self.post_run(cell_model=cell_model, sim=sim)

		# NOTE: all protocol and model data is released here!
		self.destroy(sim=sim)
		cell_model.destroy(sim=sim)

		cell_model.unfreeze(param_values.keys())

		return responses, self.recorded_trace_vectors

		# except:
		# 	import sys
		# 	import traceback
		# 	raise Exception(
		# 		"".join(traceback.format_exception(*sys.exc_info())))


	def pre_model_instantiate(self, cell_model=None, sim=None):
		"""
		Function executed before cell model instantiation.
		"""
		# Create keyword arguments for protocol setup functions
		kwargs_default = {
			# 'icell': icell, # icell filled in by cell_model
			'nrnsim': sim.neuron.h,
			'stim_data': {}, # synapses and netcons
			'rng_info': {'stream_indices': {}} # map from low index to current highest index
		}

		# NOTE: proto kwargs must only be instantiated once per model instantiation
		logger.debug("Instantiating protocol setup kwargs...")
		self._iproto_data = kwargs_default
		self._iproto_data.update(self._proto_setup_kwargs_const)

		# kwargs getters add keyword arguments on-the-fly
		for kwarg_name, kwarg_getter in self._proto_setup_kwargs_getters.items():
			self._iproto_data[kwarg_name] = kwarg_getter(self._iproto_data)

		cell_model.proto_setup_funcs = self._proto_setup_funcs_pre
		cell_model.proto_setup_kwargs = self._iproto_data


	def post_model_instantiate(self, cell_model=None, sim=None):
		"""
		Function executed after cell model instantiation.
		"""
		# Make stimuli (inputs)
		for proto_setup_post in self._proto_setup_funcs_post:
			proto_setup_post(**self._iproto_data)

		# Make custom recordings and store
		if self.record_contained_traces:
			# Prepare recording data containers
			iproto_rec_data = {
				'icell': cell_model.icell,
				'trace_specs': collections.OrderedDict(),
				'trace_vectors': {},
				'rec_hoc_objects': {},
				'rec_hoc_markers': {},
			}
			self._iproto_data.update(iproto_rec_data)

			# Custom recording functions (put trace specs in trace_spec_data)
			for rec_func in self._funcs_rec_traces:
				rec_func(**self._iproto_data)
				logger.debug("Executed recording function {}".format(rec_func))

			# Create recording vectors
			if not 't_global' in self._iproto_data['trace_specs']:
				self._iproto_data['trace_specs']['t_global'] = {'var':'t'}
			
			recData, markers = analysis.recordTraces(
									self._iproto_data['rec_hoc_objects'], 
									self._iproto_data['trace_specs'],
									recordStep=0.05)

			self._iproto_data['trace_vectors'] = recData
			self._iproto_data['rec_hoc_markers'] = markers

		# Initialize physiological conditions
		for init_func in self._proto_init_funcs:
			init_func(**self._iproto_data)

		# Defined in superclass, instantiates ephys.stimuli and ephys.recordings
		self.instantiate(sim=sim, icell=cell_model.icell)


	def post_run(self, cell_model=None, sim=None):
		"""
		Function executed immediately after running simulation, when model and
		protocol data have not been destroyed.

		@effect		saves recorded Hoc.Vector objects in dictionary 
					self.recorded_trace_vectors as numpy arrays, with the original
					trace names used for recording as keys.
		"""
		if self.record_contained_traces:
			specs = self._iproto_data['trace_specs']
			vecs = self._iproto_data['trace_vectors']
			self.recorded_trace_vectors = collections.OrderedDict((
				(k, vecs[k].as_numpy()) for k in specs.keys() # preserve order
			))

		if self.autoplot_contained_traces:
			self.plot_contained_traces()


	def plot_contained_traces(self):
		"""
		Plot traces recorded as self-contained traces, i.e. traces not passed as
		ephys.recordings objects in constructor.
		"""
		for func in self._funcs_plot_traces:
			func(self.recorded_trace_vectors)


	def get_contained_traces(self):
		"""
		Return contained traces recorded using recording functions.

		@return		collections.OrderedDict<str,numpy.array> : traces dict
		"""
		return self.recorded_trace_vectors


	def destroy(self, sim=None):
		"""
		Destroy protocol
		"""
		logger.debug("Destroying SelfContainedProtocol")
		
		# Release data belonging to instantiated protocol
		self._iproto_data = None

		# Make sure stimuli are not active in next protocol if cell model reused
		for stim in self.stimuli:
			if hasattr(stim, 'iclamp'):
				stim.iclamp.amp = 0
				stim.iclamp.dur = 0
			elif hasattr(stim, 'seclamp'):
				for i in range(3):
					setattr(stim.seclamp, 'amp%d' % (i+1), 0)
					setattr(stim.seclamp, 'dur%d' % (i+1), 0)

		# Calls destroy() on each stimulus
		super(SelfContainedProtocol, self).destroy(sim=sim)