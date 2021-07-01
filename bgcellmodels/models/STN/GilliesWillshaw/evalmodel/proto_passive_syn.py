"""
Test passive propagation of single EPSP & IPSP

@author	Lucas Koelman

@date	31/08/2017
"""

import neuron
h = neuron.h

import re

# Physiological parameters
import bgcellmodels.cellpopdata as cpd
from bgcellmodels.cellpopdata import (
	PhysioState,
	Populations as Pop,
	NTReceptors as NTR,
	ParameterSource as Cit
)

# Stimulation protocols
from proto_common import (
	StimProtocol, EvaluationStep, 
	register_step, pick_random_segments, extend_dictitem,
	logger
)

# Plotting
from bgcellmodels.common import analysis

################################################################################
# Interface functions
################################################################################

# Implemented stimulation protocol
impl_proto = StimProtocol.PASSIVE_SYN

@register_step(EvaluationStep.INIT_SIMULATION, impl_proto)
def init_sim(self, protocol):
	"""
	Initialize simulator to simulate background protocol
	"""
	model = self.target_model

	# Change simulation duration
	self._init_sim(dur=5000)

	# Make cell passive: disable all active channels
	gbar_active = self.model_data[model]['active_gbar_names']

	# for sec in h.allsec():
	for ref in self.all_sec_refs(model):
		for seg in ref.sec:
			for gbar in gbar_active:
				setattr(seg, gbar, 0.0)


@register_step(EvaluationStep.MAKE_INPUTS, impl_proto)
def make_inputs(self, connector=None):
	"""
	Make a number of either GLU (AMPA) or GABA (GABA-A) synapses distributed randomly
	over the dendritic tree, and firing sequentially with a fixed ISI.
	"""

	test_GLU = True # TODO: set this in user params on Evaluator object
	num_syn = 30

	if test_GLU:
		make_GLU_inputs(self, num_syn, connector=connector)
	else:
		make_GABA_inputs(self, num_syn, connector=connector)


@register_step(EvaluationStep.RECORD_TRACES, impl_proto)
def rec_traces(self, protocol, traceSpecs):
	"""
	Record all traces for this protocol.
	"""
	# record synaptic traces
	rec_EPSP(self, protocol, traceSpecs)

	# Record Vm in all recorded segments
	self.rec_Vm(protocol, traceSpecs)

	# Record input spikes
	rec_spikes(self, protocol, traceSpecs)


@register_step(EvaluationStep.PLOT_TRACES, impl_proto)
def plot_traces(self, model, protocol):
	"""
	Plot all traces for this protocol
	"""

	# Plot rastergrams
	spk_filter = lambda trace: 'AP' in trace
	self._plot_all_spikes(model, protocol, trace_filter=spk_filter, color='r')

	# Plot Vm in recorded segments
	figs = self._plot_all_Vm(model, protocol, fig_per='cell')
	for fig in figs:
		fig.suptitle('Vm (cell model {})'.format(model))

	# Get data
	rec_dict = self.model_data[model]['rec_data'][protocol]
	recData = rec_dict['trace_data']
	rec_dt = rec_dict['rec_dt']

	# Plot soma Vm separately
	traces = analysis.match_traces(recData, lambda t: 'soma' in t)
	figs = analysis.plotTraces(traces, rec_dt, oneFigPer='cell')


################################################################################
# Building block functions
################################################################################

def make_GLU_inputs(self, n_ctx_syn, connector=None):
	"""
	Make a single Glutamergic synapse on the STN neuron.
	"""
	if connector is None:
		cc = cpd.CellConnector(self.physio_state, self.rng)
	else:
		cc = connector

	model = self.target_model
		
	# Add CTX inputs using Tsodyks-Markram synapses
	# Distribute synapses over dendritic trees
	is_ctx_target = lambda seg: seg.diam <= 1.0         
	dend_secrefs = dend_secrefs = self.model_data[model]['dend_refs']
	ctx_target_segs = pick_random_segments(dend_secrefs, n_ctx_syn, is_ctx_target, rng=self.rng)

	# Make synapses
	new_inputs = {}
	for i_syn, target_seg in enumerate(ctx_target_segs):

		# Configure source NetStim
		t_transient = 350.
		t_between = 100.
		t_fire = t_transient + i_syn*t_between

		# Make spike generator
		stimsource = h.VecStim()
		spike_times = h.Vector(1,t_fire)
		stimsource.play(spike_times)

		# Synapse parameters (where to get them)
		syn_mech = 'Exp2Syn'
		syn_pre_post_pops = (Pop.CTX, Pop.STN)
		syn_source_target = (stimsource, target_seg)
		syn_NTR = (NTR.AMPA,)
		
		# Make synapse and NetCon
		syn, nc, wvecs = cc.make_synapse(syn_pre_post_pops, syn_source_target, syn_mech, 
							syn_NTR, use_sources=(Cit.Custom, Cit.Default))

		# Debug
		# nc.weight[0] = 10
		# nc.delay = 0
		logger.anal("Made {} synapse with t={}ms, PRE_POP={} , NTR={}".format(
				syn_mech, t_fire, syn_pre_post_pops[0], syn_NTR))

		# Save inputs
		extend_dictitem(new_inputs, 'syn_NetCons', nc)
		extend_dictitem(new_inputs, 'synapses', syn)
		extend_dictitem(new_inputs, 'NetStims', stimsource)
		extend_dictitem(new_inputs, 'stimtimevec', spike_times)

	self.add_inputs('ctx', model, **new_inputs)


def make_GABA_inputs(self, n_gpe_syn, connector=None):
	"""
	Make GABAergic synapses on the STN neuron.

	@param self         StnModelEvaluator object
	"""
	
	if connector is None:
		cc = cpd.CellConnector(self.physio_state, self.rng)
	else:
		cc = connector

	model = self.target_model

	# Pick random segments in dendrites for placing synapses
	is_gpe_target = lambda seg: seg.diam > 1.0 # select proximal dendrites
	dend_secrefs = self.model_data[model]['dend_refs']
	gpe_target_segs = pick_random_segments(dend_secrefs, n_gpe_syn, is_gpe_target, rng=self.rng)

	# Make synapses
	new_inputs = {}
	for i_syn, target_seg in enumerate(gpe_target_segs):

		# Synapse parameters (where to get them)
		syn_mech = 'Exp2Syn'
		syn_pre_post_pops = (Pop.GPE, Pop.STN)
		syn_source_target = (None, target_seg)
		syn_NTR = (NTR.GABAA,)
		
		# Make synapse and NetCon
		syn, nc, wvecs = cc.make_synapse(syn_pre_post_pops, syn_source_target, syn_mech, 
							syn_NTR, use_sources=(Cit.Custom, Cit.Default))

		# Configure source NetStim
		t_transient = 350.
		t_between = 15.
		t_fire = t_transient + i_syn*t_between
		
		# Queue events
		def queue_events():
			nc.event(t_fire)
		fih = h.FInitializeHandler(queue_events)

		# Save inputs
		extend_dictitem(new_inputs, 'PyInitHandlers', queue_events)
		extend_dictitem(new_inputs, 'HocInitHandlers', fih)
		extend_dictitem(new_inputs, 'syn_NetCons', nc)
		extend_dictitem(new_inputs, 'synapses', syn)

	self.add_inputs('gpe', model, **new_inputs)


def rec_EPSP(self, protocol, traceSpecs):
	"""
	Record traces at GLU synapses

	NOTE: all suffixes containing an index refer to the index of the
	      synapse object => this makes them identifiable in
	      NEURON GUI Point Process Manager
	"""
	model = self.target_model
	inputs = self.model_data[model]['inputs']
	rec_segs = self.model_data[self.target_model]['rec_segs'][protocol]

	n_syn = 15 # number of recorded synapses
	
	# Get NetCons
	nc_list = []
	for pop in ('ctx', 'gpe'):
		if pop in inputs:
			nc_list.extend(inputs[pop]['syn_NetCons'])
	
	# Record EPSP in post-synaptic segment
	for i_syn, nc in enumerate(nc_list):
		if i_syn > n_syn-1:
			break

		match = re.search(r'\[[\d]+\]', nc.syn().hname())
		suffix = match.group() if match else str(i_syn)

		# Add post-synaptic segment to recorded segments
		seg_label = 'sync' + suffix
		rec_segs[seg_label] = nc.syn().get_segment()

		# Record Vm in post-synaptic segment -> use evaluator.rec_Vm
		# traceSpecs['V_'+seg_label] = {'sec':seg_label, 'loc':seg.x, 'var':'v'}


def rec_spikes(self, protocol, traceSpecs):
	"""
	Record input spikes delivered to synapses.

	NOTE: all suffixes containing an index refer to the index of the
	      synapse object => this makes them identifiable in
	      NEURON GUI Point Process Manager
	"""
	model = self.target_model
	rec_pops = [Pop.GPE, Pop.CTX]

	for pre_pop in rec_pops:
		
		inputs = self.model_data[model]['inputs'].get(pre_pop.name.lower(), None)
		if inputs is None:
			continue
		
		nc_list = inputs['syn_NetCons']
		for i_syn, nc in enumerate(nc_list):

			# Add NetCon to list of recorded objects
			match = re.search(r'\[[\d]+\]', nc.syn().hname())
			suffix = match.group() if match else ('syn' + str(i_syn))
			syn_tag = pre_pop.name + suffix
			self._add_recorded_obj(syn_tag, nc, protocol)

			# Specify trace
			traceSpecs['AP_'+syn_tag] = {'netcon':syn_tag}