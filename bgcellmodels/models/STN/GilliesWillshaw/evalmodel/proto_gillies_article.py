"""
Stimulation protocols to reproduce figures in Gillies & Willshaw (2005) paper,
taken from original code released with article.

"""

import collections

import neuron
h = neuron.h

from bgcellmodels.common import analysis

import gillies_model as gillies
from bgcellmodels.cellpopdata import StnModel
from proto_common import StimProtocol, EvaluationStep, register_step

################################################################################
# Protocol CLAMP_REBOUND
################################################################################


@register_step(EvaluationStep.INIT_SIMULATION, StimProtocol.CLAMP_REBOUND)
def init_sim_rebound(self, protocol):
	"""
	Initialize simulator to simulate background protocol
	"""
	model = self.target_model

	# Custom  version of _init_sim()
	dur = 2000
	dt = 0.025
	celsius = 35
	v_init = -60
	# self._init_sim(dur=2000, dt=0.025, celsius=25, v_init=-60)

	self.sim_dur = dur
	h.tstop = dur
	
	self.sim_dt = dt
	h.dt = dt

	h.celsius = celsius # different temp from paper
	h.v_init = v_init # paper simulations USE default v_init
	gillies.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

	# Initialize NEURON simulator
	h.init() # calls finitialize()


@register_step(EvaluationStep.MAKE_INPUTS, StimProtocol.CLAMP_REBOUND)
def make_inputs_rebound(self, connector=None):
	"""
	Make a number of either GLU (AMPA) or GABA (GABA-A) synapses distributed randomly
	over the dendritic tree, and firing sequentially with a fixed ISI.
	"""
	# Get data
	model = self.target_model
	soma_refs = self.model_data[model]['soma_refs']
	soma_sec = soma_refs[0].sec

	if model == StnModel.Gillies2005:
		stim1, stim2, stim3 = h.stim1, h.stim2, h.stim3
	else:
		stim1, stim2, stim3 = [h.IClamp(soma_sec(0.5)) for i in range(3)]

	# Set stimulation protocol (pulse amplitudes & durations)
	# Original code: stim2.amp = -0.25 (hyperpolarize to -75 mV steady state)
	stim1.delay = 0
	stim1.dur = 500
	stim1.amp = 0.0

	# stim2.delay = 200
	# stim2.dur = 500
	# stim2.amp = -0.11 # -0.25 in full model (hyperpolarize to -75 mV steady state)
	stim2.amp = 0

	# Use voltage clamp (space clamp) instead of current clamp
	clamp = h.SEClamp(soma_sec(0.5))
	clamp.dur1 = 0
	clamp.dur2 = 0
	clamp.dur3 = 500
	clamp.amp3 = -75

	stim3.delay = 1000
	stim3.dur = 1000
	stim3.amp = 0.0

	# Save inputs
	new_inputs = {
		'IClamps' :	[stim1, stim2, stim3],
		'SEClamps': [clamp],
	}
	self.add_inputs('electrodes', model, **new_inputs)


@register_step(EvaluationStep.RECORD_TRACES, StimProtocol.CLAMP_REBOUND)
def rec_traces_rebound(self, protocol, traceSpecs):
	"""
	Record all traces for this protocol.
	"""

	model = self.target_model
	rec_segs = self.model_data[model]['rec_segs'][protocol]

	# Trace specs for membrane voltages
	for seclabel, seg in rec_segs.iteritems():
		traceSpecs['V_'+seclabel] = {'sec':seclabel, 'loc':seg.x, 'var':'v'}

	# Trace specs for recording ionic currents, channel states
	analysis.rec_currents_activations(traceSpecs, 'soma', 0.5)
	
	# Ca and K currents in distal dendrites
	dendloc = rec_segs['dist_dend'].x
	analysis.rec_currents_activations(traceSpecs, 'dist_dend', dendloc, ion_species=['ca','k'])


@register_step(EvaluationStep.PLOT_TRACES, StimProtocol.CLAMP_REBOUND)
def plot_traces_rebound(self, model, protocol):
	"""
	Plot all traces for this protocol
	"""
	model = self.target_model
	recData = self.model_data[model]['rec_data'][protocol]['trace_data']
	recordStep = self.model_data[model]['rec_data'][protocol]['rec_dt']

	# Plot Vm in recorded segments
	self._plot_all_Vm(model, protocol)

	# Plot ionic currents, (in)activation variables
	figs_soma, cursors_soma = analysis.plot_currents_activations(recData, recordStep, 
											sec_tag='soma')
	
	figs_dend, cursors_dend = analysis.plot_currents_activations(recData, recordStep, 
											sec_tag='dist_dend')


################################################################################
# Protocol CLAMP_PLATEAU
################################################################################


@register_step(EvaluationStep.INIT_SIMULATION, StimProtocol.CLAMP_PLATEAU)
def init_sim_plateau(self, protocol):
	"""
	Initialize simulator to simulate background protocol
	"""

	# Test model changes
	# gleak_orig = 7.84112e-05
	# gleak_fit = 12.43169e-5 # fit to match Zin_DC (see praxis_passive.py)
	# cm_orig = 1.0
	# cm_fit = cm_orig * (gleak_fit / gleak_orig)
	# for sec in h.allsec():
	# 	if not ('soma' in sec.name()):
	# 		for seg in sec:
	# 			setattr(seg, gillies.gleak_name, gleak_fit)
	# 			setattr(seg, 'cm', cm_fit)

	# Custom  version of _init_sim()
	dur = 2000
	dt = 0.025
	celsius = 30 # different temp from paper
	v_init = -60
	# self._init_sim(dur=2000, dt=0.025, celsius=25, v_init=-60)

	self.sim_dur = dur
	h.tstop = dur
	
	self.sim_dt = dt
	h.dt = dt

	h.celsius = celsius # different temp from paper
	h.v_init = v_init # paper simulations sue default v_init
	gillies.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

	# Initialize NEURON simulator
	h.init() # calls finitialize()


@register_step(EvaluationStep.MAKE_INPUTS, StimProtocol.CLAMP_PLATEAU)
def make_inputs(self, connector=None):
	"""
	Make inputs for this protocol
	"""

	# Get data
	model = self.target_model
	proto = StimProtocol.CLAMP_PLATEAU

	if model == StnModel.Gillies2005:
		stim1, stim2, stim3 = h.stim1, h.stim2, h.stim3
	else:
		soma_refs = self.model_data[model]['soma_refs']
		soma_sec = soma_refs[0].sec
		stim1, stim2, stim3 = [h.IClamp(soma_sec(0.5)) for i in range(3)]

	dur = self.sim_dur

	# Set up stimulation (5 mA/cm2 for 80 ms)
	I_hyper = -0.17 # hyperpolarize to -70 mV (see fig. 10C)
	I_depol = I_hyper + 0.2 # see fig. 10D: 0.2 nA (=stim.amp) over hyperpolarizing current
	dur_depol = 50 # see fig. 10D, top right
	del_depol = 1000
	burst_time = [del_depol-50, del_depol+200] # empirical
	self.model_data[model]['proto_vars'][proto]['burst_time'] = burst_time

	stim1.delay = 0
	stim1.dur = del_depol
	stim1.amp = I_hyper

	stim2.delay = del_depol
	stim2.dur = dur_depol
	stim2.amp = I_depol

	stim3.delay = del_depol + dur_depol
	stim3.dur = dur - (del_depol + dur_depol)
	stim3.amp = I_hyper

	# Save inputs
	new_inputs = {
		'IClamps' :	[stim1, stim2, stim3],
	}
	self.add_inputs('electrodes', model, **new_inputs)


@register_step(EvaluationStep.RECORD_TRACES, StimProtocol.CLAMP_PLATEAU)
def rec_traces_plateau(self, protocol, traceSpecs):
	"""
	Record all traces for this protocol.
	"""
	# Get data
	model = self.target_model
	
	rec_segs = self.model_data[model]['rec_segs'][protocol]
	# rec_segs = {
	# 	'soma': h.SThcell[0].soma, # middle of soma
	# 	'dist_dend': h.SThcell[0].dend1[7] # approximate location along dendrite in fig. 5C
	# }

	# Trace specs for membrane voltages
	for seclabel, seg in rec_segs.iteritems():
		traceSpecs['V_'+seclabel] = {'sec':seclabel, 'loc':seg.x, 'var':'v'}

	# Record Ca and Ca-activated currents in dendrite
	dendloc = rec_segs['dist_dend'].x

	# K currents (dendrite)
	traceSpecs['I_KCa_d'] = {'sec':'dist_dend','loc':dendloc,'mech':'sKCa','var':'isKCa'}
	
	# Ca currents (dendrite)
	traceSpecs['I_CaL_d'] = {'sec':'dist_dend','loc':dendloc,'mech':'HVA','var':'iLCa'}
	traceSpecs['I_CaN_d'] = {'sec':'dist_dend','loc':dendloc,'mech':'HVA','var':'iNCa'}
	traceSpecs['I_CaT_d'] = {'sec':'dist_dend','loc':dendloc,'mech':'CaT','var':'iCaT'}


@register_step(EvaluationStep.PLOT_TRACES, StimProtocol.CLAMP_PLATEAU)
def plot_traces_plateau(self, model, protocol):
	"""
	Plot all traces for this protocol
	"""
	# Get data
	proto = StimProtocol.CLAMP_PLATEAU
	model = self.target_model
	recData = self.model_data[model]['rec_data'][protocol]['trace_data']
	recordStep = self.model_data[model]['rec_data'][protocol]['rec_dt']

	# Plot membrane voltages
	self._plot_all_Vm(model, protocol)

	# Plot ionic currents, (in)activation variables
	figs, cursors = analysis.plot_currents_activations(recData, recordStep)

	# Dendrite currents during burst
	recDend = collections.OrderedDict([(k,v) for k,v in recData.iteritems() if k.endswith('_d')])

	burst_time = self.model_data[model]['proto_vars'][proto]['burst_time']
	analysis.cumulPlotTraces(recDend, recordStep, timeRange=burst_time)


if __name__ == '__main__':
	print("""
	Run this protocol using stn_model_evaluation.py:

	cd bgcellmodels/GilliesWillshaw
	from evalmodel.stn_model_evaluation import *
	run_protocol(StimProtocol.CLAMP_PLATEAU, StnModel.Gillies2005)
	""")