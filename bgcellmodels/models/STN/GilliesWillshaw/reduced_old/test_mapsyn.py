"""
Test suite for mapsyn.py (algorithms for mapping synaptic input locations).
"""

import collections
import logging
logging.basicConfig(format='%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logger = logging.getLogger(__name__) # create logger for this module

import numpy as np

import neuron
from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

from bgcellmodels.models.STN.GilliesWillshaw import gillies_model as gillies
import marasco_foldbased as marasco
import mapsyn
from bgcellmodels.common import analysis


def test_map_synapses(export_locals=False):
	"""
	Test for function get_syn_info()

	@param export_locals	bring function local variables into global namespace
							for interactive testing
	"""
	# Make STN cell
	soma, dends, stims = gillies.stn_cell_gillies()
	somaref, dendL_refs, dendR_refs = gillies.get_stn_refs()
	allsecrefs = [somaref] + dendL_refs + dendR_refs

	# Connect some synapses
	import numpy as np
	rng = np.random.RandomState(15203008)
	test_syns = []
	test_nc = []

	# Connect a few synapses per dendritic section
	n_syn_sec = 3
	for dend_ref in (dendL_refs + dendR_refs):
		dend_sec = dend_ref.sec
		# make synapses
		for i in range(n_syn_sec):
			syn_loc = rng.rand()
			syn = h.ExpSyn(dend_sec(syn_loc))
			syn.tau = rng.rand() * 20.
			syn.e = rng.rand() * 10.
			test_syns.append(syn)

			# make NetCon
			nc = h.NetCon(None, syn)
			nc.threshold = 0.
			nc.delay = 5.
			nc.weight[0] = rng.rand() * 10.
			test_nc.append(nc)

	# Make function to run simulation and bring cell to desired state
	def stn_setstate():
		# h.celsius = 35
		# h.set_aCSF(4)
		# h.v_init = -60
		# h.tstop = 460
		# h.dt = 0.025
		# h.init()
		# h.run()
		pass

	# Frequency for computing impedances
	Z_freq = 25.

	# Get synapse info
	save_ref_attrs = ['table_index', 'tree_index', 'gid']
	syn_info = mapsyn.get_syn_info(soma, allsecrefs, Z_freq=Z_freq, 
						init_cell=stn_setstate, save_ref_attrs=save_ref_attrs)

	# Create reduced cell
	eq_refs, newsecrefs = marasco.reduce_gillies_incremental(n_passes=7, 
									zips_per_pass=100)
	logger.debug("STN cell reduction complete.")

	# Map synapses to reduced cell
	rootref = next((ref for ref in newsecrefs if 'soma' in ref.sec.name()))
	mapsyn.map_synapses(rootref, newsecrefs, syn_info, stn_setstate, Z_freq)

	if export_locals:
		globals().update(locals())


def map_run_single_input(reduced=True):
	"""
	Test single input on full model, mapping to reduced model,
	and compare responses.
	"""

	# Make STN cell
	soma, dends, stims = gillies.stn_cell_gillies()
	somaref, dendL_refs, dendR_refs = gillies.get_stn_refs()
	dend_refs = dendL_refs + dendR_refs
	allsecrefs = [somaref] + dendL_refs + dendR_refs

	# Create one synapse
	i_tree = 0 # see diagram in marasco_reduction.pptx / slide "Indices"
	i_sec = 18  # FIXME: enter cell number here
	loc_syn = 0.5
	target_ref = [dendL_refs, dendR_refs][i_tree][i_sec-1]
	assert (target_ref.tree_index==i_tree and target_ref.table_index==i_sec)

	# asyn = h.AlphaSynapse(target_ref.sec(loc_syn))
	# asyn.onset = 325.0 # FIXME: enter spike time, (AlphaSynapse has no NET_RECEIVE)
	# asyn.tau = 5.0
	# asyn.gmax = 0.0005
	# asyn.e = 0.0
	asyn = h.Exp2Syn(target_ref.sec(loc_syn))
	asyn.tau1 = 5.0
	asyn.tau2 = 5.0
	asyn.e = 0.0
	gmax = 0.0015 # set using NetCon weight

	####################################
	# Simulate full model

	nc = h.NetCon(None, asyn)
	nc.delay = 1.0
	nc.weight[0] = gmax # weight is max conductance
	def queue_events():
		nc.event(325.0) # FIXME: enter spike time
	fih = h.FInitializeHandler(queue_events)

	# Run simulation
	secs = {'soma': soma, 'asyn': asyn}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)

	# Membrane voltages
	traceSpecs['V_soma'] = {'sec':'soma', 'loc':0.5, 'var':'v'}
	traceSpecs['i_syn'] = {'pointp':'asyn', 'var':'i'}

	# Set up recording vectors
	recordStep = 0.05
	recData, markers = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Change conductances
	def block_Na_channels():
		for sec in h.allsec():
			for seg in sec:
				seg.gna_NaL = 0.0 # prevent spontaneousdepolarization
				seg.gna_Na = 0.0 # prenvent spiking

	block_Na_channels()

	# Simulate
	h.dt = 0.025
	h.celsius = 35
	h.v_init = -60
	h.set_aCSF(4)
	h.tstop = 500
	h.init()
	h.run()

	# Plot
	figs_vm = analysis.plotTraces(recData, recordStep, traceSharex=True,
									yRange={'V_soma':(-80,0)})

	####################################
	# Reduce cell

	# Restore conductances
	gillies.reset_channel_gbar()

	def stn_setstate():
		""" Initialize cell to measure electrical properties """
		h.celsius = 35
		h.v_init = -68.0
		h.set_aCSF(4)
		h.init()

	# Get synapse info
	save_ref_attrs = ['table_index', 'tree_index', 'gid'] # properties to save
	Z_freq = 25.
	syn_info = mapsyn.get_syn_info(soma, allsecrefs, Z_freq=Z_freq, 
						init_cell=stn_setstate, save_ref_attrs=save_ref_attrs)

	# Create reduced cell
	marasco.logger.setLevel(logging.WARNING) # ignore lower than warning
	eq_refs, newsecrefs = marasco.reduce_gillies_incremental(n_passes=7, 
									zips_per_pass=100)

	# Map synapses to reduced cell
	rootref = next((ref for ref in newsecrefs if 'soma' in ref.sec.name()))
	mapsyn.map_synapses(rootref, newsecrefs, syn_info, stn_setstate, Z_freq,
						method='Vratio')

	####################################
	# Simulate reduced cell
	recData['i_syn'] = h.Vector() # prevents recording

	# Run simulation
	msyn = syn_info[0].mapped_syn
	secs = {'soma': soma, 'asyn': msyn}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)

	# Membrane voltages
	traceSpecs['V_soma'] = {'sec':'soma', 'loc':0.5, 'var':'v'}
	traceSpecs['i_syn'] = {'pointp':'asyn', 'var':'i'}

	# Set up recording vectors
	recordStep = 0.05
	recDataB, markers = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Simulate
	block_Na_channels()
	h.dt = 0.025
	h.celsius = 35
	h.v_init = -60
	h.set_aCSF(4)
	h.tstop = 500
	h.init()
	h.run()

	# Plot
	figs_B = analysis.plotTraces(recDataB, recordStep, traceSharex=True,
									yRange={'V_soma':(-80,0)})
	globals().update(locals())


def map_run_multiple_inputs():
	"""
	Test multiple inputs on full model, mapping to reduced model,
	and compare responses.
	"""

	# Make STN cell
	soma, dends, stims = gillies.stn_cell_gillies()
	somaref, dendL_refs, dendR_refs = gillies.get_stn_refs()
	dend_refs = dendL_refs + dendR_refs

	# Connect a few synapses per dendritic section
	rng = np.random.RandomState(15203008)
	n_syns = 10
	n_syns_sec = 1
	i_target_secs = rng.choice(len(dend_refs), n_syns, replace=False)

	test_syns = []
	test_nc = []
	n_syns_placed = 0

	for i_sec in i_target_secs:
		if n_syns_placed >= n_syns:
			break
		dend_sec = dend_refs[i_sec].sec

		# make N synapses per Section
		for i in range(n_syns_sec):

			# Make synapse with random properties
			syn_loc = rng.rand()
			asyn = h.AlphaSynapse(target_ref.sec(loc_syn))
			asyn.onset = rng.rand() * 5.0
			asyn.tau = rng.rand() * 10.0
			asyn.gmax = rng.rand() * 0.5
			asyn.e = [0., -85.][rng.randint(0,2)]

			n_syns_placed += 1
			test_syns.append(asyn)

			# make NetCon
			nc = h.NetCon(None, syn)
			nc.threshold = 0.
			nc.delay = 5.
			nc.weight[0] = rng.rand() * 10.
			test_nc.append(nc)

	raise NotImplementedError('TODO: finish this method')


if __name__ == '__main__':
	# test_map_synapses(export_locals=True)
	map_run_single_input()