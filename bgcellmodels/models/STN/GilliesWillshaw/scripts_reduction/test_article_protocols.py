"""
Run experiments using the Gillies & Willshaw (2006) STN neuron model

The experiments are designed to discover which currents are responsible
for the different features of STN cell dynamics

@author Lucas Koelman
@date	28-10-2016
@note	must be run from script directory or .hoc files not found

"""

import collections
import math
import sys, os.path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Make sure other modules are on Python path
scriptdir, scriptfile = os.path.split(__file__)
repo_root = os.path.normpath(os.path.join(scriptdir, '..'))
sys.path.append(repo_root)

# Our own modules
import gillies_model as gillies
from gillies_model import gillies_gdict, gillies_mechs, gillies_glist

# Load NEURON
from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

from bgcellmodels.common import analysis
from bgcellmodels.common.analysis import rec_currents_activations, plot_currents_activations

from reducemodel import (
	analyze_reduction,
	reduce_bush_sejnowski as bush,
	marasco_clusterbased as marasco_cluster,
	marasco_foldbased as marasco_fold,
	reduce_cell
)
from reducemodel.redutils import lambda_AC, ExtSecRef, getsecref
from reducemodel.interpolation import *


# Global variables
soma = None


def stn_cell_gillies(resurgent=False):
	""" Initialize full Gillies & Willshaw cell model with two dendritic trees
	@param resurgent	if True, resurgent Na current mechanism
						from Do & Bean (2001) is added to the cell
	"""
	# Create cell and three IClamp in soma
	global soma
	if soma is None:
		h.xopen("gillies_cell_singleton.hoc")
		soma = h.SThcell[0].soma
	dends = h.SThcell[0].dend0, h.SThcell[0].dend1
	stims = h.stim1, h.stim2, h.stim3

	# Add resurgent current if requested
	if resurgent:
		default_gNarsg_soma = 0.016 # same as in .mod file and Akeman papeer
		default_gNarsg_dend = 1.0e-7 # default value for Na in Gillies model
		soma.insert('Narsg')
		soma.gbar_Narsg = default_gNarsg_soma
		alldendsecs = list(dends[0]) + list(dends[1])
		for sec in alldendsecs:
			sec.insert('Narsg')
			sec.gbar_Narsg = default_gNarsg_dend

	return soma, dends, stims

def stn_cell_Rall(resurgent=False, oneseg=False):
	""" Initialise reduced Gilliew & Willshaw model with each dendritic tree
	reduced to a single equivalent cylinder/section according to Rall's
	reduction method and the number of segments/points determined by
	the electrotonic length.

	@param oneseg	if True, a single segment is usef for each section,
					otherwise the number of segments is determined by 
					dL < 0.1*lambda(100)

	OBSERVATIONS:
	- (nseg from lambda)
		- same behavior as original model except spontaneous firing

	- (nseg=1)
		- cannot get same behavour as in original model or reduced model
		  with large amount of compartments. Processes/interaction of membrane
		  mechanisms in distal dendrite is not sufficiently isolated from
		  processes in soma.
	"""
	# Properties from SThprotocell.hoc
	all_Ra = 150.224
	all_cm = 1.0
	soma_L = 18.8
	soma_diam = 18.3112

	# List of mechanisms/conductances
	stn_mechs = list(gillies_mechs) # List of mechanisms to insert
	stn_glist = list(gillies_glist) # List of channel conductances to modify
	if resurgent: # add Raman & bean (2001) resurgent channel model
		stn_mechs.append('Narsg')
		stn_glist.append('gbar_Narsg')

	# Create soma
	soma = h.Section()
	soma.nseg = 1
	soma.Ra = all_Ra
	soma.diam = soma_diam
	soma.L = soma_L
	soma.cm = all_cm
	for mech in stn_mechs: soma.insert(mech)
	setconductances(soma, -1, glist=stn_glist)
	gillies.setionstyles_gillies(soma)
	
	# Right dendritic tree (Fig. 1, small one)
	dend1 = h.Section()
	dend1.connect(soma(0))
	dend1.diam = 1.9480 # diam of topmost parent section
	dend1.L = 549.86 # equivalent length using Rall model
	dend1.Ra = all_Ra
	dend1.cm = all_cm
	# Configure ion channels (distribution,...)
	for mech in stn_mechs: dend1.insert(mech)
	if oneseg:
		dend1.nseg = 1
		setconductances(dend1, 1, fixbranch=8, fixloc=0.98)
	else:
		opt_nseg1 = int(np.ceil(dend1.L/(0.1*lambda_AC(dend1,100))))
		dend1.nseg = opt_nseg1
		setconductances(dend1, 1, glist=stn_glist)
	gillies.setionstyles_gillies(dend1) # set ion styles

	# Left dendritic tree (Fig. 1, big one)
	dend0 = h.Section()
	dend0.connect(soma(1))
	dend0.diam = 3.0973 # diam of topmost parent section
	dend0.L = 703.34 # equivalent length using Rall model
	dend0.Ra = all_Ra
	dend0.cm = all_cm
	# Configure ion channels (distribution,...)
	for mech in stn_mechs: dend0.insert(mech) # insert mechanisms
	if oneseg:
		dend0.nseg = 1
		setconductances(dend0, 0, fixbranch=10, fixloc=0.98)
	else:
		opt_nseg0 = int(np.ceil(dend0.L/(0.1*lambda_AC(dend0,100))))
		dend0.nseg = opt_nseg0
		setconductances(dend0, 0, glist=stn_glist)
	gillies.setionstyles_gillies(dend0) # set ion styles

	# Create stimulator objects
	stim1 = h.IClamp(soma(0.5))
	stim2 = h.IClamp(soma(0.5))
	stim3 = h.IClamp(soma(0.5))

	return soma, (dend0, dend1), (stim1, stim2, stim3)

def stn_cell_onedend():
	""" 
	Soma with a single equivalent dendritic tree reduced to two compartments: 
	one equivalent trunk/smooth sections + one equivalent spiny sections.

	The single tree is equivalent to three instances
	of the small (fig. 1, right) dendritic tree connected to the soma.
	The number og segments is determined empirically.

	[0..soma..1]-[0..smooth..1]-[0..spiny..1]

	TODO: use surface equivalence to scale conductances/cm and conserve axial resistance?
	"""
	# Passive properties for all segments (from SThprotocell.hoc)
	all_Ra = 150.224
	all_cm = 1.0

	# List of mechanisms/conductances
	stn_mechs = list(gillies_mechs) # List of mechanisms to insert
	stn_glist = list(gillies_glist) # List of channel conductances to modify

	# Create sections
	secnames = ['soma', 'trunk', 'spiny']
	for name in secnames:
		# ensures references are not destroyed and names are not memory addresses
		h("create %s" % name)

	# Create soma
	soma = getattr(h, secnames[0])
	soma.nseg = 1
	# Passive electrical ppties (except Rm/gpas)
	soma.L = 18.8
	soma.diam = 18.3112
	soma.Ra = all_Ra
	soma.cm = all_cm
	# Ion channels
	for mech in stn_mechs:
		soma.insert(mech)
	setconductances(soma, -1, glist=stn_glist)
	gillies.setionstyles_gillies(soma)
	
	# Primary neurite, equiv. to trunk/smooth dendrites
	trunk = getattr(h, secnames[1])
	trunk.connect(soma, 1, 0) # connect parent@0.0 to child@1.0
	# Passive electrical ppties (except Rm/gpas)
	trunk.L = 80 # sum of diameters branch 1 & 2
	trunk.diam = (1.9480+1.2272)/2.0 * 1.0 # 3 times average diameter of branch 1 & 2
	trunk.Ra = all_Ra
	trunk.cm = all_cm
	opt_nseg1 = int(np.ceil(trunk.L/(0.1*lambda_AC(trunk,100))))
	trunk.nseg = 9
	# Ion channels
	for mech in stn_mechs:
		trunk.insert(mech)
	interpconductances(trunk, 1, path=[1,2], glist=stn_glist) # set channel conductances
	gillies.setionstyles_gillies(trunk) # set ion styles

	# Secondary neurite, equivalent to spiny dendrites
	spiny = getattr(h, secnames[2])
	spiny.connect(trunk, 1, 0)
	# Passive electrical ppties (except Rm/gpas)
	spiny.L = 289 # length of long branch in full model
	spiny.diam = 0.7695 * 1.0 # 3 times diameter of long branch (branch 5)
	spiny.Ra = all_Ra
	spiny.cm = all_cm
	opt_nseg0 = int(np.ceil(spiny.L/(0.1*lambda_AC(spiny,100))))
	spiny.nseg = 9
	# Ion channels
	for mech in stn_mechs:
		spiny.insert(mech)
	interpconductances(spiny, 1, path=[5], glist=stn_glist) # set channel conductances
	gillies.setionstyles_gillies(spiny) # set ion styles

	for secname in secnames:
		sec = getattr(h, secname)
		print("Created section '{}'' with {} segments".format(sec.name(), sec.nseg))

	# Create stimulator objects
	stim1 = h.IClamp(soma(0.5))
	stim2 = h.IClamp(soma(0.5))
	stim3 = h.IClamp(soma(0.5))

	return soma, (trunk, spiny), (stim1, stim2, stim3)

def stn_cell(cellmodel):
	""" 
	Create stn cell with given identifier 

	ID		model description
	-------------------------
	1		Original Gillies & Willshaw STN cell model

	2		Rall reduction of Gillies & Willshaw STN model
			to one equivalent section per dendritic tree.
			Number of segments is determined using lambda rule.

	3		Bush & Sejnowski reduction of Gillies & Willshaw STN model

	4		Marasco reduction of Gillies & Willshaw STN model
			using custom clustering criterion based on diameter

	5		Marasco reduction of Gillies & Willshaw STN model using
			modified reduction algorithm where axial resistance is not
			averaged over cluster subtrees but compounded and 
			input resistance is conserved.

	6		Manual reduction of Gillies & Willshaw STN model with a
			single equivalent denritic tree containing one equivalent
			smooth and one spiny section next to the soma section.
			Parameters are fitted or empirically determined to produce
			qualitatively similar behaviour to the full model.

	8		Lucas' algorithm: incremental zipping of tree branches.
	"""
	stims = None

	if cellmodel==1: # Full model
		soma, dends, stims = stn_cell_gillies()
		# Load section indicated with arrow in fig. 5C
		# If you look at tree1-nom.dat it should be the seventh entry 
		# (highest L and nseg with no child sections of which there are two instances)
		dendsec = h.SThcell[0].dend1[7]
		dendloc = 0.8 # approximate location along dendrite in fig. 5C
		allsecs = [soma] + list(dends[0]) + list(dends[1])

	elif cellmodel==2: # Rall - optimal nseg
		soma, dends, stims = stn_cell_Rall()
		dendsec = dends[1]
		dendloc = 0.9
		allsecs = [soma] + list(dends)

	elif cellmodel==3: # Bush & Sejnowski method
		clusters, eq_secs = bush.reduce_bush_sejnowski()
		# filepath = "C:\\Users\\lkoelman\\Downloads\\stn_reduced_bush.p"
		# bush.save_clusters(clusters, filepath)
		# clusters = bush.load_clusters(filepath)
		# eq_secs = bush.rebuild_sections(clusters)
		soma = next(sec for sec in eq_secs if sec.name().startswith('soma'))
		dendsec = next(sec for sec in eq_secs if sec.name().startswith('spiny'))
		dendloc = 0.9
		allsecs = eq_secs

	elif cellmodel==6: # Optimized Bush model
		reduced_model_path = "C:\\Users\\lkoelman\\cloudstore_m\\simdata\\bush_sejnowski\\stn_reduced_bush.p"
		stn_controller = STNCellController(Protocol.SPONTANEOUS, reduced_model_path)
		# Fittest candidate for SPONTANEOUS optimization
		# stn_controller.gbar_adjust_allsec = True # whether gbar scaling wil apply to all sections
		# fittest_candidate = [
		# 	'soma_diam_factor': 0.563915223257359,
		# 	'soma_cm_factor': 1.4561344744499622,
		# 	'dend_cm_factor': 2.0,
		# 	'dend_Rm_factor': 2.0,
		# 	'dend_Ra': 228.0626727668688, # 150 in full model
		# 	'dend_diam_factor': 0.6617425835741709,
		# 	'gbar_sca_gna_NaL': 1.0,
		# ]
		# Fittest candidate for REBOUND optimization
		stn_controller.gbar_adjust_allsec = False # don't adjust gbar in soma
		fittest_candidate = {
			# soma properties
			'soma_diam_factor': 0.5166409889357284,
			'soma_cm_factor': 1.0,
			# dendrite passive properties
			'dend_cm_factor': 1.0046280491823951,
			'dend_Rm_factor': 1.9313327747357878,
			'dend_Ra': 193.88153564188488,
			'dend_diam_factor': 0.7810045138119329,
			# dendrite active properties
			'gbar_sca_gk_sKCa': 1.0052884459348173, # distal double step dist
			'gbar_sca_gcaL_HVA': 1.4765306846140291, # distal linear dist
			'gbar_sca_gcaT_CaT': 1.4500118248990195, # distal linear dist
		}
		soma, dendsecs = stn_controller.build_candidate(fittest_candidate)
		dendsec = next(sec for sec in dendsecs if sec.name().startswith('spiny'))
		dendloc = 0.9
		allsecs = [soma] + dendsecs

	elif cellmodel==4 or cellmodel==5: # Marasco - custom clustering

		marasco_method = True # whether trees will be averaged (True, as in paper) or Rin conserved (False)
		if cellmodel==5:
			marasco_method = Fals
		clusters, eq_secs, eq_refs = marasco_cluster.reduce_gillies_pathRi(
			customclustering=True, average_Ri=marasco_method)
		soma, dendLsecs, dendRsecs = eq_secs
		dendsec = dendRsecs[-1] # last/most distal section of small dendrite
		dendloc = 0.9
		allsecs = [soma] + dendLsecs + dendRsecs

	elif cellmodel==8: # Lucas ZipBranch algorithm

		# eq_refs, newsecrefs = marasco_fold.reduce_gillies_incremental(n_passes=7, zips_per_pass=100)

		soma_refs, dend_refs = reduce_cell.fold_gillies_marasco(False)
		newsecrefs = soma_refs + dend_refs
		
		soma = next(ref.sec for ref in newsecrefs if ref.sec.name().endswith('soma'))
		dendsec = next(ref.sec for ref in newsecrefs if (len(ref.sec.children())==0))
		dendloc = 0.9
		print("Distal segment for recording is {0}".format(dendsec(dendloc)))
		allsecs = [ref.sec for ref in newsecrefs]

	# Insert stimulation electrodes
	if stims is None:
		stim1 = h.IClamp(soma(0.5))
		stim2 = h.IClamp(soma(0.5))
		stim3 = h.IClamp(soma(0.5))
		stims = stim1, stim2, stim3
	
	return soma, [(dendsec, dendloc)], stims, allsecs

################################################################################
# Functions for reduced model
################################################################################



################################################################################
# Experiments
################################################################################

def test_spontaneous(soma, dends_locs, stims, resurgent=False):
	""" Run rest firing experiment from original Hoc file

	@param soma			soma section
	@param dends_locs	list of tuples (sec, loc) containing a section
						and x coordinate to place recording electrode
	@param stims		list of electrodes (IClamp)

	PROTOCOL
	- spontaneous firing: no stimulation

	OBSERVATIONS (GILLIES MODEL)

	- during ISI:
		- INaP (INaL) is significant (-0.02) in ISI and looks like main mechanism responsible
		  for spontaneous depolarization/firing
		- IKCa (repolarizing) is also significant (+0.004) during ISI
			- sKCa is responsible for most of the AHP, as it should
		- ICaT slowly activates in dendrites but is small (0.0008)
			- might help spontaneous depolarization

	- during AP
		- the HVA currents ICaL and ICaN contribute to depolarization
		- peak IKv3 is twice as high as IKDR (contributes more to repolarization)
	

	"""
	# Set simulation parameters
	dur = 2000
	h.dt = 0.025
	h.celsius = 35 # different temp from paper (fig 3B: 25degC, fig. 3C: 35degC)
	h.v_init = -60 # paper simulations use default v_init
	gillies.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

	# Recording: trace specification
	secs = {'soma': soma}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)

	# Membrane voltages
	traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
	traceSpecs['t_global'] = {'var':'t'}
	for i, (dend,loc) in enumerate(dends_locs):
		dendname = 'dend%i' % i
		secs[dendname] = dend
		traceSpecs['V_'+dendname] = {'sec':dendname,'loc':loc,'var':'v'}

	# Record ionic currents, open fractions, (in)activation variables
	rec_currents_activations(traceSpecs, 'soma', 0.5)

	# Set up recording vectors
	recordStep = 0.05
	recData, _ = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Simulate
	h.tstop = dur
	h.init() # calls finitialize() and fcurrent()
	t0 = h.startsw()
	h.run()
	t1 = h.startsw()
	h.stopsw() # or t1=h.startsw(); runtime = t1-t0
	print("Simulation ran for {:.6f} seconds".format(t1-t0))

	# Plot membrane voltages
	recV = collections.OrderedDict([(k,v) for k,v in recData.items() if k.startswith('V_')]) # preserves order
	figs_vm = analysis.plotTraces(recV, recordStep, yRange=(-80,40), traceSharex=True)
	vm_fig = figs_vm[0]
	vm_ax = figs_vm[0].axes[0]

	# Plot ionic currents, (in)activation variables
	figs, cursors = plot_currents_activations(recData, recordStep)

	# Save trace to file
	# V_soma = np.array(recData['V_soma'], ndmin=2)
	# T_soma = np.array(recData['t_global'], ndmin=2)
	# TV_soma = np.concatenate((T_soma, V_soma), axis=0) * 1e-3 # pyelectro expects SI units: seconds, Volts
	# fpath = 'C:\\Users\\lkoelman\\cloudstore_m\\simdata\\fullmodel\\cstim_fullmodel_Vm_dt25e-3_0ms_2000ms.csv'
	# np.savetxt(fpath, TV_soma.T, delimiter=',', fmt=['%.3E', '%.7E'])

	plt.show(block=False)
	return recData, figs, cursors

def test_plateau(soma, dends_locs, stims):
	"""
	Test plateau potential evoked by applying depolarizing pulse 
	at hyperpolarized level of membrane potential (Gillies 2006, Fig. 10C-D)

	GILLIES CURRENTS

	OTSUKA CURRENTS
	
	- KA peak decreases during burst (height of peaks during AP), as KA decreases the firing frequency also decreases
		- hence KA seems to responsible for rapid repolarization and maintenance of high-frequency firing
	
	- KCa peak increases over about 25 ms (height of peaks during AP), and decreases during last 100 ms of burst
	
	- CaT is the first depolarizing current that rises after release from hyperpolarization and seems to be
	  responsible for initiation of the rebound burst
		- CaT bootstraps burst (bootstraps pos feedback of CaL entry)
		- it runs out of fuel during burst and thus may contribute to ending the burst
			- this is contradicted by burst at regular Vm: there drop in ICaL clearly ends burst
	
	- CaL reaces steady maximum peak after approx. 70 ms into the burst, after CaT is already past its peak
		- hypothesis that it seems responsible for prolonging th burst seems plausible
		- burst seems to go on as long as CaT+CaL remains approx. constant, and burst ends as long as CaT too low

	"""
	# Get electrodes and sections to record from
	dendsec = dends_locs[0][0]
	dendloc = dends_locs[0][1]
	stim1, stim2, stim3 = stims[0], stims[1], stims[2]

	# Set simulation parameters
	dur = 2000
	h.dt = 0.025
	h.celsius = 30 # different temp from paper
	h.v_init = -60 # paper simulations sue default v_init
	gillies.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

	# Set up stimulation (5 mA/cm2 for 80 ms)
	cellarea = np.pi*soma.diam*soma.L # (micron^2)
	I_hyper = -0.17 # hyperpolarize to -70 mV (see fig. 10C)
	I_depol = I_hyper + 0.2 # see fig. 10D: 0.2 nA (=stim.amp) over hyperpolarizing current
	dur_depol = 50 # see fig. 10D, top right
	del_depol = 1000
	burst_time = [del_depol-50, del_depol+200] # empirical

	stim1.delay = 0
	stim1.dur = del_depol
	stim1.amp = I_hyper

	stim2.delay = del_depol
	stim2.dur = dur_depol
	stim2.amp = I_depol

	stim3.delay = del_depol + dur_depol
	stim3.dur = dur - (del_depol + dur_depol)
	stim3.amp = I_hyper

	# Record
	secs = {'soma': soma, 'dend': dendsec}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)

	# Membrane voltages
	traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
	for i, (dend,loc) in enumerate(dends_locs):
		dendname = 'dend%i' % i
		secs[dendname] = dend
		traceSpecs['V_'+dendname] = {'sec':dendname,'loc':loc,'var':'v'}

	# Record ionic currents, open fractions, (in)activation variables
	rec_currents_activations(traceSpecs, 'soma', 0.5)

	# K currents (dendrite)
	traceSpecs['I_KCa_d'] = {'sec':'dend','loc':dendloc,'mech':'sKCa','var':'isKCa'}
	# Ca currents (dendrite)
	traceSpecs['I_CaL_d'] = {'sec':'dend','loc':dendloc,'mech':'HVA','var':'iLCa'}
	traceSpecs['I_CaN_d'] = {'sec':'dend','loc':dendloc,'mech':'HVA','var':'iNCa'}
	traceSpecs['I_CaT_d'] = {'sec':'dend','loc':dendloc,'mech':'CaT','var':'iCaT'}

	# Start recording
	recordStep = 0.05
	recData, _ = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Simulate
	h.tstop = dur
	h.init() # calls finitialize() and fcurrent()
	t0 = h.startsw()
	h.run()
	t1 = h.startsw()
	h.stopsw() # or t1=h.startsw(); runtime = t1-t0
	print("Simulation ran for {:.6f} seconds".format(t1-t0))


	# Plot membrane voltages
	recV = collections.OrderedDict([(k,v) for k,v in recData.items() if k.startswith('V_')]) # preserves order
	figs_vm = analysis.plotTraces(recV, recordStep, yRange=(-80,40), traceSharex=True)
	vm_fig = figs_vm[0]
	vm_ax = figs_vm[0].axes[0]

	# Plot ionic currents, (in)activation variables
	figs, cursors = plot_currents_activations(recData, recordStep)

	# # Soma currents
	# recSoma = collections.OrderedDict([(k,v) for k,v in recData.items() if not k.endswith('_d')])
	# Set fixed ranges for comparison
	# traceYLims = {'V_soma': (-80, 40), 'I_Na': (-0.7, 0.1), 'I_NaL': (-2e-3, 0),
	# 			'I_KDR': (0, 0.1), 'I_Kv3': (-1e-2, 0.1), 'I_KCa': (-1e-3, 8e-3),
	# 			'I_h': (-2e-3, 6e-3), 'I_CaL': (-1.5e-2, 1e-3), 'I_CaN': (-2e-2, 1e-3),
	# 			'I_CaT': (-6e-2, 6e-2)}
	# analysis.plotTraces(recSoma, recordStep, timeRange=burst_time, yRange=traceYLims)

	# # Soma currents (relative)
	# recSoma.pop('V_soma')
	# analysis.cumulPlotTraces(recSoma, recordStep, cumulate=False, timeRange=burst_time)

	# Dendrite currents during burst
	recDend = collections.OrderedDict([(k,v) for k,v in recData.items() if k.endswith('_d')])
	analysis.cumulPlotTraces(recDend, recordStep, timeRange=burst_time)
	
	return recData, figs, cursors

def test_reboundburst(soma, dends_locs, stims):
	"""
	Run rebound burst experiment from original Hoc file (Gillies 2006, Fig. 3-4)

	GILLIES CURRENTS
	
	- same 'CaT bootstraps CaL' mechanism: 
		- small peak in CaT at beginning of burst triggers sharp rise in ICaL with 9x higher peak
		- successive CaL peaks decline in magnitude during burst\
	
	- Ih/HCN sharply declines in peak magnitude over burst to insignificance (approx negexp)
		- burst ends when it died out
		- there is no similar current in Otsuka model

	OTSUKA COMPARISON

		- In Gillies model CaT bootstrap is a small single peak at start of burst, 
		  while in Otsuka model CaT is exp declining peaks
		
		- In Gillies model ICaL is a declining ramp of peaks (approx linear), while 
		  in Otsuka it is slowly activated and inactivated (bugle of peaks)
		
		- Ih/HCN and IKCa which determine shape/recovery of AP are different, resulting
		  in a different evolution of AP shape within a burst

	"""
	# Get electrodes and sections to record from
	dendsec = dends_locs[0][0]
	dendloc = dends_locs[0][1]
	stim1, stim2, stim3 = stims[0], stims[1], stims[2]

	# Set simulation parameters
	dur = 2000
	h.dt = 0.025
	h.celsius = 35 # different temp from paper
	h.v_init = -60 # paper simulations sue default v_init
	gillies.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

	# Set up stimulation
	stim1.delay = 0
	stim1.dur = 500
	stim1.amp = 0.0

	# stim2.delay = 200
	# stim2.dur = 500
	# stim2.amp = -0.11 # -0.25 in full model (hyperpolarize to -75 mV steady state)
	stim2.amp = 0
	
	clamp = h.SEClamp(soma(0.5))
	clamp.dur1 = 0
	clamp.dur2 = 0
	clamp.dur3 = 500
	clamp.amp3 = -75

	stim3.delay = 1000
	stim3.dur = 1000
	stim3.amp = 0.0

	# Record
	secs = {'soma': soma, 'dend': dendsec}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)

	# Membrane voltages
	traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
	traceSpecs['t_global'] = {'var':'t'}
	for i, (dend,loc) in enumerate(dends_locs):
		dendname = 'dend%i' % i
		secs[dendname] = dend
		traceSpecs['V_'+dendname] = {'sec':dendname,'loc':loc,'var':'v'}

	# Record ionic currents, open fractions, (in)activation variables
	rec_currents_activations(traceSpecs, 'soma', 0.5)
	rec_currents_activations(traceSpecs, 'dend', dendloc, ion_species=['ca','k'])

	# Start recording
	recordStep = 0.025
	recData, _ = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Simulate
	h.tstop = dur
	h.init() # calls finitialize() and fcurrent()
	t0 = h.startsw()
	h.run()
	t1 = h.startsw()
	h.stopsw() # or t1=h.startsw(); runtime = t1-t0
	print("Simulation ran for {:.6f} seconds".format(t1-t0))

	# Plot membrane voltages
	recV = collections.OrderedDict([(k,v) for k,v in recData.items() if k.startswith('V_')]) # preserves order
	figs_vm = analysis.plotTraces(recV, recordStep, yRange=(-80,40), traceSharex=True)
	vm_fig = figs_vm[0]
	vm_ax = figs_vm[0].axes[0]

	# Save trace to file
	V_soma = np.array(recData['V_soma'], ndmin=2)
	T_soma = np.array(recData['t_global'], ndmin=2)
	TV_soma = np.concatenate((T_soma, V_soma), axis=0) * 1e-3 # pyelectro expects SI units: seconds, Volts
	fpath = 'C:\\Users\\lkoelman\\cloudstore_m\\simdata\\fullmodel\\rebound_full_SEClamp_Vm_dt25e-3_0-2000_ms.csv'
	np.savetxt(fpath, TV_soma.T, delimiter=',', fmt=['%.3E', '%.7E'])

	# Plot ionic currents, (in)activation variables
	figs_soma, cursors_soma = plot_currents_activations(recData, recordStep, sec_tag='soma')
	figs_dend, cursors_dend = plot_currents_activations(recData, recordStep, sec_tag='dend')
	figs = figs_soma + figs_dend
	cursors = cursors_soma + cursors_dend
	plt.show(block=False)

	# Dendrite currents during burst
	# burst_time = [980, 1120]
	# recDend = collections.OrderedDict([(k,v) for k,v in recData.items() if k.endswith('_d')])
	# analysis.cumulPlotTraces(recDend, recordStep, timeRange=burst_time)

	return recData, figs, cursors

def test_slowbursting(soma, dends_locs, stims):
	""" Test slow rhythmic bursting mode under conditions of constant 
		hyperpolarizing current injection and lower sKCa conductance

	PROTOCOL
	- lower sKCa conductance by 90percent to promote bursting
	- inject constant hyperpolarizing current

	PAPER OBSERVATIONS

	- CaL itself is responsible for burst initiation
			- CaL bootstraps itself in the dendrites
	
	- intra-burst: inter-AP Vm gradually hyperpolarizes due to current injection
	
	- inter-burst: in the dendrites, CaL and CaT (very small due to relatively high Vm) gradually depolarize the membrane
		- slow depolarization continues until majority of CaL channels activated
	"""
	# Get electrodes and sections to record from
	dendsec = dends_locs[0][0]
	dendloc = dends_locs[0][1]
	stim1, stim2, stim3 = stims[0], stims[1], stims[2]

	# Set simulation parameters
	dur = 2000
	h.dt = 0.025
	h.celsius = 37 # for slow bursting experiment
	h.v_init = -60 # paper simulations sue default v_init
	gillies.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

	# Set up stimulation
	stim1.delay = 0
	stim1.dur = dur
	stim1.amp = -0.25

	stim2.delay = 0
	stim2.dur = 0
	stim2.amp = 0.0

	stim3.delay = 0
	stim3.dur = 0
	stim3.amp = 0.0 

	# Record
	secs = {'soma': soma, 'dend': dendsec}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)
	traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
	# Na currents
	traceSpecs['I_Na'] = {'sec':'soma','loc':0.5,'mech':'Na','var':'ina'}
	traceSpecs['I_NaL'] = {'sec':'soma','loc':0.5,'mech':'NaL','var':'inaL'}
	# K currents
	traceSpecs['I_KDR'] = {'sec':'soma','loc':0.5,'mech':'KDR','var':'ik'}
	traceSpecs['I_Kv3'] = {'sec':'soma','loc':0.5,'mech':'Kv31','var':'ik'}
	traceSpecs['I_KCa'] = {'sec':'soma','loc':0.5,'mech':'sKCa','var':'isKCa'}
	traceSpecs['I_h'] = {'sec':'soma','loc':0.5,'mech':'Ih','var':'ih'}
	# K currents (dendrite)
	traceSpecs['dI_KCa'] = {'sec':'dend','loc':dendloc,'mech':'sKCa','var':'isKCa'}
	# Ca currents (soma)
	traceSpecs['I_CaL'] = {'sec':'soma','loc':0.5,'mech':'HVA','var':'iLCa'}
	traceSpecs['I_CaN'] = {'sec':'soma','loc':0.5,'mech':'HVA','var':'iNCa'}
	traceSpecs['I_CaT'] = {'sec':'soma','loc':0.5,'mech':'CaT','var':'iCaT'}
	# Ca currents (dendrite)
	traceSpecs['dI_CaL'] = {'sec':'dend','loc':dendloc,'mech':'HVA','var':'iLCa'}
	traceSpecs['dI_CaN'] = {'sec':'dend','loc':dendloc,'mech':'HVA','var':'iNCa'}
	traceSpecs['dI_CaT'] = {'sec':'dend','loc':dendloc,'mech':'CaT','var':'iCaT'}
	# Start recording
	recordStep = 0.05
	recData, _ = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Simulate
	h.tstop = dur
	if fullmodel:
		h.applyApamin() # lower sKCa conductance by 90 percent
	else:
		gillies.applyApamin(soma, dends)
	h.init()
	h.run()
	if fullmodel:
		h.washApamin() # restore sKCa conductance to original level

	# Analyze
	burst_time = [] # enter burst time
	# Soma currents
	recSoma = collections.OrderedDict()
	for k, v in recData.items():
		if not k.startswith('d'): recSoma[k] = recData[k]
	analysis.plotTraces(recSoma, recordStep)
	# Soma currents (relative)
	recSoma.pop('V_soma')
	analysis.cumulPlotTraces(recSoma, recordStep, cumulate=False)
	# Dendrite currents (relative)
	recDend = collections.OrderedDict()
	for k, v in recData.items():
		if k.startswith('d'): recDend[k] = recData[k]
	analysis.cumulPlotTraces(recDend, recordStep, cumulate=False)
	return recData

def compare_conductance_dist(gnames):
	"""
	Compare conductance distribution betweel full and reduced model

	NOTE: better use ModelView in NEURON GUI instead
	"""

	# Get cells of original model
	or_somaref, or_treeL, or_treeR = analyze_reduction.get_gillies_cells()
	or_secrefs = [or_somaref] + or_treeL + or_treeR
	or_spinysec = next(ref for ref in or_treeR if ref.table_index==8)
	seg = or_spinysec.sec(0.9)
	plateau_gnames = ['gk_sKCa', 'gcaL_HVA', 'gcaN_HVA', 'gcaT_CaT']
	baseline = getattr(seg, 'gk_sKCa')
	for gname in plateau_gnames:
		print("OR ({}) Relative {} = {}".format(seg.sec.name(), gname, getattr(seg, gname)/baseline))

	# Plot conductances in original model
	filter_func = lambda secref: (secref.tree_index==1) and (secref.table_index in (1,3,8))
	label_func = lambda secref: "sec {}".format(secref.table_index)
	for gname in gnames:
		analyze_reduction.plot_tree_ppty(or_somaref, or_secrefs, gname, 
					filter_func, label_func)

	# Get cells of reduced model
	clusters, eq_secs = bush.reduce_bush_sejnowski(delete_old_cells=False)
	eq_refs = [ExtSecRef(sec=sec) for sec in eq_secs]
	eq_somaref = next(ref for ref in eq_refs if ref.sec.name().startswith('soma'))
	eq_spinyref = next(ref for ref in eq_refs if ref.sec.name().startswith('spiny'))
	seg = eq_spinyref.sec(0.9)
	baseline = getattr(seg, 'gk_sKCa')
	for gname in plateau_gnames:
		print("EQ ({}) Relative {} = {}".format(seg.sec.name(), gname, getattr(seg, gname)/baseline))

	# Plot conductances in reduced model
	filter_func = lambda ref: True # plot all sections
	label_func = lambda ref: ref.sec.name() + '_bush'
	for gname in gnames:
		analyze_reduction.plot_tree_ppty(eq_somaref, eq_refs, gname, 
					filter_func, label_func)

# if __name__ == '__main__':
def run_experimental_protocol(reduction_method):
	"""
	Run one of the experiments using full or reduced STN model
	"""
	# Make cell
	soma, dends_locs, stims, allsecs = stn_cell(cellmodel=reduction_method)

	# Manual cell adjustments
	if reduction_method == 4: # manual adjustments to Marasco reduction
		soma.Ra = 2*soma.Ra # correct incorrect calculation for Ra soma cluster
		for sec in h.allsec():
			for seg in sec:
				seg.gpas_STh = 0.75 * seg.gpas_STh
				seg.cm = 3.0 * seg.cm
				seg.gna_NaL = 0.6 * seg.gna_NaL

	# Adjustments to Bush & Sejnowski method
	if reduction_method == 3:
		Rm_factor = math.sqrt(1.)
		Cm_factor = math.sqrt(1.)
		rc_factor = Rm_factor*Cm_factor
		for sec in h.allsec():
			print("Scaling RC of section {} by factor {:.3f}".format(sec.name(), rc_factor))
			for seg in sec:
				seg.gpas_STh = seg.gpas_STh / Rm_factor # multiply Rm is divide gpas
				seg.cm = seg.cm * Cm_factor
				if sec.same(soma):
					seg.gna_NaL = 2.0 * seg.gna_NaL
				else:
					seg.gna_NaL = 1.25 * seg.gna_NaL

	# Adjustments to incremental reduction method
	if reduction_method == 8:
		n_adjusted = 0
		
		for sec in h.allsec():
			
			if sec.name().endswith('soma'):
				print("Skipping soma")
				continue
			
			for seg in sec:
				# seg.gna_NaL = 1.075 * seg.gna_NaL
				seg.gna_NaL = 8e-6 * 1.3 # full model value = uniform 8e-6
				n_adjusted += 1
		
		print("Adjusted parameters of {} segments".format(n_adjusted))


	# Attach duplicate of one tree
	# from redutils import dupe_subtree
	# copy_mechs = {'STh': ['gpas']} # use var gillies_gdict for full copy
	# trunk_copy = duplicate_subtree(h.trunk_0, copy_mechs	, [])
	# trunk_copy.connect(soma, h.trunk_0.parentseg().x, 0)

	# Run experimental protocol
	recData = test_spontaneous(soma, dends_locs, stims)
	# recData = test_reboundburst(soma, dends_locs, stims)
	# recData = test_plateau(soma, dends_locs, stims)
	# recData = test_slowbursting()

	# If run as function, uncomment to make variables available
	globals().update(locals())
	

if __name__ == '__main__':
	run_experimental_protocol(reduction_method=8)
	# compare_conductance_dist(gillies_glist)
	# stn_cell_gillies(resurgent=False)