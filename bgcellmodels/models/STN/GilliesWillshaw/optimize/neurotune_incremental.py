"""
Incremental reduction & optimization of Gillies & Willshaw STN neuron model.

@author Lucas Koelman
@date	03-05-2017
"""
# Python modules
import re
import collections

# Load NEURON
import neuron
from neuron import h

# Third party modules
import numpy as np
import pyelectro, neurotune
from pyelectro import analysis
from neurotune import optimizers
from neurotune import evaluators
from neurotune import controllers

# Own modules
from bgcellmodels.common import analysis as recording
import marasco_foldbased as marasco


# Print which modules we are using
for module in [neuron, pyelectro, neurotune]:
	mod_path = os.path.abspath(module.__file__)
	mod_name = module.__name__
	mod_version = 'N/A' if not hasattr(module, '__version__') else module.__version__
	print("Using module {} - version {} at {}".format(mod_name, mod_version, mod_path))


class Protocol:
	""" Experimental protocols """
	SPONTANEOUS = 0
	REBOUND = 2
	PLATEAU = 3

# Saved voltage trace for each protocol
protocol_Vm_paths = { # ADJUSTME
	Protocol.SPONTANEOUS:"C:\\Users\\lkoelman\\cloudstore_m\\simdata\\fullmodel\\spont_fullmodel_Vm_dt25e-3_0ms_2000ms.csv",
	Protocol.REBOUND:"C:\\Users\\lkoelman\\cloudstore_m\\simdata\\fullmodel\\rebound_full_SEClamp_Vm_dt25e-3_0-2000_ms.csv",
}

# Time interval where voltage traces should be compared for each protocol
protocol_intervals = { # ADJUSTME
	Protocol.SPONTANEOUS:	(200.0, 1000.0), # representative interval
	Protocol.REBOUND:		(500.0, 700.0), # interval of rebound burst
	Protocol.PLATEAU:		(200.0, 1000.0),
}

class Simulation(object):
	""" 
	Helper class for the Controller object to run simulations and generate
	raw simulation data
	"""

	def __init__(self, rec_section, dt=0.025):
		self.sec = rec_section
		self.dt = dt

	def set_aCSF(self, req):
		""" Set global initial ion concentrations (artificial CSF properties) """
		if req == 3: # Beurrier et al (1999)
			h.nai0_na_ion = 15
			h.nao0_na_ion = 150

			h.ki0_k_ion = 140
			h.ko0_k_ion = 3.6

			h.cai0_ca_ion = 1e-04
			h.cao0_ca_ion = 2.4

			h('cli0_cl_ion = 4') # self-declared Hoc var
			h('clo0_cl_ion = 135') # self-declared Hoc var

		if req == 4: # Bevan & Wilson (1999)
			h.nai0_na_ion = 15
			h.nao0_na_ion = 128.5

			h.ki0_k_ion = 140
			h.ko0_k_ion = 2.5

			h.cai0_ca_ion = 1e-04
			h.cao0_ca_ion = 2.0

			h('cli0_cl_ion = 4')
			h('clo0_cl_ion = 132.5')

		if req == 0: # NEURON's defaults
			h.nai0_na_ion = 10
			h.nao0_na_ion = 140

			h.ki0_k_ion = 54
			h.ko0_k_ion = 2.5

			h.cai0_ca_ion = 5e-05
			h.cao0_ca_ion = 2

			h('cli0_cl_ion = 0')
			h('clo0_cl_ion = 0')

	def set_recording(self, recordStep=0.025):
		"""
		Set up recording Vectors to record from relevant pointers
		"""
		# Named sections to record from
		secs = {'soma': self.sec}

		# Specify traces
		traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)
		traceSpecs['V_soma'] = {'sec':'soma', 'loc':0.5, 'var':'v'}
		traceSpecs['t_global'] = {'var':'t'}

		# Make trace record Vectors
		self.rec_dt = recordStep
		self.recData = recording.recordTraces(secs, traceSpecs, self.rec_dt)

	def show(self, candidate):
		"""
		Plot the result of the simulation once it's been intialized
		"""
		from matplotlib import pyplot as plt

		v_rec = np.array(self.recData['V_soma'])
		t_rec = np.array(self.recData['t_global'])

		plt.plot(t_rec, v_rec)
		cand_params = ["{:.2f}".format(par) for par in candidate.values()]
		plt.title("Simulation for cand: [{}]".format(';'.join(cand_params)))
		plt.xlabel("Time [ms]")
		plt.ylabel("Voltage [mV]")

		plt.show(block=True)

	def simulate_protocol(self, protocol):
		""" Simulate the given experimental protocol """

		if protocol == Protocol.SPONTANEOUS:
			# Spontaneous firing, no stimulation
			h.dt = self.dt
			h.celsius = 37 # different temp from paper (fig 3B: 25degC, fig. 3C: 35degC)
			h.v_init = -60 # paper simulations use default v_init
			self.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

		elif protocol == Protocol.PLATEAU:
			raise NotImplementedError("Implement PLATEAU protocol")

		elif protocol == Protocol.REBOUND:
			# Rebound burst: response to brief hyperpolarizing pulse
			h.dt = self.dt
			h.celsius = 35 # different temp from paper
			h.v_init = -60 # paper simulations sue default v_init
			self.set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

			# Stimulation: hyperpolarize to -75, wait for steady state, then release
			clamp = h.SEClamp(self.sec(0.5)) # NOTE: use SEClamp instead of IClamp to hyperpolarize to same level independent of passive parameters
			clamp.dur1 = 0
			clamp.dur2 = 0
			clamp.dur3 = 500
			clamp.amp3 = -75

		else:
			raise Exception("Unrecognized protocol {0}".format(protocol))

		# Set up recording vectors
		self.set_recording()

		# Run simulation
		h.init() # calls finitialize() and fcurrent()
		t_interval = protocol_intervals[protocol]
		neuron.run(t_interval[1] + 5.0) # do not need to simulate further than end of protocol

		# Return time and voltage traces
		# NOTE: do not extract segment, time axes must match with target trace!
		#       (extraction is done via Evaluator.analysis_start/end_time)
		v_rec = np.array(self.recData['V_soma'])
		t_rec = np.array(self.recData['t_global'])
		return t_rec, v_rec


class STNCellController(object):
	""" 
	The Controller class maps candidate parameter sets to raw data
	by setting up the model for the candidate parameter set and 
	simulating it

	The only interface it must adhere to is providing a run()
	method with the following signature:

	run(candidates: list(list(float)), parameters:list(string)) -> (t_trace, v_trace)
	"""

	def __init__(self, protocol, reduced_model_path):
		""" Make new Controller that builds and runs model starting from
			the given base model.

		@type	protocol			Protocol enum instance
		@param	protocol			Experimental protocol that should be simulated

		@type	reduced_model_path	string
		@param	reduced_model_path	Path to pickle file containing the serialized
									clusters of the reduced model.
		"""
		self.target_protocol = protocol
		self.allsecrefs = [] # SectionRefs to existing model sections
		self.n_passes_current = 0
		self.consistent_seg_props = False

	def run(self, candidates, parameters):
		"""
		Run simulation for each candidate

		NOTE: this is the only function where its signature is fixed
			  by an interface/protoype of the neurotune package

		@type	candidates	list(list(float))
		@param	candidates	list of candidate parameter sets (each candidate
							is a list of parameter values)

		@type	parameters	list of parameter names
		@param	parameters	our self-defined list of parameter names: these
							will be used to build the model in some way

		@return				list of (t_trace, v_trace) for each candidate,
							i.e. a list(tuple(list(float), list(float)))
		"""
		traces = []
		for candidate in candidates:
			# Build simulation for candidate
			cand_params = dict(zip(parameters,candidate))
			soma_sec, dend_secs = self.build_candidate(cand_params)
			sim = Simulation(soma_sec)

			# Simulate candidate
			t, v = sim.simulate_protocol(self.target_protocol)
			traces.append([t,v])
			# sim.show(cand_params)

		return traces

	def initialize_model(self):
		"""
		Initialize model before reduction.
		"""
		# Load full model to be reduced (Gillies & Willshaw STN)
		for sec in h.allsec():
			if not sec.name().startswith('SThcell'):
				h.delete_section() # delete existing cells
		if not hasattr(h, 'SThcells'):
			h.xopen("gillies_cell_singleton.hoc")

		# Make sections accesible by name and index
		somaref = ExtSecRef(sec=h.SThcell[0].soma)
		dendLrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend0] # 0 is left tree
		dendRrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend1] # 1 is right tree
		allsecrefs = [somaref] + dendLrefs + dendRrefs
		orsecrefs = allsecrefs # SectionRef to original sections

		# Get references to root sections of the 3 identical trees
		dendR_root = getsecref(h.SThcell[0].dend1[0], dendRrefs)
		dendL_juction = getsecref(h.SThcell[0].dend0[0], dendLrefs)
		dendL_upper_root = getsecref(h.SThcell[0].dend0[1], dendLrefs)
		dendL_lower_root = getsecref(h.SThcell[0].dend0[2], dendLrefs)

		# Assign indices used in Gillies code (sth-data folder)
		somaref.tree_index = -1
		somaref.table_index = 0
		for j, dendlist in enumerate((dendLrefs, dendRrefs)):
			for i, secref in enumerate(dendlist):
				secref.tree_index = j # left tree is 0, right is 1
				secref.table_index = i+1 # same as in /sth-data/treeX-nom.dat

		# Calculate section path properties for entire tree
		for secref in allsecrefs:
			# Assign path length, path resistance, electrotonic path length to each segment
			redtools.sec_path_props(secref, f_lambda, gleak_name)

		# Store segment properties along interpolation path
		interp_tree_id = 1
		interp_table_ids = (1,3,8)
		path_secs = [secref for secref in orsecrefs if (secref.tree_index == interp_tree_id and 
														secref.table_index in interp_table_ids)]
		sec_assigned = ['pathL0', 'pathL1', 'pathri0', 'pathri1', 'pathLelec0', 'pathLelec1']
		seg_assigned = ['pathL_seg', 'pathri_seg', 'pathL_elec']
		path_props = [redtools.get_sec_props_obj(ref, mechs_chans, seg_assigned, sec_assigned) for ref in path_secs]

		# Save model information
		self.n_passes_current = 0

		# Hack to turn local variables into instance variables (add them to the object's namespace)
		loc_vars = dict(locals())
		loc_vars.pop('self')
		self.__dict__.update(locals())

	def reduce_model(self, n_passes):
		"""
		Do N incremental reduction passes.
		"""
		# Hack to bring instance variables (self's namespace) into local namespace
		locals().update(self.__dict__)
		raise Exception("TESTME: check if self.__dict__ has variabled that need to be purged before updating locals()")

		# Start iterative collapsing procedure
		for i_pass in range(n_passes):
			# Check that we haven't collapsed too many levels
			if not (dendL_upper_root.exists() and dendL_lower_root.exists()):
				logger.warning("Maximum number of collapses reached: Cannot collapse past trunk sections.")
				break

			############################################################################
			# 0. Pre-clustering: calculate properties

			# Assign attributes used in calculation
			assign_sth_indices(somaref, allsecrefs)
			assign_attributes(somaref, allsecrefs, {'max_passes': 100})

			# Assign Strahler numbers
			clutools.assign_topology_attrs(dendR_root, allsecrefs)
			clutools.assign_topology_attrs(dendL_upper_root, allsecrefs)
			clutools.assign_topology_attrs(dendL_lower_root, allsecrefs)
			dendL_juction.order = 1
			dendL_juction.level = 0
			dendL_juction.strahlernumber = dendL_upper_root.strahlernumber+1
			somaref.order = 0
			somaref.level = 0
			somaref.strahlernumber = dendL_juction.strahlernumber

			############################################################################
			# 1. Merge a number of Y-sections

			# Find collapsable branch points
			target_Y_secs = marasco.find_collapsable(allsecrefs, i_pass, 'highest_level', zips_per_pass)

			# Do collapse operation at each branch points
			clusters = marasco.calc_collapses(target_Y_secs, i_pass, allsecrefs)

			############################################################################
			# 2. Create equivalent sections
			# Mark Sections
			for secref in allsecrefs:
				secref.is_substituted = False
				secref.is_deleted = False

			eq_refs, newsecrefs = marasco.sub_equivalent_Y_sections(clusters, allsecrefs, path_props,
											interp_prop='path_L', interp_method='linear_neighbors', 
											gbar_scaling='area')

			allsecrefs = newsecrefs # prepare for next iteration

			############################################################################
			# 3. Finalize
			
			# Set ion styles
			for sec in h.allsec(): # makes each section the CAS
				h.ion_style("na_ion",1,2,1,0,1)
				h.ion_style("k_ion",1,2,1,0,1)
				h.ion_style("ca_ion",3,2,1,1,1)

		# Save current model state
		self.allsecrefs = newsecrefs
		self.n_passes_current += n_passes

		# Print tree structure
		logger.info("Equivalent tree topology:")
		if logger.getEffectiveLevel() >= logging.DEBUG:
			h.topology()

	def reset_model(self):
		"""
		Reset segment properties to those stored in secref.or_seg_props
		"""
		for ref in self.allsecrefs:
			if not ref.exists():
				continue
			for j_seg, seg in enumerate(ref.sec):
				for pname, pval in ref.or_seg_props[j_seg]:
					seg.__setattr__(pname, pval)
		self.consistent_seg_props = True

	def save_model_state(self):
		"""
		Save current segment properties for later resetting
		"""
		for ref in self.allsecrefs:
			redtools.store_seg_props(ref, mechs_chans, attr_name='or_seg_props')
		self.consistent_seg_props = True

	def build_candidate(self, cand_params):
		"""
		Build simulation for given candidate parameter set

		@post		self.eq_secs is updated to contain the NEURON Sections
					representing the given candidate

		@return		tuple (soma section, list(dendrite sections))
		"""
		# Initialize the model
		n_passes_target = cand_params['n_reduction_passes']
		if not any(self.allsecrefs):
			self.initialize_model()
			self.reduce_model(n_passes_target)
			self.save_model_state()
			
		# Reduce further if required
		if n_passes_target >= self.n_passes_current:
			self.reset_model()
			self.reduce_model(n_passes_target - self.n_passes_current)
			self.save_model_state()

		# Reset the model
		if not self.consistent_seg_props:
			self.reset_model()

		# Pattern matching for gbar factors
		scale_prefix = r'^gbar_sca_'
		scale_pattern = re.compile(scale_prefix)
		gmin_prefix = r'^gbar_min_'
		gmin_pattern = re.compile(gmin_prefix)
		gmax_prefix = r'^gbar_max_'
		gmax_pattern = re.compile(gmax_prefix)

		# Adapt model according to candidate parameters
		for par_name, par_value in cand_params.items():
			scale_match = re.search(scale_pattern, par_name)
			gmin_match = re.search(gmin_pattern, par_name)
			gmax_match = re.search(gmax_pattern, par_name)

			if par_name == 'soma_cm_factor':
				for seg in soma_sec:
					seg.cm = seg.cm * par_value

			elif par_name == 'soma_Rm_factor':
				for seg in soma_sec:
					gleak_val = getattr(seg, gleak_name) / par_value
					setattr(seg, gleak_name, gleak_val)

			elif par_name == 'soma_Ra':
				soma_sec.Ra = par_value

			elif par_name == 'soma_diam_factor':
				for seg in soma_sec:
					seg.diam = seg.diam * par_value

			elif par_name == 'dend_cm_factor':
				for sec in dend_secs:
					for seg in sec:
						seg.cm = seg.cm * par_value

			elif par_name == 'dend_Rm_factor':
				for sec in dend_secs:
					for seg in sec:
						gleak_val = getattr(seg, gleak_name) / par_value
						setattr(seg, gleak_name, gleak_val)

			elif par_name == 'dend_Ra':
				for sec in dend_secs:
					sec.Ra = par_value

			elif par_name == 'dend_diam_factor':
				for sec in dend_secs:
					for seg in sec:
						seg.diam = seg.diam * par_value

			elif scale_match:
				prefix_suffix = re.split(scale_prefix, par_name)
				gname = prefix_suffix[1]
				secs = self.eq_secs if self.gbar_adjust_allsec else dend_secs
				for sec in secs:
					for seg in sec:
						gval = getattr(seg, gname) * par_value
						setattr(seg, gname, gval)
			else:
				raise Exception("Unrecognized parameter '{}' with value <{}>".format(
								par_name, par_value))

		self.consistent_seg_props = False

		return soma_sec, dend_secs

class CustomIClampEvaluator(evaluators.__Evaluator):
	"""
	Evaluate the fitness value of candidates by calculating and comparing
	metrics on the simulation results.

	Based on neurotune.IClampEvaluator
	"""
	def __init__(self,
				 analysis_start_time,
				 controller,
				 analysis_end_time,
				 target_data_path,
				 parameters,
				 analysis_var,
				 weights,
				 targets=None,
				 automatic=False):

		super(CustomIClampEvaluator, self).__init__(parameters,
											  weights,
											  targets,
											  controller)
	  
		self.analysis_start_time = analysis_start_time
		self.analysis_end_time = analysis_end_time
		self.target_data_path = target_data_path
		self.analysis_var = analysis_var

		print('target data path in evaluator:' + target_data_path)
		
		if automatic == True:
			t , v_raw = analysis.load_csv_data(target_data_path)
			v = np.array(v_raw)

			v_smooth = list(analysis.smooth(v))

			ic_analysis = analysis.IClampAnalysis(
							v_smooth,
							t,
							analysis_var,
							start_analysis = analysis_start_time,
							end_analysis = analysis_end_time,
							show_smoothed_data = False
						) 

			ic_analysis.analyse()

			self.targets = ic_analysis.analysis_results

			print('Obtained targets are:')
			print(self.targets)
		
	def evaluate(self,candidates,args):
		
		print("\n>>>>>  Evaluating: ")
		for cand in candidates: print(">>>>>       %s"%cand)
		
		simulations_data = self.controller.run(candidates,
											   self.parameters)

		fitness = []
		
		for times, samples in simulations_data:
			# Calculate metrics for each trace
			data_analysis = analysis.IClampAnalysis(
								samples,
								times,
								self.analysis_var,
								start_analysis = self.analysis_start_time,
								end_analysis = self.analysis_end_time,
								target_data_path = self.target_data_path,
								show_smoothed_data = False,
							)

			
			try:
				data_analysis.analyse()
			except:
				data_analysis.analysable_data = False
				
				
			fitness_value = self.evaluate_fitness(
								data_analysis,
								self.targets,
								self.weights,
								cost_function=evaluators.normalised_cost_function
							)
			fitness.append(fitness_value)

			print('Fitness: %s\n'%fitness_value)
			
		return fitness
	

	def evaluate_fitness(self,
						 data_analysis,
						 target_dict={},
						 target_weights=None,
						 cost_function=evaluators.normalised_cost_function):
		"""
		Return the estimated fitness of the data, based on the cost function being used.
			:param data_analysis:	IClampAnalysis instance
			:param target_dict:		key-value pairs for targets
			:param target_weights: 	key-value pairs for target weights
			:param cost_function: 	cost function (callback) to assign individual targets sub-fitness.
		"""
	
		# calculate max fitness value (TODO: there may be a more pythonic way to do this..)
		worst_cumulative_fitness=0
		for target in target_dict.keys():
			if target_weights == None: 
				target_weight = 1
			else:
				if target in target_weights.keys():
					target_weight = target_weights[target]
				else:
					target_weight = 1.0
	
			worst_cumulative_fitness += target_weight

		#if we have 1 or 0 peaks we won't conduct any analysis
		if data_analysis.analysable_data == False:
			print('Data is non-analysable')
			return worst_cumulative_fitness

		else:
			fitness = 0

			for target in target_dict.keys():

				target_value=target_dict[target]

				if target_weights == None: 
					target_weight = 1
				else:
					if target in target_weights.keys():
						target_weight = target_weights[target]
					else:
						target_weight = 1.0
				if target_weight > 0:
					value = data_analysis.analysis_results[target]
					#let function pick Q automatically
					inc = target_weight*cost_function(value,target_value)
					fitness += inc

					print('Target %s (weight %s): target val: %s, actual: %s, fitness increment: %s'%(target, target_weight, target_value, value, inc))

			return fitness

def optimization_routine():
	"""
	Main method for the optimization routine.

	Performs the role of the Optimizer (e.g. see neurotune.CustomOptimizerA).

	HOWTO adjust for optimization:
		- Simulate protocol with full model and save trace with same time step
		- Adjust protocol for reduced model (stim timing + adjust stim current to polarize to same level)
		- Adjust optimization parameters (see ADJUSTME tags):
			- adjust parameters to optimize ('chromosome')
			- adjust seed candidates
			- adjust targets and their weights based on protocol
	"""

	# Make a controller to simulate candidates
	reduced_model_path = "C:\\Users\\lkoelman\\cloudstore_m\\simdata\\bush_sejnowski\\stn_reduced_bush.p"
	# Voltage trace of original model for comparison
	target_protocol = Protocol.REBOUND # ADJUSTME
	target_Vm_path = protocol_Vm_paths[target_protocol]

	# Create controller to run simulations
	stn_controller = STNCellController(target_protocol, reduced_model_path)
	stn_controller.gbar_adjust_allsec = False # whether gbar scaling wil apply to all sections (ADJUSTME)

	# Parameters that constitute a candidate ('DNA') and their bounds
	# NOTE: based on fitting routine described in Gillies & Willshaw (2006)
	passive_params_bounds = {
		# soma properties
		'soma_cm_factor': 		(1.0,2.0),
		'soma_Rm_factor': 		(0.5,5.0),
		'soma_Ra':				(100.,300.), # 150 in full model
		'soma_diam_factor':		(0.5,1.0),
		# dendrite properties
		'dend_cm_factor':		(1.0,2.0),
		'dend_Rm_factor':		(0.5,2.0),
		'dend_Ra':				(100.,300.), # 150 in full model
		'dend_diam_factor':		(0.5,2.0),
	}

	active_params_bounds = {
		# soma properties properties
		'soma_gNaL_factor':		(0.5,2.0), # scales constant gNaL distribtuion
		# dendrite properties
		'gbar_min_gk_Ih':		(0.5,2.0), # distal linear dist
		'gbar_max_gk_Ih':		(0.5,2.0),
		'gbar_min_gk_KDR':		(0.5,2.0), # distal step dist
		'gbar_max_gk_KDR':		(0.5,2.0),
		'gbar_min_gk_Kv31':		(0.5,2.0), # proximal linear dist
		'gbar_max_gk_Kv31':		(0.5,2.0),
		'gbar_min_gk_sKCa':		(0.5,2.0), # distal double step dist
		'gbar_max_gk_sKCa':		(0.5,2.0),
		'gbar_min_gcaL_HVA':	(0.5,2.0), # distal linear dist
		'gbar_max_gcaL_HVA':	(0.5,2.0),
		'gbar_min_gcaN_HVA':	(0.5,2.0), # proximal linear dist
		'gbar_max_gcaN_HVA':	(0.5,2.0),
		'gbar_min_gcaT_CaT':	(0.5,2.0), # distal linear dist
		'gbar_max_gcaT_CaT':	(0.5,2.0),
		# dendrite properties: scaling factors
		'gbar_sca_gk_Ih':		(0.5,2.0), # distal linear dist
		'gbar_sca_gk_KDR':		(0.5,2.0), # distal step dist
		'gbar_sca_gk_Kv31':		(0.5,2.0), # proximal linear dist
		'gbar_sca_gk_sKCa':		(0.5,2.0), # distal double step dist
		'gbar_sca_gcaL_HVA':	(0.5,2.0), # distal linear dist
		'gbar_sca_gcaN_HVA':	(0.5,2.0), # proximal linear dist
		'gbar_sca_gcaT_CaT':	(0.5,2.0), # distal linear dist
		'gbar_sca_gna_NaL':		(0.5,1.0), # constant dist
		# 'gbar_sca_gna_Na_factor', # constant dist (negligibly small)
	}

	all_params_bounds = {}
	all_params_bounds.update(passive_params_bounds)
	all_params_bounds.update(active_params_bounds)

	# Parameters to tune spontaneous firing
	spont_params = [
		# soma passive properties
		'soma_diam_factor',
		'soma_cm_factor',
		# 'soma_Ra', # 150 in full model
		# dendrite passive properties
		'dend_cm_factor',
		'dend_Rm_factor',
		'dend_Ra', # 150 in full model
		'dend_diam_factor',
		# active properties
		'gbar_sca_gna_NaL',
	]
	# Parameters to tune bursts
	burst_params = [
		# soma properties
		'soma_diam_factor',
		'soma_cm_factor',
		# dendrite passive properties
		'dend_cm_factor',
		'dend_Rm_factor',
		'dend_Ra',
		'dend_diam_factor',
		# dendrite active properties
		'gbar_sca_gk_sKCa', # distal double step dist
		'gbar_sca_gcaL_HVA', # distal linear dist
		'gbar_sca_gcaT_CaT', # distal linear dist
	]
	

	# Seeds (initial candidates), e.g. from previous optimization
	spont_seeds = [ # soma_diam, soma_cm, dend_cm, dend_Rm, dend_Ra, dend_diam, gna_NaL
		[1.0, 1.0, 1.0, 1.0, 150., 1.0, 1.0], # default after model reduction
		[0.5573267195761666, 1.252340068894955, 1.8412707861542312, 2.0, 228.62780124436583, 0.6617425835741709, 1.0],
		[0.563915223257359, 1.4561344744499622, 2.0, 2.0, 228.0626727668688, 0.6617425835741709, 1.0],
	]

	burst_seeds = [ # soma_diam, soma_cm, dend_cm, dend_Rm, dend_Ra, dend_diam, gsKCa, gcaL, gcaT
		[1.0, 1.0, 1.0, 1.0, 150., 1.0, 1.0, 1.0, 1.0], # default after model reduction
		[0.5573267195761666, 1.252340068894955, 1.8412707861542312, 2.0, 228.62780124436583, 0.6617425835741709, 1.0, 1.0, 1.0],
		[0.563915223257359, 1.4561344744499622, 2.0, 2.0, 228.0626727668688, 0.6617425835741709, 1.0, 1.0, 1.0],
		[0.5040251966091115, 1.282989976, 1.30888976, 1.143723054, 195.2667889, 0.782315645, 1.063518297, 1.644093284, 1.4511723121050117],
		[0.5166409889357284, 1.0, 1.0046280491823951, 1.9313327747357878, 193.88153564188488, 0.7810045138119329, 1.0052884459348173, 1.4765306846140291, 1.4500118248990195] 
	]

	# Parameters for calculation of voltage trace metrics
	trace_analysis_params = {
		'peak_delta':		1e-4, # the value by which a peak or trough has to exceed its neighbours to be considered outside of the noise
		'baseline':			-10., # voltage at which AP width is measured
		'dvdt_threshold':	0, # used in PPTD method described by Van Geit 2007
	}

	# Weight for components of error measure (= targets)
	# NOTE: see metrics in http://pyelectro.readthedocs.io/en/latest/pyelectro.html
	# NOTE: default targets are all keys of in IClampAnalysis.analysis_results dict
	spont_error_weights = {
		# Spike timing/frequency related
		'first_spike_time': 1.0,		# time of first AP
		'max_peak_no': 2.0,				# number of AP peaks
		'min_peak_no': 1.0,				# number of AP throughs
		'spike_frequency_adaptation': 1.0,	# slope of exp fit to initial & final frequency
		'trough_phase_adaptation': 1.0,	# slope of exp fit to phase of first and last through
		'mean_spike_frequency': 2.0,	# mean AP frequency
		'interspike_time_covar': 1.0,	# coefficient of variation of ISIs
		'average_maximum': 1.0,			# average value of AP peaks
		'average_minimum': 1.0,			# average value of AP throughs
		# Spike shape related
		'spike_broadening': 1.0,		# ratio first AP width to avg of following AP widths
		'spike_width_adaptation': 1.0,	# slope of exp fit to first & last AP width
		'peak_decay_exponent': 1.0,		# Decay of AP peak height
		'trough_decay_exponent': 1.0,	# Decay of AP through height
		'pptd_error': 2.0				# Phase-plane trajectory density (see Van Geit (2008))
	}

	burst_error_weights = {
		# Spike timing/frequency related
		'first_spike_time': 0.0,		# time of first AP
		'max_peak_no': 1.0,				# number of AP peaks
		'min_peak_no': 1.0,				# number of AP throughs
		'spike_frequency_adaptation': 1.0,	# slope of exp fit to initial & final frequency
		'trough_phase_adaptation': 1.0,	# slope of exp fit to phase of first and last through
		'mean_spike_frequency': 0.0,	# mean AP frequency
		'interspike_time_covar': 1.0,	# coefficient of variation of ISIs
		'average_maximum': 1.0,			# average value of AP peaks
		'average_minimum': 1.0,			# average value of AP throughs
		# Spike shape related
		'spike_broadening': 0.0,		# ratio first AP width to avg of following AP widths
		'spike_width_adaptation': 1.0,	# slope of exp fit to first & last AP width
		'peak_decay_exponent': 1.0,		# Decay of AP peak height
		'trough_decay_exponent': 1.0,	# Decay of AP through height
		'pptd_error': 1.0				# Phase-plane trajectory density (see Van Geit (2008))
	}

	# Final parameters for current optimization (ADJUSTME)
	final_params = burst_params
	final_error_weights = burst_error_weights
	final_seeds = burst_seeds

	# Make evaluator to map candidate parameter sets to fitness values
	stn_evaluator = CustomIClampEvaluator(
						controller = stn_controller,
						analysis_start_time = protocol_intervals[target_protocol][0],
						analysis_end_time = protocol_intervals[target_protocol][1],
						target_data_path = target_Vm_path,
						parameters = final_params,
						analysis_var = trace_analysis_params,
						weights = final_error_weights,
						targets = None, # if not None: provide self-computed metrics
						automatic = True # if automatic: metrics computed from target trace
					)

	#make an optimizer
	min_constraints = [all_params_bounds[par][0] for par in final_params]
	max_constraints = [all_params_bounds[par][1] for par in final_params]
	# The optimizer creates an inspyred.ec.EvolutionaryComputation() algorithm
	# and calls its evolve() method (see docs at http://pythonhosted.org/inspyred/reference.html#inspyred.ec.EvolutionaryComputation)
	my_optimizer = optimizers.CustomOptimizerA(
						max_constraints,
						min_constraints,
						stn_evaluator,
						population_size = 15, # initial number of individuals/candidates
						max_evaluations = 400,
						num_selected = 3, # how many individuals should become parents
						num_offspring = 6, # total number of offspring (default = pop size)
						num_elites = 1,
						seeds = final_seeds #  iterable collection of candidate solutions to include in the initial population
					)

	#run the optimizer
	my_optimizer.optimize()

if __name__ == '__main__':
	optimization_routine()
