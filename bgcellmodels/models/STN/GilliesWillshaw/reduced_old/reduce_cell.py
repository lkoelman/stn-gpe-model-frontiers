"""
Object-oriented interface for various compartmental cell reduction methods.

@author Lucas Koelman
@data	24-08-2017
"""

from enum import Enum, unique

from bgcellmodels.common.nrnutil import ExtSecRef, getsecref
from neuron import h

import marasco_foldbased as marasco
import stratford_folding as stratford
import mapsyn

# logging of DEBUG/INFO/WARNING messages
import logging
logging.basicConfig(format='%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logname = 'folding'
logger = logging.getLogger(logname) # create logger for this module


################################################################################
# Reduction classes
################################################################################


@unique
class ReductionMethod(Enum):
	Rall = 0
	Stratford = 1			# Stratford, K., Mason, A., Larkman, A., Major, G., and Jack, J. J. B. (1989) - The modelling of pyramidal neurones in the visual cortex
	BushSejnowski = 2		# Bush, P. C. & Sejnowski, T. J. Reduced compartmental models of neocortical pyramidal cells. Journal of Neuroscience Methods 46, 159-166 (1993).
	Marasco = 3				# Marasco, A., Limongiello, A. & Migliore, M. Fast and accurate low-dimensional reduction of biophysically detailed neuron models. Scientific Reports 2, (2012).


class FoldReduction(object):
	"""
	Class grouping methods and data used for reducing
	a compartmental cable model of a NEURON cell.
	"""

	# For each step in reduction process: a dict mapping ReductionMethod -> (func, arg_names)
	
	_PREPROC_FUNCS = {
		ReductionMethod.Marasco:	(marasco.preprocess_impl, []),
		ReductionMethod.Stratford:	(stratford.preprocess_impl, []),
	}

	_PREP_FOLD_FUNCS = {
		ReductionMethod.Marasco:	(marasco.prepare_folds_impl, []),
		ReductionMethod.Stratford:	(stratford.prepare_folds_impl, []),
	}

	_CALC_FOLD_FUNCS = {
		ReductionMethod.Marasco:	(marasco.calc_folds_impl, []),
		ReductionMethod.Stratford:	(stratford.calc_folds_impl, []),
	}

	_MAKE_FOLD_EQ_FUNCS = {
		ReductionMethod.Marasco:	(marasco.make_folds_impl, []),
		ReductionMethod.Stratford:	(stratford.make_folds_impl, []),
	}

	_POSTPROC_FUNCS = {
		ReductionMethod.Marasco:	(marasco.postprocess_impl, []),
		ReductionMethod.Stratford:	(stratford.post_process_impl, []),
	}

	# Make accessible by step
	_REDUCTION_STEPS_FUNCS = {
		'preprocess':		_PREPROC_FUNCS,
		'prepare_fold':		_PREP_FOLD_FUNCS,
		'calculate_fold':	_CALC_FOLD_FUNCS,
		'make_fold':		_MAKE_FOLD_EQ_FUNCS,
		'postprocess':		_POSTPROC_FUNCS,
	}


	def register_reduction_step(self, method, step, func, kwarg_list=None):
		"""
		Register function to perform reduction step 'step' for reduction method 'method'/

		@param	method			ReductionMethod member

		@param	step			ReductionStep member
		"""
		if kwarg_list is None:
			kwarg_list = []
		self._REDUCTION_STEPS_FUNCS[step][method] = (func, kwarg_list)


	def __init__(self, soma_secs, dend_secs, fold_root_secs, method):
		"""
		Initialize reduction of NEURON cell with given root Section.

		@param	soma_secs		list of root Sections for the cell (up to first branch points).
								This list must not contain any Section in dend_secs

		@param	dend_secs		list of dendritic section for each trunk section,
								i.e. the lists are non-overlapping / containing
								unique Sections

		@param	method			ReductionMethod instance
		"""

		# Parameters for reduction method (set by user)
		self._REDUCTION_PARAMS = dict((method, {}) for method in list(ReductionMethod))

		# Reduction method
		self.active_method = method
		self._mechs_gbars_dict = None

		# Find true root section
		first_root_sec = soma_secs[0]
		first_root_ref = ExtSecRef(sec=soma_secs[0])
		root_sec = first_root_ref.root # pushes CAS
		
		# Save unique sections
		self._soma_refs = [ExtSecRef(sec=sec) for sec in soma_secs]
		self._dend_refs = [ExtSecRef(sec=sec) for sec in dend_secs]

		# Save ion styles
		ions = ['na', 'k', 'ca']
		self._ion_styles = dict(((ion, h.ion_style(ion+'_ion')) for ion in ions))
		h.pop_section() # pops CAS

		# Save root sections
		self._root_ref = getsecref(root_sec, self._soma_refs)
		allsecrefs = self.all_sec_refs
		self._fold_root_refs = [getsecref(sec, allsecrefs) for sec in fold_root_secs]

		# Set NetCon to be mapped
		self._syns_tomap = []
		self._map_syn_info = []

	@property
	def all_sec_refs(self):
		"""
		Get list of SectionRef to all sections.
		"""
		return list(self._soma_refs) + list(self._dend_refs)


	@property
	def soma_refs(self):
		"""
		Get list of SectionRef to somatic sections.
		"""
		return self._soma_refs
	

	@property
	def dend_refs(self):
		"""
		Get list of SectionRef to dendritic sections.
		"""
		return self._dend_refs

	def set_syns_tomap(self, syns):
		"""
		Set synapses to map.

		@param syns		list(SynInfo)
		"""
		self._syns_tomap = syns

	@property
	def map_syn_info(self):
		"""
		Synapse properties before and after mapping (electrotonic properties etc.)
		"""
		return self._map_syn_info
	

	def get_mechs_gbars_dict(self):
		"""
		Get dictionary of mechanism names and their conductances.
		"""
		return self._mechs_gbars_dict
	

	def set_mechs_gbars_dict(self, val):
		"""
		Set mechanism names and their conductances
		"""
		self._mechs_gbars_dict = val
		self.gbar_names = [gname+'_'+mech for mech,chans in val.iteritems() for gname in chans]
		self.active_gbar_names = list(self.gbar_names)
		self.active_gbar_names.remove(self.gleak_name)

	# make property
	mechs_gbars_dict = property(get_mechs_gbars_dict, set_mechs_gbars_dict)


	def update_refs(self, soma_refs=None, dend_refs=None):
		"""
		Update Section references after sections have been created/destroyed/substituted.

		@param soma_refs	list of SectionRef to at least all new soma sections
							(may also contain existing sections)

		@param dend_refs	list of SectionRef to at least all new dendritic sections
							(may also contain existing sections)
		"""
		# Destroy references to deleted sections
		self._soma_refs = [ref for ref in self._soma_refs if ref.exists()]
		self._dend_refs = [ref for ref in self._dend_refs if ref.exists()]
		
		# Add newly created sections
		if soma_refs is not None:
			self._soma_refs = list(set(self._soma_refs + soma_refs)) # get unique references

		if dend_refs is not None:
			self._dend_refs = list(set(self._dend_refs + dend_refs)) # get unique references


	def set_reduction_params(self, method, params):
		"""
		Set parameters for given reduction method.
		"""
		self._REDUCTION_PARAMS[method] = params


	def get_reduction_param(self, method, param):
		"""
		Get reduction parameter for given method.
		"""
		return self._REDUCTION_PARAMS[method][param]


	def destroy(self):
		"""
		Release references to all stored data
		"""
		# Parameters for reduction method (set by user)
		self._REDUCTION_PARAMS = None
		self._mechs_gbars_dict = None

		self._soma_refs = None
		self._dend_refs = None
		self._root_ref = None
		self._fold_root_refs = None

		self._syns_tomap = None
		self._map_syn_info = None


	def _exec_reduction_step(self, step, method, step_args=None):
		"""
		Execute reduction step 'step' using method 'method'
		"""
		try:
			func, arg_names = self._REDUCTION_STEPS_FUNCS[step][method]
		
		except KeyError:
			raise NotImplementedError("{} function not implemented for "
									  "reduction method {}".format(step, method))
		
		else:
			user_params = self._REDUCTION_PARAMS[method]
			user_kwargs = dict((kv for kv in user_params.iteritems() if kv[0] in arg_names)) # get required args
			
			if step_args is None:
				step_args = []

			logger.anal("Executing reduction step {}".format(step))
			
			func(*step_args, **user_kwargs)


	def preprocess_cell(self, method):
		"""
		Pre-process cell: calculate properties & prepare data structures
		for reduction procedure

		@param	method		ReductionMethod instance: the reduction method that we
								should preprocess for.

		@pre		The somatic and dendritic sections have been set

		@post		Computed properties will be available as attributes
					on Section references in _soma_refs and _dend_refs,
					in addition to other side effects specified by the
					specific preprocessing function called.
		"""
		# Execute custom preprocess function
		self._exec_reduction_step('preprocess', method, step_args=[self])

		# Compute synapse mapping info
		if any(self._syns_tomap):

			# Existing synapse attributes to save (SectionRef attributes)
			save_ref_attrs = ['table_index', 'tree_index', 'gid']

			# Mapping parameters
			Z_freq			= self.get_reduction_param(method, 'Z_freq')
			init_func		= self.get_reduction_param(method, 'Z_init_func')
			linearize_gating= self.get_reduction_param(method, 'Z_linearize_gating')
			
			# Compute mapping data
			syn_info = mapsyn.get_syn_info(self.soma_refs[0].sec, self.all_sec_refs, 
								syn_tomap=self._syns_tomap, Z_freq=Z_freq, 
								init_cell=init_func, linearize_gating=linearize_gating,
								save_ref_attrs=save_ref_attrs)

			self._map_syn_info = syn_info


	def prepare_folds(self, method):
		"""
		Prepare next fold operation.
		"""
		self._exec_reduction_step('prepare_fold', method, step_args=[self])


	def calc_folds(self, method, i_pass):
		"""
		Fold branches at branch points identified by given criterion.
		"""
		self._exec_reduction_step('calculate_fold', method, step_args=[self, i_pass])


	def make_fold_equivalents(self, method):
		"""
		Make equivalent Sections for branches that have been folded.
		"""
		self._exec_reduction_step('make_fold', method, step_args=[self])


	def postprocess_cell(self, method=None):
		"""
		Post-process cell: perform any additional modifications for reduction
		algorithm.
		"""
		if method is None:
			method = self.active_method

		# Execute custom postprocess function
		self._exec_reduction_step('postprocess', method, step_args=[self])
		

	def map_synapses(self, method=None):
		"""
		Map any synapses if present

		@see	set_syns_tomap() for setting synapses.

		@pre	Any synapses provided through syns_tomap must be preprocessed
				for mapping (electronic properties computed), with results 
				stored in map_syn_info attribute.
		"""
		if method is None:
			method = self.active_method

		# Map synapses to reduced cell
		if any(self._map_syn_info):
			logger.debug("Mapping synapses...")

			# Mapping parameters
			Z_freq			= self.get_reduction_param(method, 'Z_freq')
			init_func		= self.get_reduction_param(method, 'Z_init_func')
			linearize		= self.get_reduction_param(method, 'Z_linearize_gating')
			map_method		= self.get_reduction_param(method, 'syn_map_method')
			
			# Map synapses (this modifies syn_info objects)
			mapsyn.map_synapses(self.soma_refs[0], self.all_sec_refs, self._map_syn_info, 
								init_func, Z_freq, linearize_gating=linearize,
								method=map_method)
		else:
			logger.debug("No synapse data available for mapping.")



	def reduce_model(self, num_passes, method=None, map_synapses=True):
		"""
		Do a fold-based reduction of the compartmental cell model.

		@param	num_passes		number of 'folding' passes to be done. One pass corresponds to
								folding at a particular node level (usually the highest).
		"""
		if method is None:
			method = self.active_method
		self.active_method = method # indicate what method we are using

		# Start reduction process
		self.preprocess_cell(method)

		# Fold one pass at a time
		for i_pass in xrange(num_passes):
			self.prepare_folds(method)
			self.calc_folds(method, i_pass)
			self.make_fold_equivalents(method)
			logger.debug('Finished folding pass {}'.format(i_pass))

		# Finalize reduction process
		self.postprocess_cell(method)

		# Map synapses
		if map_synapses:
			self.map_synapses(method=method)

		logger.debug('Finished cell reduction with method {}'.format(method))


################################################################################
# Cell model-specific tweaks
################################################################################

def adjust_gbar_spontaneous(reduction):
	"""
	Adjust gbar (NaL) to fix spontaneous firing rate.
	"""
	# Apply correction (TODO: remove this after fixing reduction)
	for ref in reduction.all_sec_refs:
		sec = ref.sec
		
		if sec.name().endswith('soma'):
			print("Skipping soma")
			continue
		
		logger.anal("Scaled gna_NaL in section {}".format(sec))
		for seg in sec:
			# seg.gna_NaL = 1.075 * seg.gna_NaL
			seg.gna_NaL = 8e-6 * 1.3 # full model value = uniform 8e-6


################################################################################
# Reduction Experiments
################################################################################

def gillies_stratford_reduction():
	"""
	Make FoldReduction object with Stratford method.
	"""
	from bgcellmodels.models.STN.GilliesWillshaw import gillies_model
	if not hasattr(h, 'SThcell'):
		gillies_model.stn_cell_gillies()

	# Make sections accesible by name and index
	soma = h.SThcell[0].soma
	dendL = list(h.SThcell[0].dend0) # 0 is left tree
	dendR = list(h.SThcell[0].dend1) # 1 is right tree
	dends = dendL + dendR

	# Get references to root sections of the 3 identical trees
	dendR_root			= h.SThcell[0].dend1[0]
	dendL_juction		= h.SThcell[0].dend0[0]
	dendL_upper_root	= h.SThcell[0].dend0[1] # root section of upper left dendrite
	dendL_lower_root	= h.SThcell[0].dend0[2] # root section of lower left dendrite
	fold_roots = [dendR_root, dendL_upper_root, dendL_lower_root]

	# Parameters for reduction
	def stn_setstate():
		# Initialize cell for analyzing electrical properties
		h.celsius = 35
		h.v_init = -68.0
		h.set_aCSF(4)
		h.init()

	# Reduce model
	red_method = ReductionMethod.Stratford
	reduction = FoldReduction([soma], dends, fold_roots, red_method)
	
	# Reduction parameters
	reduction.gleak_name = gillies_model.gleak_name
	reduction.mechs_gbars_dict = gillies_model.gillies_gdict
	reduction.set_reduction_params(red_method, {
		'fold_dX' :				0.25,
		'Z_freq' :				25.,
		'Z_init_func' :			stn_setstate,
		'Z_linearize_gating' :	False,
		'gbar_scaling' :		None,
		'syn_map_method' :		'Ztransfer',
	})

	return reduction


def fold_gillies_stratford(export_locals=True):
	"""
	Fold Gillies STN model using given reduction method
	
	@param	export_locals		if True, local variables will be exported to the global
								namespace for easy inspection
	"""
	# Make reduction object
	reduction = gillies_stratford_reduction()

	# Do reduction
	reduction.reduce_model(num_passes=1)

	if export_locals:
		globals().update(locals())

	return reduction._soma_refs, reduction._dend_refs


def gillies_marasco_reduction(tweak=True):
	"""
	Make FoldReduction object with Marasco method.
	"""

	from bgcellmodels.models.STN.GilliesWillshaw import gillies_model
	if not hasattr(h, 'SThcell'):
		gillies_model.stn_cell_gillies()

	# Make sections accesible by name and index
	soma = h.SThcell[0].soma
	dendL = list(h.SThcell[0].dend0) # 0 is left tree
	dendR = list(h.SThcell[0].dend1) # 1 is right tree
	dends = dendL + dendR

	# Get references to root sections of the 3 identical trees
	dendR_root			= h.SThcell[0].dend1[0]
	dendL_upper_root	= h.SThcell[0].dend0[1] # root section of upper left dendrite
	dendL_lower_root	= h.SThcell[0].dend0[2] # root section of lower left dendrite
	fold_roots = [dendR_root, dendL_upper_root, dendL_lower_root]

	# Parameters for reduction
	def stn_setstate():
		# Initialize cell for analyzing electrical properties
		h.celsius = 35
		h.v_init = -68.0
		h.set_aCSF(4)
		h.init()

	# Reduce model
	red_method = ReductionMethod.Marasco
	reduction = FoldReduction([soma], dends, fold_roots, red_method)

	# Reduction parameters
	reduction.gleak_name = gillies_model.gleak_name
	reduction.mechs_gbars_dict = gillies_model.gillies_gdict
	reduction.set_reduction_params(red_method, {
		'Z_freq' :				25.,
		'Z_init_func' :			stn_setstate,
		'Z_linearize_gating' :	False,
		'gbar_scaling' :		'area',
		'syn_map_method' :		'Ztransfer',
		'post_tweak_funcs' :	[adjust_gbar_spontaneous] if tweak else [],
	})

	return reduction


def fold_gillies_marasco(export_locals=True):
	"""
	Fold Gillies STN model using given reduction method
	
	@param	export_locals		if True, local variables will be exported to the global
								namespace for easy inspection
	"""
	# Make reduction object
	reduction = gillies_marasco_reduction()
	
	# Do reduction
	reduction.reduce_model(num_passes=7)
	reduction.update_refs()

	if export_locals:
		globals().update(locals())

	return reduction._soma_refs, reduction._dend_refs


if __name__ == '__main__':
	fold_gillies_marasco()