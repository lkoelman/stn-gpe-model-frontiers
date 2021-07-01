# Gillies & Willshaw model mechanisms
import math
from neuron import h

# Custom modules
from bgcellmodels.models.STN.GilliesWillshaw import gillies_model
from bgcellmodels.models.STN.GilliesWillshaw.reduced_old import reduce_cell

from bgcellmodels.cersei.collapse.optimize.cellmodels import CollapsableCell # Cell as MockCell
from bgcellmodels.models.STN.GilliesWillshaw.reduced import cersei_reduce

# Adjust verbosity of loggers
import logging
logging.getLogger('redops').setLevel(logging.WARNING)
logging.getLogger('folding').setLevel(logging.WARNING)
logging.getLogger('marasco').setLevel(logging.WARNING)

# Module globals
gleak_name = gillies_model.gleak_name


def optimize_passive(export_locals=True):
	"""
	Optimize reduced cell model to have same passive response.

	Algorithm

		- make cell passive by setting all active conductances to zero
		
		- use Brent's PRAXIS optimization method to match Rin measured at soma
		
		- scale sec.cm to preserve tau = 1/sec.gleak * sec.cm
	

	Notes
	
		- explanation PRAXIS optimization method: http://reference.wolfram.com/language/tutorial/UnconstrainedOptimizationPrincipalAxisMethod.html
		
		- NEURON PRAXIS optimization: see https://www.neuron.yale.edu/neuron/static/new_doc/analysis/programmatic/optimization.html
		
		- for example usage of fit_praxis(): search for "fit_praxis" on ModelDB
	"""

	# Create full model
	somaref, dendLrefs, dendRrefs = gillies_model.get_stn_refs()
	all_refs = [somaref] + dendLrefs + dendRrefs
	gillies_model.stn_init_physiology()
	gillies_model.make_passive(all_refs)

	# Measure Zin_DC
	imp_full = h.Impedance()
	linearize_gating = False
	imp_full.loc(0.5, sec=h.SThcell[0].soma) # injection site
	imp_full.compute(0.0, int(linearize_gating))
	Zin_target = imp_full.input(0.5, sec=h.SThcell[0].soma)
	del imp_full
	print("Zin_DC in full passive model is {}".format(Zin_target))

	# Create reduced model
	reduction = reduce_cell.gillies_marasco_reduction()
	reduction.reduce_model(num_passes=7)
    
    # Make reduction object
    reduction = make_reduction(method=ReductionMethod.BushSejnowski)
    reduction.reduce_model(num_passes=7)
    
    # Disable active conductances
	gillies_model.make_passive(reduction.all_sec_refs)

	# Optimization cost function
	imp_red = h.Impedance()
	imp_red.loc(0.5, sec=h.SThcell[0].soma) # injection site

	def cost_fun(param_hvec):
		"""
		Cost function (error measure) to minimize
		"""
		# Adjust leak conducance
		gleak_val = param_hvec.x[0]
		for sec in h.allsec():
			if 'soma' in sec.name():
				print("Skipping soma")
				continue
			for seg in sec:
				setattr(seg, gleak_name, gleak_val)

		# Measure input inpedance
		imp_red.compute(0.0, int(linearize_gating))
		Zin_red = imp_red.input(0.5, sec=h.SThcell[0].soma)
		return (Zin_red - Zin_target)**2

	# Parameters for optimization
	Zin_tolerance = 1e-6 # MOhm (full model: 81.88 MOhm)
	gleak_step = .1e-5 # original gleak = 7.84112e-05
	praxis_verbosity = 2
	h.attr_praxis(Zin_tolerance, gleak_step, praxis_verbosity)
	
	# Optimize
	gleak_orig = 7.84112e-05
	param_vec = h.Vector([gleak_orig])
	Zin_err_sqrd = h.fit_praxis(cost_fun, param_vec)

	# Report results
	Zin_final = imp_red.input(0.5, sec=h.SThcell[0].soma)
	Zin_error = math.sqrt(Zin_err_sqrd)
	gleak_fit = param_vec.x[0]
	g_ratio = gleak_fit / gleak_orig # use to scale cm to tau_m is conserved

	print("""Optimization finished:
		gleak_orig	= {}
		gleak_fit	= {}
		gleak_ratio	= {}
		Zin_orig	= {}
		Zin_fit		= {}
		Zin error	= {}""".format(gleak_orig, gleak_fit, g_ratio, Zin_target, Zin_final, Zin_error))

	if export_locals:
		globals().update(locals())

	return gleak_fit, g_ratio

if __name__ == '__main__':
	optimize_passive(export_locals=True)