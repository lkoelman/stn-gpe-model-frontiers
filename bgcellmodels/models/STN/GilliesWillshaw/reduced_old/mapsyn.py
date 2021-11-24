"""
Algorithms for mapping synaptic input locations from detailed compartmental neuron model
to reduced morphology model.
"""

import re
from textwrap import dedent
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
logger = logging.getLogger('redops') # create logger for this module

from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

from bgcellmodels.models.STN.GilliesWillshaw import gillies_model as gillies
import redutils as redtools
from bgcellmodels.common import configutil
from redutils import getsecref, seg_index

class SynInfo(object):
	"""
	Struct-like object containing synapse properties.
	

	It may contain following attributes:

		mod_name			str
		
		sec_name			str
		
		sec_hname			str
		
		sec_loc				float

		'mech_attr_i'		float
							one attribute for each synaptic mechanism parameter

		afferent_netcons	list(NetCon)
		
		afferent_weights	list(list(float))
							weight vector for each incoming NetCon

		'secref_attr_i'		object
							saved SectionRef attributes


		path_ri				float
							summed seg.ri() up to synapse segment
		
		max_path_ri			float
							= max(syn_sec.pathri_seg) # max path resistance in Section
		
		min_path_ri			float
							= min(syn_sec.pathri_seg) # min path resistance in Section

		Zc					float
							Ztransfer, i.e. |v(soma)/i(syn)| or |v(syn)/i(soma)|
		
		Zin					float
							Zinput, i.e. v(x)/i(x) measured at synapse

		k_syn_soma			float
							voltage transfer ratio, i.e. |v(soma)/v(syn|

		mapped_syn			HocObject
							the new synapse that the original was mapped to
	"""
	def __init__(self, **kwds):
		self.__dict__.update(kwds)

	# def __repr__(self):
	# 	return "{} at {}({})".format(self.mod_name, self.sec_hname, self.sec_loc)


def subtree_has_node(crit_func, noderef, allsecrefs):
	"""
	Check if the given function applies to (returns True) any node
	in subtree of given node
	"""
	if noderef is None:
		return False
	elif crit_func(noderef):
		return True
	else:
		childsecs = noderef.sec.children()
		childrefs = [getsecref(sec, allsecrefs) for sec in childsecs]
		for childref in childrefs:
			if subtree_has_node(crit_func, childref, allsecrefs):
				return True
		return False

# Parameter names for synaptic mechanisms in .MOD files
synmech_parnames = {
	'ExpSyn': ['tau', 'e'], # gmax stored in NetCon weight vector
	
	'Exp2Syn': ['tau1', 'tau2', 'e'], # gmax stored in NetCon weight vector
	
	'AlphaSynapse': ['onset', 'tau', 'gmax', 'e'],
	
	'GABAsyn': ['tau_r_GABAA', 'tau_d_GABAA', 'tau_r_GABAB', 'tau_d_GABAB', 
				'gmax_GABAA', 'gmax_GABAB', 'Erev_GABAA', 'Erev_GABAB', 
				'tau_rec', 'tau_facil', 'U1', 'use_stdp_A', 'use_stdp_B'],
	
	'GLUsyn': ['tau_r_AMPA', 'tau_d_AMPA', 'tau_r_NMDA', 'tau_d_NMDA', 'mg', 
				'gmax_AMPA', 'gmax_NMDA', 'tau_rec', 'tau_facil', 'U1', 'e'],
}

def get_syn_info(rootsec, allsecrefs, syn_mod_pars=None, Z_freq=25., linearize_gating=False,
					init_cell=None, save_ref_attrs=None, sever_netcons=True,
					attr_mappers=None, nc_tomap=None, syn_tomap=None):
	"""
	For each synapse on the neuron, calculate and save information for placing an equivalent
	synaptic input on a morphologically simplified neuron.

	
	ARGUMENTS

		@param rootsec			any section of the cell

		@param syn_mod_pars		dict of <synaptic mechanism name> : list(<attribute names>)
								containing attributes that need to be stored

		@param init_cell		function to bring the cell to the desired state to measure transfer
								impedances, e.g. simulating under a particular input

		@param sever_netcons	Set NetCon target to None for each connection, this is useful
								when the cell or synapses are destroyed


	RETURN VALUES

		@return					list(SynInfo) containing information about each synapse
								to be mapped


	USAGE

		- either provide syn_tomap (list(SynInfo)), or nc_tomap, or neither.
		  In the latter case, all NetCon targetting the cell are used.
	"""
	if syn_mod_pars is None:
		syn_mod_pars = synmech_parnames

	if save_ref_attrs is None:
		save_ref_attrs = []

	if attr_mappers is None:
		attr_mappers = {}

	# Calculate section path properties for entire tree
	for secref in allsecrefs:
		redtools.sec_path_props(secref, 100., gillies.gleak_name)

	# Measure transfer impedance and filter parameters
	logger.debug("Initializing cell for electrotonic measurements...")
	init_cell()
	
	imp = h.Impedance() # imp = h.zz # previously created
	imp.loc(0.5, sec=rootsec) # injection site

	# NOTE: linearize_gating corresponds to checkbox 'include dstate/dt contribution' in NEURON GUI Impedance Tool
	#		- (see equations in Impedance doc) 
	#		- 0 = calculation with current values of gating vars
	#		- 1 = linearize gating vars around V
	# imp.compute(Z_freq, int(linearize_gating)) # compute A(x->loc) for all x where A is Vratio/Zin/Ztransfer
	last_freq = None


	# Find all Synapses on cell (all Sections in whole tree)
	if syn_tomap is not None:
		cell_syns = syn_tomap

	else:
		# Get references to all synapses ourselves
		if nc_tomap is None:
			# Get all NetCon targetting the cell
			dummy_syn = h.Exp2Syn(rootsec(0.5))
			dummy_nc = h.NetCon(None, dummy_syn)

			# unique synapses targeting the same cell
			cell_ncs = [nc for nc in list(dummy_nc.postcelllist()) if not nc.same(dummy_nc)] # all NetCon targeting same tree as dummy
			syns = set([nc.syn() for nc in cell_ncs]) # NOTE: this works, since set uses == for uniqueness, which compares using syn.hocobjptr()
		
		else:
			cell_ncs = nc_tomap
			syns = set([nc.syn() for nc in cell_ncs]) # NOTE: this works, since set uses == for comparison, which uses syn.hocobjptr()

		# Make a list of SynInfo ourselves
		cell_syns = [SynInfo(orig_syn=syn) for syn in syns]
		for syn in cell_syns:
			syn.PSP_median_frequency = Z_freq
			syn.afferent_netcons = [nc for nc in cell_ncs if nc.syn().same(syn.orig_syn)]


	if len(cell_syns) == 0:
		logger.warn("No synapses found on tree of Section {}".format(rootsec))


	# Save synapse properties
	logger.debug("Getting synapse properties...")
	for syn_info in sorted(cell_syns, key=lambda syn: syn.PSP_median_frequency):

		# get synapse HocObject
		syn = syn_info.orig_syn

		# Compute transfer function at representative frequency
		syn_freq = syn_info.PSP_median_frequency
		if (last_freq is None) or (syn_freq != last_freq):
			imp.compute(syn_freq, int(linearize_gating))
		
		# Get synaptic mechanism name
		match_mechname = re.search(r'^[a-zA-Z0-9]+', syn.hname())
		synmech = match_mechname.group()
		if synmech not in syn_mod_pars:
			raise Exception("Synaptic mechanism '{}' not in given mechanism list".format(synmech))

		# Get its segment and location
		syn_seg = syn.get_segment()
		syn_sec = syn_seg.sec
		syn_secref = getsecref(syn_sec, allsecrefs)
		syn_loc = syn.get_loc() # changes CAS
		h.pop_section()

		# Create struct to save synapse information
		syn_info.mod_name = synmech
		syn_info.sec_name = syn_sec.name()
		syn_info.sec_hname = syn_sec.hname()
		syn_info.sec_loc = syn_loc # can also use nc.postcell() and nc.postloc()

		# Save synaptic mechanism parameters
		mech_params = syn_mod_pars[synmech]
		for par in mech_params:
			setattr(syn_info, par, getattr(syn, par))

		# Save all NetCon objects targetting this synapse
		syn_info.afferent_weights = [[nc.weight[i] for i in range(int(nc.wcnt()))] for nc in syn_info.afferent_netcons]

		# Save requested properties of synapse SectionRef
		syn_info.saved_ref_attrs = save_ref_attrs
		for attr in save_ref_attrs:
			setattr(syn_info, attr, getattr(syn_secref, attr))

		# Save other computed properties
		for attr, mapper in attr_mappers.items():
			setattr(syn_info, attr, mapper(syn))

		# Get axial path resistance to synapse
		syn_info.path_ri = syn_secref.pathri_seg[seg_index(syn_seg)] # summed seg.ri() up to synapse segment
		syn_info.max_path_ri = max(syn_secref.pathri_seg) # max path resistance in Section
		syn_info.min_path_ri = min(syn_secref.pathri_seg) # min path resistance in Section

		# Get transfer impedances
		syn_info.Zc = imp.transfer(syn_loc, sec=syn_sec) # query transfer impedanc,e i.e.  |v(loc)/i(x)| or |v(x)/i(loc)|
		syn_info.Zin = imp.input(syn_loc, sec=syn_sec) # query input impedance, i.e. v(x)/i(x)
		syn_info.k_syn_soma = imp.ratio(syn_loc, sec=syn_sec) # query voltage transfer ratio, i.e. |v(loc)/v(x)|

		# Point all afferent NetCon connections to nothing
		if sever_netcons:
			for nc in syn_info.afferent_netcons:
				nc.setpost(None)

		# Delete reference to synapse about to be deleted
		syn_info.orig_syn = None

	return cell_syns



def map_synapse(noderef, allsecrefs, syn_info, imp, method, passed_synsec=False):
	"""
	Map synapse to a section in subtree of noderef
	"""
	cur_sec = noderef.sec

	if method == 'Ztransfer':
		Zc = syn_info.Zc
		elec_fun = imp.transfer
	
	elif method == 'Vratio':
		Zc = syn_info.k_syn_soma
		elec_fun = imp.ratio
	
	else:
		raise Exception("Unknown synapse placement method '{}'".format(method))

	Zc_0 = elec_fun(0.0, sec=cur_sec)
	Zc_1 = elec_fun(1.0, sec=cur_sec)

	logger.anal("Entering section with Zc(0.0)={} , Zc(1.0)={} (target Zc={})".format(Zc_0, Zc_1, Zc))

	# Assume monotonically decreasing Zc away from root section
	if (Zc_0 <= Zc <= Zc_1) or (Zc_1 <= Zc <= Zc_0) or (Zc >= Zc_0 and Zc >= Zc_1):

		# Calculate Zc at midpoint of each internal segment
		locs_Zc = [(seg.x, elec_fun(seg.x, sec=cur_sec)) for seg in cur_sec]
		Zc_diffs = [abs(Zc-pts[1]) for pts in locs_Zc]

		# Map synapse with closest Zc at midpoint
		seg_index = Zc_diffs.index(min(Zc_diffs))
		x_map = locs_Zc[seg_index][0]

		err_Zc = (locs_Zc[seg_index][1] - Zc) / Zc
		logger.debug("Arrived: Zc relative error is {:.2f} %".format(err_Zc))
		
		return noderef, x_map	

	else:
		logger.anal("No map: descending further...")
		
		# Get child branches
		childrefs = [getsecref(sec, allsecrefs) for sec in noderef.sec.children()]

		# If we are in correct tree but Zc smaller than endpoints, return endpoint
		if not any(childrefs):
			assert Zc < Zc_0 and Zc < Zc_1
			return noderef, 1.0

		# Else, recursively search child nodes
		# Function for checking that section maps to original synapse section
		mapsto_synsec = lambda ref: (ref.tree_index==syn_info.tree_index) and (
									ref.gid==syn_info.gid or (
										hasattr(ref, 'zipped_sec_gids') and 
										(syn_info.gid in ref.zipped_sec_gids)
										)
									)
		
		# Did we pass the synapse's original section?
		passed_synsec = passed_synsec or mapsto_synsec(noderef)
		
		# child_mapsto_synsec = [subtree_has_node(mapsto_synsec, ref, allsecrefs) for ref in childrefs]
		# passed_synsec = not any(child_mapsto_synsec)
		# if passed_synsec:
		# 	assert mapsto_synsec(noderef) # assume that synapse was positioned on this Section

		for childref in childrefs:

			# Only descend if original synapse section already passed, or in subtree
			if passed_synsec or subtree_has_node(mapsto_synsec, childref, allsecrefs):
				return map_synapse(childref, allsecrefs, syn_info, imp, method, passed_synsec)
		
		raise Exception("The synapse did not map onto any segment in this subtree.")



def map_synapses(rootref, allsecrefs, orig_syn_info, init_cell, Z_freq, 
					syn_mod_pars=None, method='Ztransfer', linearize_gating=False):
	"""
	Map synapses to equivalent synaptic inputs on given morphologically
	reduced cell.

	@param rootsec			root section of reduced cell

	@param orig_syn_info	SynInfo for each synapse on original cell

	@param method			Method for positioning synapses and scaling
							synaptic conductances:

							'Ztransfer' positions synapses at loc with 
										~= transfer impedance

							'Vratio' positions them at loc with 
										~= Voltage attenuation

	@effect					Create one synapse for each original synapse in 
							orig_syn_info. A reference to this synapse is saved 
							as as attribute 'mapped_syn' on the SynInfo object.
	"""
	# Synaptic mechanisms
	if syn_mod_pars is None:
		syn_mod_pars = synmech_parnames

	# Compute transfer impedances
	logger.debug("Initializing cell for electrotonic measurements...")
	init_cell()
	logger.debug("Placing impedance measuring electrode...")
	
	imp = h.Impedance() # imp = h.zz # previously created
	imp.loc(0.5, sec=rootref.sec) # injection site

	# NOTE: linearize_gating corresponds to checkbox 'include dstate/dt contribution' in NEURON GUI Impedance Tool
	#		- (see equations in Impedance doc) 
	#		- 0 = calculation with current values of gating vars
	#		- 1 = linearize gating vars around V
	# imp.compute(Z_freq, int(linearize_gating)) # compute A(x->loc) for all x where A is Vratio/Zin/Ztransfer
	last_freq = None
	for syn in orig_syn_info:
		if not hasattr(syn, 'PSP_median_frequency'):
			syn.PSP_median_frequency = Z_freq

	# Loop over all synapses
	logger.debug("Mapping synapses to reduced cell...")
	for syn_info in sorted(orig_syn_info, key=lambda syn: syn.PSP_median_frequency):

		# Compute transfer function at representative frequency
		syn_freq = syn_info.PSP_median_frequency
		if (last_freq is None) or (syn_freq != last_freq):
			imp.compute(syn_freq, int(linearize_gating))


		# Find the segment with same tree index and closest Ztransfer match,
		map_ref, map_x = map_synapse(rootref, allsecrefs, syn_info, imp, method)
		map_sec = map_ref.sec
		logger.anal("Synapse was in {}({}) -> mapped to {}\n".format(
						syn_info.sec_name, syn_info.sec_loc, map_sec(map_x)))

		# Make the synapse
		syn_mod = syn_info.mod_name
		synmech_ctor = getattr(h, syn_mod) # constructor function for synaptic mechanism
		mapped_syn = synmech_ctor(map_sec(map_x))
		syn_info.mapped_syn = mapped_syn

		# Copy synapse properties
		mech_params = syn_mod_pars[syn_mod]
		for par in mech_params:
			val = getattr(syn_info, par)
			setattr(mapped_syn, par, val)

		# Change target of afferent connections
		for aff_nc in syn_info.afferent_netcons:
			aff_nc.setpost(mapped_syn)
			aff_nc.active(1) # NetCons are turned off if target was set to None!

		# Measure electrotonic properties for scaling synapses
		map_Zc = imp.transfer(map_x, sec=map_sec) # query transfer impedanc,e i.e.  |v(loc)/i(x)| or |v(x)/i(loc)|
		map_Zin = imp.input(map_x, sec=map_sec) # query input impedance, i.e. v(x)/i(x)
		map_k = imp.ratio(map_x, sec=map_sec) # query voltage transfer ratio, i.e. |v(loc)/v(x)|

		# Ensure synaptic conductances scaled correctly to produce same effect at soma
		# using relationship `Vsoma = Zc*gsyn*(V-Esyn) = k_{syn->soma}*Zin*gsyn*(V-Esyn)`

		# Report discrepancy
		logger.anal(dedent("""
			Original synapse:\nZc={}\tk={}\tZin={}
			Mapped synapse:\nZc={}\tk={}\tZin={}
			Discrepancy: Zc_old/Zc_new={} = scale factor for gmax
		""".format(syn_info.Zc, syn_info.k_syn_soma, syn_info.Zin,
					map_Zc, map_k, map_Zin, syn_info.Zc/map_Zc)))

		# Calculate scale factor
		#	`Vsoma = k_{syn->soma} * Vsyn`			(Vsyn is local EPSP)
		#	`Vsoma = k_{syn->soma} * Zin * Isyn`	(Isyn is synaptic current source)
		#	`Vsoma = (Zc / Zin) * Zin * Isyn`
		#	`Vsoma = Zc * Isyn`
		#	`Vsoma = Zc * gsyn * (V-Esyn)`
		if method == 'Ztransfer':
			# method 1: (don't conserve local Vsyn):
			#			- place at x with +/- same Zc 
			#			- ensure Isyn is the same
			#				- correct gsyn using exact Zc measurement (see below)
			scale_g = syn_info.Zc/map_Zc
		
		elif method == 'Vratio':
			# Method 2: (conserve local Vsyn):
			#			- place at x with +/- same k_{syn->soma}
			#			- ensure that Vsyn is the same (EPSP, local depolarization)
			#				- to maintain Vsyn = Zin * gsyn * (V-Esyn) : scale gsyn
			#				  to compensate for discrepancy in Zin
			scale_g = syn_info.Zin/map_Zin
		
		else:
			raise Exception("Unknown synapse placement method '{}'".format(method))

		# Scale conductances (weights) of all incoming connections
		scale_netcon_weights = False

		if len(getattr(syn_info, 'gbar_param_specs', [])) > 0:
			for gbar_spec in syn_info.gbar_param_specs:

				mechtype, mech_parname, mech_paridx = configutil.interpretParamSpec(gbar_spec)
				
				if mechtype == 'syn':
					target = mapped_syn
					if mech_paridx is None:
						old_val = getattr(target, mech_parname)
						setattr(target, mech_parname, old_val*scale_g)
					else:
						old_val = getattr(target, mech_parname)[int(mech_paridx)]
						getattr(target, mech_parname)[int(mech_paridx)] = old_val*scale_g
				
				elif mechtype == 'netcon':
					scale_netcon_weights = True

		else:
			scale_netcon_weights = True
					

		if scale_netcon_weights:
			# Default action: scale NetCon weights
			for i_nc, nc in enumerate(syn_info.afferent_netcons):
				
				# Get original weights (weights are reset when disconnecting)
				orig_weights = syn_info.afferent_weights[i_nc]
				assert len(orig_weights) == int(nc.wcnt())
				
				# Scale original weights
				for i_w in range(int(nc.wcnt())):
					nc.weight[i_w] = orig_weights[i_w] * scale_g
					logger.anal("Scaled weight {} by factor {}".format(orig_weights[i_w], scale_g))

			
		

		
