"""
Reduce Gillies & Willshaw (2006) STN neuron model using the method
described in Marasco & Migliore (2012)


@author Lucas Koelman
@date	5-12-2016
"""

# Python modules
import re
import math
PI = math.pi

# NEURON modules
import neuron
h = neuron.h
h.load_file("stdlib.hoc") # Load the standard library

# Own modules
import redutils as redtools
from bgcellmodels.common.nrnutil import ExtSecRef, getsecref, seg_index # for convenience
import cluster as clutools
from cluster import Cluster
import interpolation as interp
from marasco_merging import merge_seg_subtree
import folding

# Gillies STN model
from bgcellmodels.models.STN.GilliesWillshaw.gillies_model import (
	gillies_gdict, gillies_glist)

mechs_chans = gillies_gdict
glist = gillies_glist
gleak_name = 'gpas_STh'
f_lambda = 100.0

# logging of DEBUG/INFO/WARNING messages
import logging
logging.basicConfig(format='%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logname = "marasco" # __name__
logger = logging.getLogger(logname) # create logger for this module

# Log to file
# fmtr = logging.Formatter('%(levelname)s:%(message)s @%(filename)s:%(lineno)s')
# fh = logging.FileHandler('reduce_marasco.log')
# fh.setFormatter(fmtr)
# logger.addHandler(fh)
# Log to stream
# ch = logging.StreamHandler(sys.stdout)
# ch.setFormatter(fmtr)
# logger.addHandler(ch)



################################################################################
# Merging / Folding segments
################################################################################


def calc_collapses(target_Y_secs, i_pass, allsecrefs):
	"""
	Do collapse operations: calculate equivalent Section properties for each collapse.

	@return			list of Cluster objects with properties of equivalent
					Section for each set of collapsed branches.
	"""

	# Flag all sections as unvisited
	for secref in allsecrefs:

		if not hasattr(secref, 'absorbed'):
			secref.absorbed = [False] * secref.sec.nseg
			secref.visited = [False] * secref.sec.nseg
		
		secref.zip_labels = [None] * secref.sec.nseg

	# Create Clusters: collapse (zip) each Y section up to length of first branch point
	clusters = []
	for j_zip, par_ref in enumerate(target_Y_secs):
		par_sec = par_ref.sec
		child_secs = par_sec.children()

		# Function to determine which segment can be 'zipped'
		min_child_L = min(sec.L for sec in child_secs) # Section is unbranched cable
		eligfunc = lambda seg, jseg, ref: (ref.parent.same(par_sec)) and (seg.x*seg.sec.L <= min_child_L)
		
		# Name for equivalent zipped section
		# name_sanitized = par_sec.name().replace('[','').replace(']','').replace('.','_')
		name_sanitized = re.sub(r"[\[\]\.]", "", par_sec.name())
		alphabet_uppercase = [chr(i) for i in range(65,90+1)] # A-Z are ASCII 65-90
		zip_label = "zip{0}_{1}".format(alphabet_uppercase[i_pass], name_sanitized)
		zip_id = 1000*i_pass + j_zip

		# Function for processing zipped SectionRefs
		zipped_sec_gids = set()
		def process_zipped_seg(seg, jseg, ref):
			# Tag segment with label of current zip operation
			ref.zip_labels[jseg] = zip_label
			# Save GIDs of original cell that are zipped/absorbed into equivalent section
			if hasattr(ref, 'table_index') and ref.table_index >= 1:
				zipped_sec_gids.add(ref.gid)
			else:
				zipped_sec_gids.update(ref.zipped_sec_gids)
		
		# Perform 'zip' operation
		far_bound_segs = [] # last (furthest) segments that are zipped
		eq_seq, eq_br = merge_seg_subtree(par_sec(1.0), allsecrefs, eligfunc, 
							process_zipped_seg, far_bound_segs)

		logger.debug("Target Y-section for zipping: {0} (tree_id:{1} , table_id:{2})".format(
						par_sec.name(), par_ref.tree_index, par_ref.table_index))
		bounds_info = ["\n\tsegment {0} [{1}/{2}]".format(seg, seg_index(seg)+1, seg.sec.nseg) for seg in far_bound_segs]
		logger.debug("Zipping up to {0} boundary segments:{1}".format(len(far_bound_segs), "\n".join(bounds_info)))

		# Make Cluster object that represents collapsed segments
		cluster = Cluster(zip_label)
		cluster.eqL = eq_br.L_eq
		cluster.eqdiam = eq_br.diam_eq
		cluster.eq_area_sum = eq_br.L_eq * PI * eq_br.diam_eq
		cluster.eqri = eq_br.Ri_eq
		cluster.zipped_sec_gids = zipped_sec_gids
		cluster.zip_id = zip_id
		cluster.parent_seg = par_sec(1.0)
		cluster.bound_segs = far_bound_segs # Save boundaries (for substitution)

		# Calculate cluster statistics #########################################
		
		# Gather all cluster sections & segments
		clu_secs = [secref for secref in allsecrefs if (cluster.label in secref.zip_labels)]
		clu_segs = [seg for ref in clu_secs for jseg,seg in enumerate(ref.sec) if (
						ref.zip_labels[jseg]==cluster.label)]

		# Calculate max/min path length
		clu_path_L = [ref.pathL_seg[j] for ref in clu_secs for j,seg in enumerate(ref.sec) if (
						ref.zip_labels[j]==cluster.label)]
		cluster.orMaxPathL = max(clu_path_L)
		cluster.orMinPathL = min(clu_path_L)

		# Calculate min/max axial path resistance
		clu_path_ri = [ref.pathri_seg[j] for ref in clu_secs for j,seg in enumerate(ref.sec) if (
						ref.zip_labels[j]==cluster.label)]
		cluster.orMaxpathri = max(clu_path_ri)
		cluster.orMinpathri = min(clu_path_ri)

		# Calculate area, capacitance, conductances
		cluster.or_area = sum(seg.area() for seg in clu_segs)
		cluster.or_cmtot = sum(seg.cm*seg.area() for seg in clu_segs)
		cluster.or_cm = cluster.or_cmtot / cluster.or_area
		cluster.or_gtot = dict((gname, 0.0) for gname in glist)
		for gname in glist:
			cluster.or_gtot[gname] += sum(getattr(seg, gname)*seg.area() for seg in clu_segs)

		# Equivalent axial resistance
		clu_segs_Ra = [seg.sec.Ra for seg in clu_segs]
		if min(clu_segs_Ra) == max(clu_segs_Ra):
			cluster.eqRa = clu_segs_Ra[0]
		else:
			logger.warning("Sections have non-uniform Ra, calculating average "
							"axial resistance per unit length, weighted by area")
			cluster.eqRa = PI*(cluster.eqdiam/2.)**2*cluster.eqri*100./cluster.eqL # eq. Ra^eq
		clusters.append(cluster)

		# Calculate electrotonic path length
		cluster.or_L_elec = sum(redtools.seg_L_elec(seg, gleak_name, f_lambda) for seg in clu_segs)
		cluster.eq_lambda = redtools.calc_lambda_AC(f_lambda, cluster.eqdiam, cluster.eqRa, cluster.or_cmtot/cluster.eq_area_sum)
		cluster.eq_L_elec = cluster.eqL/cluster.eq_lambda
		eq_min_nseg = redtools.calc_min_nseg_hines(f_lambda, cluster.eqL, cluster.eqdiam, 
											cluster.eqRa, cluster.or_cmtot/cluster.eq_area_sum)

		# Debug
		print_attrs = ['eqL', 'eqdiam', 'or_area', 'eq_area_sum', 'or_L_elec', 'eq_L_elec']
		clu_info = ("- {0}: {1}".format(prop, getattr(cluster, prop)) for prop in print_attrs)
		logger.debug("Equivalent section for zipped Y-section has following properties:\n\t{0}".format(
						"\n\t".join(clu_info)))
		logger.debug("Zip reduces L/lambda by {0:.2f} %; number of segments saved is {1} (Hines rule)\n".format(
						cluster.eq_L_elec/cluster.or_L_elec, len(clu_segs)-eq_min_nseg))

	return clusters


################################################################################
# Equivalent Section creation
################################################################################


def sub_equivalent_Y_sections(clusters, orsecrefs, interp_path, interp_prop='path_L', 
								interp_method='linear_neighbors', 
								gbar_scaling='area'):
	"""
	Substitute equivalent Section for each cluster into original cell.

	@see	docstring of function `reduce_bush_sejnowski.equivalent_sections()`

	@param	interp_prop		property used for calculation of path length (x of interpolated
							x,y values), one of following:
							
							'path_L': path length (in micrometers)

							'path_ri': axial path resistance (in Ohms)

							'path_L_elec': electrotonic path length (L/lambda, dimensionless)

	@param	interp_method	how numerical values are interpolated, one of following:	
							
							'linear_neighbors':
								linear interpolation of 'adjacent' segments in full model 
								(i.e. next lower and higher electrotonic path length). 
							
							'linear_dist':
								estimate linear distribution and interpolate it
							
							'left_neighbor', 'right_neighbor', 'nearest_neighbor':
								extrapolation of 'neighoring' segments in terms of L/lambda
	"""
	# List with equivalent section for each cluster
	eq_secs = []
	eq_refs = []

	# Create equivalent sections and passive electric structure
	for i, cluster in enumerate(clusters):

		# Create equivalent section
		if cluster.label in [sec.name() for sec in h.allsec()]:
			raise Exception('Section named {} already exists'.format(cluster.label))
		created = h("create %s" % cluster.label)
		if created != 1:
			raise Exception("Could not create section with name '{}'".format(cluster.label))
		eqsec = getattr(h, cluster.label)
		eqref = ExtSecRef(sec=eqsec)

		# Set geometry 
		eqsec.L = cluster.eqL
		eqsec.diam = cluster.eqdiam
		cluster.eq_area = sum(seg.area() for seg in eqsec) # should be same as cluster eqSurf

		# Passive electrical properties (except Rm/gleak)
		eqsec.Ra = cluster.eqRa

		# Save metadata
		eqref.zipped_sec_gids = cluster.zipped_sec_gids
		eqref.zip_id = cluster.zip_id

		# Append to list of equivalent sections
		eq_secs.append(eqsec)
		eq_refs.append(eqref)

		# Connect to tree (need to trace path from soma to section)
		eqsec.connect(cluster.parent_seg, 0.0) # see help(sec.connect)

	# Set active properties and finetune
	for i_clu, cluster in enumerate(clusters):
		logger.debug("Scaling properties of cluster %s ..." % clusters[i_clu].label)
		eqsec = eq_secs[i_clu]
		eqref = eq_refs[i_clu]

		# Insert all mechanisms
		for mech in mechs_chans.keys():
			eqsec.insert(mech)

		# Scale passive electrical properties
		area_ratio = cluster.or_area / cluster.eq_area
		logger.debug("Surface area ratio is %f" % area_ratio)

		# Scale Cm
		eq_cm = cluster.or_cmtot / cluster.eq_area # more accurate than cm * or_area/eq_area

		# Scale Rm
		or_gleak = cluster.or_gtot[gleak_name] / cluster.or_area
		eq_gleak = or_gleak * area_ratio # same as reducing Rm by area_new/area_old

		# Set number of segments based on rule of thumb electrotonic length
		eqsec.nseg = redtools.calc_min_nseg_hines(f_lambda, eqsec.L, eqsec.diam, eqsec.Ra, eq_cm)

		# Save Cm and conductances for each section for reconstruction
		cluster.nseg = eqsec.nseg # save for reconstruction
		cluster.eq_gbar = dict((gname, [float('NaN')]*cluster.nseg) for gname in glist)
		cluster.eq_cm = [float('NaN')]*cluster.nseg

		# Set Cm and gleak (Rm) for each segment
		if gbar_scaling is not None:
			for j, seg in enumerate(eqsec):
				setattr(seg, 'cm', eq_cm)
				setattr(seg, gleak_name, eq_gleak)
				cluster.eq_cm[j] = eq_cm
				cluster.eq_gbar[gleak_name][j] = eq_gleak

		# Get active conductances
		active_glist = list(glist)
		active_glist.remove(gleak_name) # get list of active conductances

		# Calculate path lengths in equivalent section
		redtools.sec_path_props(eqref, f_lambda, gleak_name)
		
		if interp_prop == 'path_L':
			seg_prop = 'pathL_seg'
		elif interp_prop == 'path_ri':
			seg_prop = 'pathri_seg'
		elif interp_prop == 'path_L_elec':
			seg_prop = 'pathL_elec'
		else:
			raise ValueError("Unknown path property '{}'".format(interp_prop))

		# Find conductances at same path length (to each segment midpoint) in original cell
		for j_seg, seg in enumerate(eqsec):
			
			# Get adjacent segments along path
			path_L = getattr(eqref, seg_prop)[j_seg]
			bound_segs, bound_L = interp.find_adj_path_segs(interp_prop, path_L, interp_path)
			# bounds_info = "\n".join(("\t- bounds {0} - {1}".format(a, b) for a,b in bound_segs))
			# logger.debug("Found boundary segments at same path length x={0:.3f}:\n{1}".format(path_L, bounds_info))

			# INTERPOLATE: Set conductances by interpolating neighbors
			for gname in active_glist:
				
				if interp_method == 'linear_neighbors':
					gval = interp.interp_gbar_linear_neighbors(path_L, gname, bound_segs, bound_L)
				
				else:
					match_method = re.search(r'^[a-z]+', interp_method)
					method = match_method.group() # should be nearest, left, or right
					gval = interp.interp_gbar_pick_neighbor(path_L, gname, 
										bound_segs[0], bound_L[0], method)
				seg.__setattr__(gname, gval)
				cluster.eq_gbar[gname][j_seg] = gval

		# Re-scale gbar distribution to yield same total gbar (sum(gbar*area))
		if gbar_scaling is not None:
			for gname in active_glist:
				
				eq_gtot = sum(getattr(seg, gname)*seg.area() for seg in eqsec)
				if eq_gtot <= 0.:
					eq_gtot = 1.
				
				or_gtot = cluster.or_gtot[gname]
				
				for j_seg, seg in enumerate(eqsec):
					
					if gbar_scaling == 'area':
						# conserves ratio in each segment but not total original conductance
						scale = cluster.or_area/cluster.eq_area
					
					elif gbar_scaling == 'gbar_integral':
						# does not conserve ratio but conserves gtot_or since: sum(g_i*area_i * or_area/eq_area) = or_area/eq_area * sum(gi*area_i) ~= or_area/eq_area * g_avg*eq_area = or_area*g_avg
						scale = or_gtot/eq_gtot
					
					else:
						raise Exception("Unknown gbar scaling method'{}'.".format(gbar_scaling))
					
					# Set gbar
					gval = getattr(seg, gname) * scale
					seg.__setattr__(gname, gval)
					
					cluster.eq_gbar[gname][j_seg] = gval # save for reconstruction

		# Check gbar calculation
		for gname in active_glist:
			gtot_or = cluster.or_gtot[gname]
			gtot_eq_scaled = sum(getattr(seg, gname)*seg.area() for seg in eqsec)
			# logger.debug("conductance %s : gtot_or = %.3f ; gtot_eq = %.3f",
			# 				gname, gtot_or, gtot_eq_scaled)

		# Debugging info:
		logger.debug("Created equivalent Section '%s' with \n\tL\tdiam\tcm\tRa\tnseg"
						"\n\t%.3f\t%.3f\t%.3f\t%.3f\t%d\n", 
						clusters[i_clu].label, eqsec.L, eqsec.diam, eqsec.cm, eqsec.Ra, eqsec.nseg)
	
	# Substitute equivalent section into tree
	logger.info("\n###############################################################"
				"\nSubstituting equivalent sections...\n")
	
	for i_clu, cluster in enumerate(clusters):
		eqsec = eq_secs[i_clu]
		# Disconnect substituted segments and attach segment after Y boundary
		# Can only do this now since all paths need to be walkable before this
		logger.debug("Substituting zipped section {0}".format(eqsec))
		redtools.sub_equivalent_Y_sec(eqsec, cluster.parent_seg, cluster.bound_segs, 
					orsecrefs, mechs_chans, delete_substituted=True)
		logger.debug("Substitution complete.\n")

	# build new list of valid SectionRef
	newsecrefs = [ref for ref in orsecrefs if not (ref.is_substituted or ref.is_deleted)]
	newsecrefs.extend(eq_refs)
	return eq_refs, newsecrefs


################################################################################
# Utility functions
################################################################################


def assign_attributes(noderef, allsecrefs, attr_dict):
	"""
	Assign attributes to all SectionRefs in subtree of given section

	@param attr_dict	dictionary of key-value pairs (attribute_name, attribute_value)
	"""
	# Assign current node
	for aname, aval in attr_dict.items():
		setattr(noderef, aname, aval)

	childsecs = noderef.sec.children()
	childrefs = [getsecref(sec, allsecrefs) for sec in childsecs]
	for childref in childrefs:
		assign_attributes(childref, allsecrefs, attr_dict)


def assign_sth_indices(noderef, allsecrefs, parref=None):
	"""
	Re-assign tree index and table index as labelled in Gillies & Willshaw code
	(see sth-data folder).
	"""
	if not hasattr(noderef, 'tree_index'):
		noderef.tree_index = parref.tree_index
	
	if not hasattr(noderef, 'table_index'):
		noderef.table_index = -1 # unassigned

	# assign a unique cell GID
	if not hasattr(noderef, 'gid'):
		if noderef.table_index >= 0:
			noderef.gid = min(0,noderef.tree_index)*100 + noderef.table_index
		else:
			noderef.gid = noderef.zip_id

	childsecs = noderef.sec.children()
	childrefs = [getsecref(sec, allsecrefs) for sec in childsecs]
	for childref in childrefs:
		assign_sth_indices(childref, allsecrefs, parref=noderef)


################################################################################
# Interface implementations (reduce_cell.py)
################################################################################


def preprocess_impl(reduction):
	"""
	Preprocess cell for Marasco reduction.

	(Implementation of interface declared in reduce_cell.CollapseReduction)

	@param	reduction		reduce_cell.CollapseReduction object
	"""

	allsecrefs = reduction.all_sec_refs

	dendL_secs = list(h.SThcell[0].dend0)
	dendR_secs = list(h.SThcell[0].dend1)
	dend_lists = [dendL_secs, dendR_secs]

	# Assign indices used in Gillies code (sth-data folder)
	for somaref in reduction._soma_refs:
		somaref.tree_index = -1
		somaref.table_index = 0

	for secref in reduction._dend_refs:
		for i_dend, dendlist in enumerate(dend_lists):
			if any([sec.same(secref.sec) for sec in dendlist]):
				secref.tree_index = i_dend
				secref.table_index= next((i+1 for i,sec in enumerate(dendlist) if sec.same(secref.sec)))
		
	# Assign unique GID to each Section
	assign_sth_indices(reduction._root_ref, allsecrefs)

	# Calculate section path properties for entire tree
	for secref in allsecrefs:
		# Assign path length, path resistance, electrotonic path length to each segment
		redtools.sec_path_props(secref, f_lambda, gleak_name)


	# Choose stereotypical path for interpolation
	interp_tree_id = 1
	interp_table_ids = (1,3,8)
	path_secs = [secref for secref in allsecrefs if (secref.tree_index == interp_tree_id and 
													secref.table_index in interp_table_ids)]

	# Compute properties along this path
	sec_props = ['pathL0', 'pathL1', 'pathri0', 'pathri1', 'pathLelec0', 'pathLelec1']
	seg_props = ['pathL_seg', 'pathri_seg', 'pathL_elec']
	mechs_chans = reduction.mechs_gbars_dict
	reduction.path_props = [redtools.get_sec_props_obj(ref, mechs_chans, seg_props, sec_props) for ref in path_secs]


def prepare_folds_impl(reduction):
	"""
	Prepare next collapse operation: assign topology information
	to each Section.

	(Implementation of interface declared in reduce_cell.CollapseReduction)
	"""	

	root_ref = reduction._root_ref
	soma_refs = reduction._soma_refs
	allsecrefs = reduction.all_sec_refs

	# Calculate/assign properties used in calculation
	assign_sth_indices(root_ref, allsecrefs)
	assign_attributes(root_ref, allsecrefs, {'max_passes': 100})

	logger.info("\n###############################################################"
				"\nAssigning topology & path properties ...\n")

	# Assign topology info (order, level, strahler number)
	for fold_root in reduction._fold_root_refs:
		clutools.assign_topology_attrs(fold_root, allsecrefs)

	dendL_juction = getsecref(h.SThcell[0].dend0[0], allsecrefs)
	dendL_upper_root = getsecref(h.SThcell[0].dend0[1], allsecrefs)

	dendL_juction.order = 1
	dendL_juction.level = 0
	dendL_juction.strahlernumber = dendL_upper_root.strahlernumber+1
	
	for somaref in soma_refs:
		somaref.order = 0
		somaref.level = 0
		somaref.strahlernumber = dendL_juction.strahlernumber

	# Assign path properties
	for secref in allsecrefs:
		# Calculate path length, path resistance, electrotonic path length to each segment
		redtools.sec_path_props(secref, f_lambda, gleak_name)


def calc_folds_impl(reduction, i_pass, Y_criterion='highest_level'):
	"""
	Collapse branches at branch points identified by given criterion.
	"""
	allsecrefs = reduction.all_sec_refs

	# Find collapsable branch points
	target_Y_secs = folding.find_collapsable(allsecrefs, i_pass, Y_criterion)

	# Do collapse operation at each branch points
	clusters = calc_collapses(target_Y_secs, i_pass, allsecrefs)

	# Save results
	reduction.clusters = clusters


def make_folds_impl(reduction):
	"""
	Make equivalent Sections for branches that have been folded.
	"""

	# Mark Sections
	for secref in reduction.all_sec_refs:
		secref.is_substituted = False
		secref.is_deleted = False

	# Make new Sections
	eq_refs, newsecrefs = sub_equivalent_Y_sections(reduction.clusters, 
							reduction.all_sec_refs, reduction.path_props,
							interp_prop='path_L', interp_method='linear_neighbors', 
							gbar_scaling='area')

	# Set ion styles
	for sec in h.allsec(): # makes each section the CAS
		h.ion_style("na_ion",1,2,1,0,1)
		h.ion_style("k_ion",1,2,1,0,1)
		h.ion_style("ca_ion",3,2,1,1,1)

	reduction.update_refs(dend_refs=eq_refs) # prepare for next iteration


def postprocess_impl(reduction):
	"""
	Post-process cell for Marasco reduction.

	(interface declared in reduce_cell.CollapseReduction)

	@param	reduction		reduce_cell.CollapseReduction object
	"""
	# Tweaking
	tweak_funcs = reduction.get_reduction_param(reduction.active_method, 'post_tweak_funcs')
	for func in tweak_funcs:
		func(reduction)

	# Assign identifiers (for synapse placement etc.)
	assign_sth_indices(reduction._root_ref, reduction.all_sec_refs)

	# Assign topology info (order, level, strahler number)
	for fold_root in reduction._fold_root_refs:
		clutools.assign_topology_attrs(fold_root, reduction.all_sec_refs)


################################################################################
# Reduction Experiments
################################################################################


def reduce_gillies_incremental(n_passes, zips_per_pass):
	"""
	Reduce Gillies & Willshaw STN neuron model by incrementally
	collapsing (zipping) Y-sections in each of 3 identical trees.

	ALGORITHM
	- cluster each of 3 identical trees
	- merge/collapse segments in each cluster
	- create equivalent sections
		- conductances are interpolated according to distance L/lambda from soma
	"""
	############################################################################
	# 0. Load full model to be reduced (Gillies & Willshaw STN)
	# for sec in h.allsec():
	# 	if not sec.name().startswith('SThcell') and delete_old_cells:
	# 		h.delete_section() # delete existing cells
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
	dendL_upper_root = getsecref(h.SThcell[0].dend0[1], dendLrefs) # root section of upper left dendrite
	dendL_lower_root = getsecref(h.SThcell[0].dend0[2], dendLrefs) # root section of lower left dendrite

	############################################################################
	#  Pre-clustering: calculate properties

	# Assign indices used in Gillies code (sth-data folder)
	somaref.tree_index = -1
	somaref.table_index = 0
	for j, dendlist in enumerate((dendLrefs, dendRrefs)):
		for i, secref in enumerate(dendlist):
			secref.tree_index = j # left tree is 0, right is 1
			secref.table_index = i+1 # same as in /sth-data/treeX-nom.dat
	
	assign_sth_indices(somaref, allsecrefs) # assign cell GIDs

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

	############################################################################
	# Iterative/recursive collapsing

	eq_secrefs = []
	for i_pass in range(n_passes):
		# Check that we haven't collapsed too many levels
		if not (dendL_upper_root.exists() and dendL_lower_root.exists()):
			logger.warning("Maximum number of collapses reached: Cannot collapse past trunk sections.")
			break

		logger.info("\n###############################################################"
					"\nStarting reduction pass {0}/{1}...\n".format(i_pass+1, n_passes))

		############################################################################
		# 0. Pre-clustering: calculate properties

		# Assign attributes used in calculation
		assign_sth_indices(somaref, allsecrefs)
		assign_attributes(somaref, allsecrefs, {'max_passes': 100})

		# Assign Strahler numbers
		logger.info("\n###############################################################"
					"\nAssigning Strahler's numbers...\n")
		
		clutools.assign_topology_attrs(dendR_root, allsecrefs)
		clutools.assign_topology_attrs(dendL_upper_root, allsecrefs)
		clutools.assign_topology_attrs(dendL_lower_root, allsecrefs)
		
		dendL_juction.order = 1
		dendL_juction.level = 0
		dendL_juction.strahlernumber = dendL_upper_root.strahlernumber+1
		
		somaref.order = 0
		somaref.level = 0
		somaref.strahlernumber = dendL_juction.strahlernumber

		# Assign path properties
		for secref in allsecrefs:
			# Calculate path length, path resistance, electrotonic path length to each segment
			redtools.sec_path_props(secref, f_lambda, gleak_name)

		############################################################################
		# 1. Merge a number of Y-sections
		logger.info("\n###############################################################"
					"\nFinding Y-sections to merge...\n")

		# Find collapsable branch points
		target_Y_secs = find_collapsable(allsecrefs, i_pass, 'highest_level', zips_per_pass)

		# Do collapse operation at each branch points
		clusters = calc_collapses(target_Y_secs, i_pass, allsecrefs)

		############################################################################
		# 2. Create equivalent sections

		logger.info("\n###############################################################"
					"\nCreating equivalent sections...\n")
		# Mark Sections
		for secref in allsecrefs:
			secref.is_substituted = False
			secref.is_deleted = False

		eq_refs, newsecrefs = sub_equivalent_Y_sections(clusters, allsecrefs, path_props,
							interp_prop='path_L', interp_method='linear_neighbors', 
							gbar_scaling='area')

		allsecrefs = newsecrefs # prepare for next iteration
		eq_secrefs.extend(eq_refs)

		############################################################################
		# 3. Finalize
		
		# Set ion styles
		for sec in h.allsec(): # makes each section the CAS
			h.ion_style("na_ion",1,2,1,0,1)
			h.ion_style("k_ion",1,2,1,0,1)
			h.ion_style("ca_ion",3,2,1,1,1)

		# Print tree structure
		logger.info("Equivalent tree topology:")
		if logger.getEffectiveLevel() <= logging.DEBUG:
			h.topology()

	# Make sure all sections have identifiers
	assign_sth_indices(somaref, allsecrefs)

	eq_refs = [ref for ref in eq_secrefs if ref.exists()]
	return eq_refs, newsecrefs



if __name__ == '__main__':
	# clusters, eq_secs = reduce_gillies_partial(delete_old_cells=True)
	eq_refs, newsecrefs = reduce_gillies_incremental(n_passes=7, zips_per_pass=100)
	from neuron import gui # check in ModelView: conductance distribution, structure