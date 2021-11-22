"""
Reduce Gillies & Willshaw (2006) STN neuron model using the method
described in Marasco & Migliore (2012)


@author Lucas Koelman
@date	5-12-2016
"""

# Python modules
import sys
import re
import math
PI = math.pi

# logging of DEBUG/INFO/WARNING messages
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
logname = "reduction" # __name__
logger = logging.getLogger(logname) # create logger for this module
fmtr = logging.Formatter('%(levelname)s:%(message)s')
# Log to file
# fh = logging.FileHandler('reduce_marasco.log')
# fh.setFormatter(fmtr)
# logger.addHandler(fh)
# Log to stream
# ch = logging.StreamHandler(sys.stdout)
# ch.setFormatter(fmtr)
# logger.addHandler(ch)

import os.path
scriptdir, scriptfile = os.path.split(__file__)
modulesbase = os.path.normpath(os.path.join(scriptdir, '..'))
sys.path.append(modulesbase)

# NEURON modules
import neuron
h = neuron.h

# Load NEURON function libraries
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# Load NEURON mechanisms
# add this line to nrn/lib/python/neuron/__init__.py/load_mechanisms()
# from sys import platform as osplatform
# if osplatform == 'win32':
# 	lib_path = os.path.join(path, 'nrnmech.dll')
NRN_MECH_PATH = os.path.normpath(os.path.join(scriptdir, 'nrn_mechs'))
neuron.load_mechanisms(NRN_MECH_PATH)

# Own modules
import redutils as redtools
from redutils import ExtSecRef, EqProps, getsecref, lambda_AC, prev_seg, seg_index, same_seg # for convenience
import cluster as clutools
from cluster import Cluster
import interpolation as interp
import reduce_bush_sejnowski as redbush
import analyze_reduction as analysis

# Global variables (convert to class members in future)
gillies_mechs_chans = {'STh': ['gpas'], # passive/leak channel
				'Na': ['gna'], 'NaL': ['gna'], # Na channels
				'KDR': ['gk'], 'Kv31': ['gk'], 'sKCa':['gk'], # K channels
				'Ih': ['gk'], # nonspecific channels
				'CaT': ['gcaT'], 'HVA': ['gcaL', 'gcaN'], # Ca channels
				'Cacum': []} # No channels

mechs_chans = gillies_mechs_chans
gleak_name = 'gpas_STh'
glist = [gname+'_'+mech for mech,chans in mechs_chans.items() for gname in chans]
f_lambda = 100.0

def merge_parallel(childrefs, allsecrefs):
	"""
	Merge parallel sections at branch point using <br> equations in Marasco (2012)

	ALGORITHM
	- `L_br = sum(S_i*L_i)/sum(S_i)` where S_I is the area of branch i
	- `diam_br = sqrt(sum(diam_i^2))`
	- `r_a,br = prod(r_a,i)/sum(r_a,i)` where r_a is the axial resistance Ri
	"""
	# Initialize combined properties of branched sections (children)
	L_br = 0.
	diam_br = 0.
	Ra_br = 0.
	rin_br = 1.
	eqsurf_sum = 0.
	ri_sum = 0.
	gtot_br = dict((gname, 0.0) for gname in glist) # sum of gbar multiplied by area
	cmtot_br = 0. # sum of cm multiplied by area

	# Update combined properties using the dimensions/electrical properties
	# of each child according to <br> expressions in Marasco (2012) 
	for childref in childrefs:
		# get equivalent child properties
		L_child, diam_child, Ra_child, ri_child, cmtot_child, gtot_child = merge_sequential(childref, allsecrefs)
		eqsurf_child = PI*diam_child*L_child

		# combine according to <br> (parallel) expressions
		eqsurf_sum += eqsurf_child
		L_br += eqsurf_child*L_child # LENGTH: eq (1) - weight by area
		diam_br += diam_child**2 # RADIUS: eq (2) - 2-norm of radii
		Ra_br += Ra_child # SPECIFIC AXIAL RESISTANCE - average Ra
		rin_br *= ri_child # ABSOLUTE AXIAL RESISTANCE - parallel conductances
		ri_sum += ri_child # need sum in final calculation

		# Distributed properties
		cmtot_br += cmtot_child
		for gname in glist:
			gtot_br[gname] += gtot_child[gname]

		# mark child as absorbed
		childref.absorbed = True
		childref.visited = True

	# Finalize <br> calculation (MUST BE VALID IF ONLY ONE CHILDREF)
	L_br /= eqsurf_sum # eq. L_br
	diam_br = math.sqrt(diam_br) # eq. rho_br
	Ra_br = Ra_br/len(childrefs) # average Ra (NOTE: unused, cluster Ra calculated from equation Ra^{eq})
	cross_br = PI*diam_br**2/4. # cross-section area
	rax_br = Ra_br*L_br/cross_br/100. # absolute axial resistance of section with equivalent dimensions and Ra
	if len(childrefs) > 1:
		rin_br /= ri_sum # eq. r_a,br: product/sum

	return L_br, diam_br, Ra_br, rin_br, rax_br, cmtot_br, gtot_br

def merge_sequential(rootref, allsecrefs):
	""" 
	Merge sequential sections into one equivalent sections
	"""
	# Get references
	rootsec = rootref.sec
	allchildrefs = [getsecref(sec, allsecrefs) for sec in rootsec.children()]
	childrefs = [ref for ref in allchildrefs if ref.cluster_label==rootref.cluster_label]

	# Gather properties for root sec
	L_root = rootsec.L
	diam_root = rootsec.diam
	Ra_root = rootsec.Ra
	ri_root = sum(seg.ri() for seg in rootsec) # absolute axial resistance between 0-1 ends

	# Distributed properties
	cmtot_seq = sum(seg.cm*seg.area() for seg in rootsec) # sum of cm multiplied by area
	gtot_seq = dict((gname, 0.0) for gname in glist) # sum of gbar multiplied by area
	for gname in glist:
		gtot_seq[gname] += sum(getattr(seg, gname)*seg.area() for seg in rootsec)

	# Handle leaf sections
	if not any(childrefs):
		return L_root, diam_root, Ra_root, ri_root, cmtot_seq, gtot_seq

	# Combine properties of parallel branched sections (children)
	L_br, diam_br, Ra_br, rin_br, rax_br, cmtot_br, gtot_br = merge_parallel(childrefs, allsecrefs)

	# use <seq> expressions in Marasco (2012) to merge equivalent child into parent
	L_seq = L_root + L_br # L_seq equation
	Ra_seq = (Ra_root + Ra_br)/2.
	ri_seq = ri_root + rin_br # var 'newri2' in Marasco code used for ri calculation
	rax_seq = ri_root + rax_br # var 'newri' in Marasco code used for diam calculation
	diam_seq = math.sqrt(Ra_seq*L_seq*4./PI/rax_seq/100.) # rho_seq equation (conserves ri_seq)
	
	# Keep track of total conductance/capacitance
	cmtot_seq += cmtot_br
	for gname in glist:
		gtot_seq[gname] += gtot_br[gname]

	return L_seq, diam_seq, Ra_seq, ri_seq, cmtot_seq, gtot_seq

def collapse_subtree(rootref, allsecrefs):
	""" 
	Recursively merge within-cluster connected sections in subtree
	of the given node using <br> and <seq> expressions in Marasco (2012).
	"""
	# Collapse is equal to sequential merge of the root and equivalent parallel circuit of its children
	return merge_sequential(rootref, allsecrefs)

def collapse_seg_subtree(rootseg, allsecrefs, eligfunc=None, modfunc=None, bound_segs=None):
	""" 
	Recursively merge within-cluster connected sections in subtree
	of the given segment using equations <br> and <seq> in Marasco (2012).

	@pre		in the tree, all connections between sections must be made
				according to the convention that the 0-end of the child
				connects to the 1-end of the parent

	@param eligfunc		function(seg,jseg,secref) -> bool indicating whether
						the segment is eligible to be absorbed (seg is the segment,
						jseg is its index, secref is a SectionRef to its section)

	ALGORITHM
	- for each child segment: recursively call `equivalent_properties = collapse_subtree(child)`
	- then combine the equivalent properties: absorb into current rootseg and return
	"""
	# Get segment info
	rootsec = rootseg.sec
	rootref = getsecref(rootsec, allsecrefs)
	co_segs = [seg for seg in rootsec]
	i_rootseg = seg_index(rootseg)
	# logger.debug("Collapsing children of segment {0} (segment {1}/{2})".format(rootseg, i_rootseg+1, rootsec.nseg))

	# Calculate properties of root segment
	L_root = rootsec.L/rootsec.nseg
	diam_root = rootseg.diam
	Ra_root = rootsec.Ra
	Ri_root = rootseg.ri() # axial resistance between start-end (0-1)

	# 1. Gather children to absorb #############################################
	child_refs = [getsecref(sec, allsecrefs) for sec in rootsec.children()]

	# Function for checking if segments are absorbable
	if eligfunc is None:
		rootlabel = rootref.cluster_labels[i_rootseg]
		absorbable = lambda seg, jseg, ref: ref.cluster_labels[jseg]==rootlabel
	else:
		absorbable = eligfunc
	if modfunc is None:
		modfunc = lambda seg, jsef, ref: None

	# Get aborbable child segments
	child_segs = []
	root_is_bound = False # whether current root is a boundary segment
	if i_rootseg == rootsec.nseg-1:
		# IF END SEGMENT: get child segments that are in same cluster
		for secref in child_refs:
			childseg = next(seg for seg in secref.sec) # get the first segment
			if absorbable(childseg, 0, secref):
				child_segs.append(childseg)
				# mark segment
				secref.visited[0] = True
				secref.absorbed[0] = True
				modfunc(childseg, 0, secref)
			else:
				root_is_bound = True
		if not any(child_segs):
			root_is_bound = True

	else:
		# IF NOT END SEGMENT: get adjacent segment if in same cluster
		j_seg = i_rootseg+1
		childseg = co_segs[j_seg]
		if absorbable(childseg, j_seg, rootref):
			child_segs.append(childseg)
			# mark segment
			rootref.visited[j_seg] = True
			rootref.absorbed[j_seg] = True
			modfunc(childseg, j_seg, rootref)
		else:
			root_is_bound = True

	# Update boundary segments (furthest collapsed segment)
	if root_is_bound and (bound_segs is not None):
		bound_segs.append(rootseg)

	# 2. Recursively call collapse on children #################################

	# Base Case: leaf segments (i.e. no child segments in same cluster)
	if not any(child_segs):
		# logger.debug("No child segments to absorb: ending recursive descent.")
		eq_props_seq = EqProps(L_eq=L_root, diam_eq=diam_root, Ra_eq=Ra_root, Ri_eq=Ri_root)
		eq_props_br = EqProps(L_eq=0., diam_eq=0., Ra_eq=0., Ri_eq=0.)
		return eq_props_seq, eq_props_br

	# Parallel merge of child properties (use <br> equations Marasco (2012))
	# NOTE: if only one child, <br> equations must yield properties of that child
	L_br = 0.
	diam_br = 0.
	Ra_br = 0.
	Ri_br = 1.
	area_br = 0.
	cross_area_br = 0.
	Ri_br_sum = 0. # sum of axial resistances of parallel child branches
	for child in child_segs:
		# Recursive call (if not absorbable, child segments were not added)
		# logger.debug("Absorbing child segment {0}...".format(child))
		cp, _ = collapse_seg_subtree(child, allsecrefs, eligfunc, modfunc, bound_segs)
		ch_area = PI*cp.diam_eq*cp.L_eq

		# Update equivalent properties of all child branches in parallel
		area_br += ch_area
		L_br += ch_area*cp.L_eq		# Length: Eq. L_br
		diam_br += cp.diam_eq**2	# Diameter: Eq. rho_br
		Ra_br += cp.Ra_eq			# Axial resistivity: average Ra
		Ri_br *= cp.Ri_eq			# Axial resistance: Eq. r_a,br
		Ri_br_sum += cp.Ri_eq		# Sum of branch axial resistances

	# Finalize <br> calculation
	L_br /= area_br
	diam_br = math.sqrt(diam_br)
	Ra_br = Ra_br/len(child_segs) # NOTE: unused, cluster Ra calculated from equation Ra^{eq}
	cross_area_br = PI*diam_br**2/4.
	if len(child_segs) > 1:
		Ri_br /= Ri_br_sum # Must be valid parallel circuit of Ri if only one child
	Ri_br_eqgeom = Ra_br*L_br/cross_area_br/100. # Axial resistance of section with geometry equal to merged properties
	eq_props_br = EqProps(L_eq=L_br, diam_eq=diam_br, Ra_eq=Ra_br, Ri_eq=Ri_br_eqgeom)

	# 3. Use solution of recursive call to solve ##############################

	# Sequential merge of root & merged child properties (use <seq> equations Marasco (2012))
	L_seq = L_root + L_br			# Eq. L_seq
	Ra_seq = (Ra_root + Ra_br)/2.	# for rho_seq calculation
	Ri_seq = Ri_root + Ri_br		# Eq. r_a,seq
	Ri_seq_eqgeom = Ri_root + Ri_br_eqgeom						# for diam_seq calculation
	diam_seq = math.sqrt(Ra_seq*L_seq*4./PI/Ri_seq_eqgeom/100.)	# Eq. rho_seq
	eq_props_seq = EqProps(L_eq=L_seq, diam_eq=diam_seq, Ra_eq=Ra_seq, Ri_eq=Ri_seq)
	return eq_props_seq, eq_props_br


def map_pathri_gbar(cluster, allsecrefs):
	"""
	Map axial path resistance values to gbar values in the cluster

	ALGORITHM:
	- for each section:
		- for each segment in section, save gbar and axial path resistance
			- axial path resistance obtained by interpolating pathri0 & pathri1
	"""
	clu_secs = [secref for secref in allsecrefs if secref.cluster_label==cluster.label]

	# Keep a dict that maps gname to a collection of data points (pathri, gbar)
	cluster.pathri_gbar = dict((gname, []) for gname in glist)
	for secref in clu_secs:
		for seg in secref.sec:
			for gname in glist:
				if not hasattr(seg, gname):
					continue # section doesn't have ion channel: skip
				# Add data point
				seg_pathri = secref.pathri0 + seg.x*(secref.pathri1-secref.pathri0)
				seg_gbar = getattr(seg, gname)
				cluster.pathri_gbar[gname].append((seg_pathri, seg_gbar))

def calc_gbar(cluster, gname, pathri):
	"""
	Calculate gbar for a point on the equivalent section of given cluster
	given the axial path resistance to that point.
	"""
	gbar_pts = cluster.pathri_gbar[gname] # list of (pathri, gbar) data points
	gbar_pts.sort(key=lambda pt: pt[0]) # sort by pahtri ascending
	eq_path_x = (pathri - cluster.pathri0) / (cluster.pathri1 - cluster.pathri0)
	if eq_path_x <= 0.:
		return gbar_pts[0]
	elif eq_path_x >= 1.:
		return gbar_pts[-1]
	else:
		# average of points that lie within pathri +/- 11% of segment axial resistance
		deltari = 0.11*(cluster.pathri1-cluster.pathri0)
		surr_pts = [pt for pt in gbar_pts if (pt[0] >= pathri-deltari) and (pt[0] <= pathri+deltari)]
		if not any(surr_pts):
			# take average of two closest points
			gbar_pts.sort(key=lambda pt: abs(pt[0]-pathri))
			surr_pts = [pt for i, pt in enumerate(gbar_pts) if i < 2]
		return sum(pt[1] for pt in surr_pts)/len(surr_pts) # take average


def has_cluster_parentseg(aseg, cluster, clu_secs):
	"""
	Utility function to check if segment has parent segment in same cluster
	"""
	parseg = prev_seg(aseg)
	if parseg is None:
		return False # no parent segment
	parref = getsecref(parseg.sec, clu_secs)
	if parref is None:
		return False # parent section has no segments in same cluster
	j_parseg = seg_index(parseg)
	return parref.cluster_labels[j_parseg]==cluster.label


def zip_fork_branches(allsecrefs, i_pass, zips_per_pass, Y_criterion):
	"""
	Find forks (Y-sections) that can be 'zipped' and calculate properties
	of equivalent zipped section.

	@return			list of Cluster objects with properties of equivalent
					Section for each 'zipped' fork.
	"""
	# Calculate section statistics
	for secref in allsecrefs:
		# Calculate path length, path resistance, electrotonic path length to each segment
		redtools.sec_path_props(secref, f_lambda, gleak_name)
		# Flag all sections as unvisited
		if not hasattr(secref, 'absorbed'):
			secref.absorbed = [False] * secref.sec.nseg
			secref.visited = [False] * secref.sec.nseg
		secref.zip_labels = [None] * secref.sec.nseg

	# Find Y-sections
	if Y_criterion=='max_electrotonic':
		# Find Y with longest electrotonic length to first Y among its children.
		# This corresponds to  Y with child section that has longest L/lambda.
		# Collapsing this Y will eliminate most compartments.
		for secref in allsecrefs:
			child_secs = secref.sec.children()
			secref.collapsable_L_elec = 0.0
			if any(child_secs):
				min_child_L = min(sec.L for sec in child_secs)
			for chi_sec in child_secs:
				# Get segment that is that distance away from child
				furthest_seg = chi_sec(min_child_L/chi_sec.L)
				furthest_L_elec = redtools.seg_path_L_elec(furthest_seg, f_lambda, gleak_name)
				secref.collapsable_L_elec += (secref.pathLelec1 - furthest_L_elec)

		# Find other Y sections that meet same criteria
		max_L_collapsable = max(ref.collapsable_L_elec for ref in allsecrefs)
		low, high = 0.95*max_L_collapsable, 1.05*max_L_collapsable
		candidate_Y_secs = [ref for ref in allsecrefs if (any(ref.sec.children())) and (
															low<ref.collapsable_L_elec<high)]
		target_level = max(ref.level for ref in candidate_Y_secs)
		target_Y_secs = [ref for ref in candidate_Y_secs if ref.level==target_level]

		# Report findings
		logger.debug("The maximal collapsable electrotonic length is %f", max_L_collapsable)
		logger.debug("Found %i sections with collapsable length within 5%% of this value", 
						len(candidate_Y_secs))
		logger.debug(("The highest level of a parent node with this value is %i, "
						"and the number of nodes at this level with the same value is %i"), 
						target_level, len(target_Y_secs))

	elif Y_criterion=='highest_level':
		max_level = max(ref.level for ref in allsecrefs)
		target_Y_secs = [ref for ref in allsecrefs if (ref.level==max_level-1) and (
						 i_pass+1 <= ref.max_passes) and (len(ref.sec.children()) >= 2)]

	else:
		raise Exception("Unknow Y-section selection criterion '{}'".format(Y_criterion))

	# Don't collapse more Y-sections than given maximum
	n_Y = len(target_Y_secs)
	logger.debug("Found {0} Y-sections that meet selection criterion. Keeping {1}/{2}\n\n".format(
					n_Y, min(n_Y, zips_per_pass), n_Y))
	target_Y_secs = [ref for i,ref in enumerate(target_Y_secs) if i<zips_per_pass]

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
		eq_seq, eq_br = collapse_seg_subtree(par_sec(1.0), allsecrefs, eligfunc, 
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


def merge_seg_cluster(cluster, allsecrefs, average_Ri):
	"""
	Merge cluster of segments (use for segment-based clustering)

	@param average_Ri	see param 'average_Ri' in merge_sec_cluster()

	@post	cluster will have following attributes, calculated from member segments:

			or_area		original surface area of all member segments
			or_cmtot	original summed capacitance (nonspecific) of all member segments
			or_gtot		original maximum conductance (nonspecific) of all member segments,
						for all inserted conductance mechanisms
			
			eqL			equivalent length
			eqdiam		equivalent diameter
			eqri		equivalent axial resistance Ri (absolute, in Mega Ohms)
			eqRa		equivalent axial cytoplasmic resistivity
			eq_area		equivalent area based on equivalent dimensions and passive properties
			eq_area_sum	sum of equilvalent surfaces of 'islands' of segments

	ALGORITHM
	
	- find root segment of next subtree in cluster 
		- subtree = set of connected segments with same cluster label
		- root = segment without a parent in the same cluster
	
	- collapse subtree using <br> and <seq> expressions Marasco (2012)
		- collapsing is a sequence of parallel-> sequential merge operations
			- starting at leaf nodes and progressing to root
		- first merge parallel branches at branch point (<br> expressions)
		- then do sequential merge of equivalent section and parent section (<seq> expressions)
	
	- then combine equivalent sections for merge subtrees using <eq> expressions
		- this yields final equivalent section properties for cluster

	"""
	# Gather all cluster sections & segments
	clu_secs = [secref for secref in allsecrefs if (cluster.label in secref.cluster_labels)]
	clu_segs = []
	for ref in allsecrefs:
		# Flag all sections as unvisited
		if not hasattr(ref, 'absorbed'):
			ref.absorbed = [False] * ref.sec.nseg
		if not hasattr(ref, 'visited'):
			ref.visited = [False] * ref.sec.nseg

		# Gather segments that are member of current cluster
		for i, seg in enumerate(ref.sec):
			if ref.cluster_labels[i] == cluster.label:
				clu_segs.append(seg)

	# Calculate axial path resistance to all segments in cluster
	clu_seg_pathri = []
	for secref in clu_secs:
		# Assign axial path resistances
		redtools.sec_path_ri(secref, store_seg_ri=True) # stores pathri on SectionRefs
		# Store pathri to start of segments in cluster
		segs_pathri = [pathri for i, pathri in enumerate(secref.pathri_seg) if (
						secref.cluster_labels[i]==cluster.label)]
		# Also store pathri to end of most distal segment
		if secref.cluster_labels[-1]==cluster.label:
			segs_pathri.append(secref.pathri1)
		clu_seg_pathri.extend(segs_pathri)

	# Get min & max path resistance in cluster (full model)
	cluster.orMaxpathri = max(clu_seg_pathri)
	cluster.orMinpathri = min(clu_seg_pathri)
	
	# Initialize equivalent properties
	cluster.eqL = 0.
	cluster.eqdiam = 0.
	cluster.eqri = 0.
	cluster.eq_area_sum = 0. # sum of surface calculated from equivalent dimensions

	# Initialize original properties
	cluster.or_area = sum(seg.area() for seg in clu_segs)
	cluster.or_cmtot = sum(seg.cm*seg.area() for seg in clu_segs)
	cluster.or_gtot = dict((gname, 0.0) for gname in glist)
	for gname in glist:
		cluster.or_gtot[gname] += sum(getattr(seg, gname)*seg.area() for seg in clu_segs)

	# Make generator that finds root segments
	root_gen = (seg for ref in clu_secs for i, seg in enumerate(ref.sec) if (
					not ref.visited[i] and not has_cluster_parentseg(seg, cluster, clu_secs)))

	# Start merging algorithm
	num_roots = 0
	while True:
		# Find next root segment in cluster (i.e. with no parent segment in same cluster)
		rootseg = next(root_gen, None)
		if rootseg is None:
			break # No more root segment left
		rootref = getsecref(rootseg.sec, allsecrefs)

		# Collapse subtree
		logger.debug("Collapsing subtree of cluster root segment %s", repr(rootseg))
		eq_props, eq_br = collapse_seg_subtree(rootseg, allsecrefs)

		# Combine properties of collapse sections using <eq> expressions
		eq_area = eq_props.L_eq * PI * eq_props.diam_eq
		cluster.eq_area_sum += eq_area
		cluster.eqL += eq_props.L_eq * eq_area	# eq. L^eq
		cluster.eqdiam += eq_props.diam_eq**2	# eq. rho^eq
		cluster.eqri += eq_props.Ri_eq			# eq. ra^eq

		# Mark as visited
		i_seg = seg_index(rootseg)
		rootref.visited[i_seg] = True
		num_roots += 1

	# Finalize <or> calculation
	cluster.or_cm = cluster.or_cmtot / cluster.or_area

	# Finalize <eq> calculation
	cluster.eqL /= cluster.eq_area_sum			# eq. L^eq
	cluster.eqdiam = math.sqrt(cluster.eqdiam)	# eq. rho^eq
	cluster.eqri /= num_roots					# eq. ra^eq
	if average_Ri:
		cluster.eqRa = PI*(cluster.eqdiam/2.)**2*cluster.eqri*100./cluster.eqL # eq. Ra^eq
	else:
		cluster.eqRa = sum(seg.sec.Ra for seg in clu_segs)/len(clu_segs) # alternative Ra^eq: average in cluster
	cluster.eq_area = cluster.eqL*PI*cluster.eqdiam # area_eq based on equivalent geometry

	# Debugging info
	logger.debug("Merged cluster '%s': equivalent properties are:\
		\n\teqL\teqdiam\teqRa\teqri\teq_area\
		\n\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		cluster.label, cluster.eqL, cluster.eqdiam, cluster.eqRa, cluster.eqri, cluster.eq_area)


def merge_sec_cluster(cluster, allsecrefs, average_Ri):
	"""
	Merge sections in cluster

	@param average_Ri	If True, Ra is calculated so that Ri (absolute axial
						resistance) of the equivalent section for each cluster is 
						the average Ri of all disconnected subtrees merged into that
						section. This does NOT conserve input resistance of the tree.

							If False, Ra is the average Ra in the cluster and the diameter 
						is calculated so that the absolute axial resistance Ri is 
						equivalent to the parallel circuit of all the unconnected 
						subtrees. This preserves input resistance

	ALGORITHM
	- find the next root of a within-cluster connected subtree
		- i.e. a section without a parent in the same cluster
		- this can be an isolated section (no mergeable children within cluster)
	- collapse subtree of that root section
	- update equivalent cluster properties using the <eq> expressions in Marasco (2012)

	"""
	# Gather sections in this cluster and mark/flag them
	clu_secs = [secref for secref in allsecrefs if secref.cluster_label==cluster.label]
	for sec in clu_secs:
		sec.absorbed = False
		sec.visited = False

	# Assign axial path resistance (sum of Ri (seg.ri()) along path)
	for sec in clu_secs:
		redtools.sec_path_ri(sec) # assigns pathri0/pathri1

	# Get min & max path resistance in cluster (full model)
	cluster.orMaxpathri = max(secref.pathri1 for secref in clu_secs)
	cluster.orMinpathri = min(secref.pathri0 for secref in clu_secs)
	
	# Initialize equivalent properties
	cluster.eqL = 0.
	cluster.eqdiam = 0.
	cluster.eqri = 0.
	cluster.eq_area_sum = 0. # sum of surface calculated from equivalent dimensions
	
	# Initialize original properties
	cluster.or_area = sum(sum(seg.area() for seg in secref.sec) for secref in clu_secs)
	cluster.or_cmtot = 0.
	cluster.or_gtot = dict((gname, 0.0) for gname in glist)

	# utility function to check if sec has parent within cluster
	def has_clusterparent(secref):
		return secref.has_parent() and (getsecref(secref.parent, clu_secs) is not None)

	# Find connected subtrees within cluster and merge/collapse them
	rootfinder = (sec for sec in clu_secs if (not sec.visited and not has_clusterparent(sec))) # compiles generator function
	for secref in rootfinder:
		rootref = clutools.clusterroot(secref, clu_secs) # make sure it is a cluster root

		# Collapse subtree
		logger.debug("Collapsing subtree of cluster root %s", repr(rootref))
		L_eq, diam_eq, Ra_eq, ri_eq, cmtot_eq, gtot_eq = collapse_subtree(rootref, allsecrefs)

		# Combine properties of collapse sections using <eq> expressions
		surf_eq = L_eq*PI*diam_eq
		cluster.eq_area_sum += surf_eq
		cluster.eqL += L_eq * surf_eq
		cluster.eqdiam += diam_eq**2
		cluster.eqri += ri_eq

		# Save distributed properties
		cluster.or_cmtot += cmtot_eq
		for gname in glist:
			cluster.or_gtot[gname] += gtot_eq[gname]

		# Mark as visited
		rootref.visited = True

	# Check each section either absorbed or rootsec
	# assert not any(not sec.absorbed and has_clusterparent(sec) for sec in clu_secs), (
	assert all(sec.absorbed or not has_clusterparent(sec) for sec in clu_secs), (
			'Each section should be either absorbed or be a root within the cluster')

	# Finalize <or> calculation
	cluster.or_cm = cluster.or_cmtot / cluster.or_area

	# Finalize <eq> calculation
	cluster.eqL /= cluster.eq_area_sum # LENGTH: equation L_eq
	cluster.eqdiam = math.sqrt(cluster.eqdiam) # RADIUS: equation rho_eq
	cluster.eqri /= sum(not sec.absorbed for sec in clu_secs) # ABSOLUTE AXIAL RESISTANCE: equation r_a,eq
	if average_Ri:
		cluster.eqRa = PI*(cluster.eqdiam/2.)**2*cluster.eqri*100./cluster.eqL # conserve eqri as absolute axial resistance
	else:
		cluster.eqRa = sum(secref.sec.Ra for secref in clu_secs)/len(clu_secs) # average Ra in cluster
	cluster.eq_area = cluster.eqL*PI*cluster.eqdiam # EQUIVALENT SURFACE

	# Debugging info
	logger.debug("Merged cluster '%s': equivalent properties are:\
		\n\teqL\teqdiam\teqRa\teqri\teq_area\
		\n\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		cluster.label, cluster.eqL, cluster.eqdiam, cluster.eqRa, cluster.eqri, cluster.eq_area)

def min_nseg_marasco(sec):
	""" Minimum number of segments based on electrotonic length """
	return int((sec.L/(0.1*lambda_AC(sec,100.))+0.9)/2)*2 + 1  

def equivalent_sections(clusters, allsecrefs, gradients):
	""" Create the reduced/equivalent cell by creating 
		a section for each cluster 

	@param clusters		list of Cluster objects containing data
						for each cluster

	@param allsecrefs	list of SectionRef (mutable) with first element
						a ref to root/soma section

	@param gradients	if True, gbar in each segment is determined from
						the average gbar at the same path resistance in the
						original model (i.e. nonuniform gbar). If False,
						a uniform gbar is used in each equivalent section.

	@return				list of SectionRef containing equivalent Section 
						for each cluster (in same order) as well as min
						and max path resistance for each cluster/section
						as properties pathri0/pathri1 on SectionRef objects
	"""
	# Create equivalent section for each clusters
	# eq_secs = [h.Section() for clu in clusters]
	for clu in clusters:
		h("create %s" % clu.label) # ensures references are not destroyed and names are clear
	eq_secs = [getattr(h, clu.label) for clu in clusters]
	eq_secrefs = [ExtSecRef(sec=sec) for sec in eq_secs]

	# Connect sections
	for i, clu_i in enumerate(clusters):
		for j, clu_j in enumerate(clusters):
			if clu_j is not clu_i and clu_j.parent_label == clu_i.label:
				eq_secs[j].connect(eq_secs[i], clu_j.parent_pos, 0)

	# Set dimensions, passive properties, active properties
	for i, secref in enumerate(eq_secrefs):
		sec = secref.sec
		sec.push() # Make section the CAS
		cluster = clusters[i]

		# Store some cluster properties on SectionRef
		secref.cluster_label = cluster.label
		secref.order = cluster.order

		# Set geometry 
		sec.L = cluster.eqL
		sec.diam = cluster.eqdiam
		sec_area = sum(seg.area() for seg in sec) # should be same as cluster eq_area
		surf_fact = cluster.or_area/cluster.eq_area # scale factor: ratio areas original/equivalent

		# Passive electrical properties (except Rm/gleak)
		sec.cm = cluster.or_cmtot / sec_area
		sec.Ra = cluster.eqRa

		# Set number of segments based on rule of thumb electrotonic length
		sec.nseg = redtools.min_nseg_hines(sec)

		# calculate min/max path resistance in equivalent section (cluster)
		pathri0, pathri1 = redtools.sec_path_ri(eq_secrefs[i])
		cluster.pathri0 = pathri0
		cluster.pathri1 = pathri1
		sec_ri = sum(seg.ri() for seg in sec)
		assert (pathri0 < pathri1), ('Axial path resistance at end of section '
									 'should be higher than at start of section')
		assert (-1e-6 < (pathri1 - pathri0) - sec_ri < 1e-6), ('absolute axial '
						'resistance not consistent with axial path resistance')

		# Insert all mechanisms and set conductances
		for mech in mechs_chans.keys():
			sec.insert(mech)
		for seg in sec:
			for gname in glist:
				if gradients:
					# Look for average gbar value at points with same path resistance
					seg_pathri = pathri0 + seg.x*(pathri1-pathri0)
					gval = calc_gbar(cluster, gname, seg_pathri)
				else:
					gval = cluster.or_gtot[gname] / sec_area # yields same sum(gbar*area) as in full model
				seg.__setattr__(gname, gval)
		
		# Re-scale gbar distribution to yield same total gbar (sum(gbar*area))
		if gradients:
			for gname in glist:
				gtot_eq = sum(getattr(seg, gname)*seg.area() for seg in sec)
				gtot_or = cluster.or_gtot[gname]
				if gtot_eq <= 0. : gtot_eq = 1.
				for seg in sec:
					# NOTE: old method, almost same but less accurate
					# seg.__setattr__(gname, getattr(seg, gname)*surf_fact)
					# NOTE: this conserves gtot since sum(g_i*area_i * gtot_or/gtot_eq) = gtot_or/gtot_eq*sum(gi*area_i) = gtot_or/gtot_eq*gtot_eq
					seg.__setattr__(gname, getattr(seg, gname)*gtot_or/gtot_eq)
				# Check calculation
				gtot_eq_scaled = sum(getattr(seg, gname)*seg.area() for seg in sec)
				logger.debug("== Conductance %s ===\nOriginal: gtot = %.3f\nEquivalent: gtot = %.3f",
								gname, gtot_or, gtot_eq_scaled)

		# Debugging info:
		logger.debug("Created equivalent section '%s' with \n\tL\tdiam\tcm\tRa\tri\tpathri0\tpathri1\tnseg\
		\n\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d", cluster.label, sec.L, sec.diam, sec.cm, sec.Ra, 
		sec_ri, pathri0, pathri1, sec.nseg)

		# Unset CAS
		h.pop_section()

	return eq_secs, eq_secrefs # return both or secs will be deleted

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

			# Set conductances by interpolating neighbors
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
	return eq_secs, newsecrefs

def label_order(label):
	""" Return order (distance from soma) based on label """
	if label.startswith('soma'):
		return 0
	elif label.startswith('trunk'):
		return 1
	elif label.startswith('smooth'):
		return 2
	elif label.startswith('spiny'):
		return 3
	else:
		return 4

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
	dendL_upper_root = getsecref(h.SThcell[0].dend0[1], dendLrefs)
	dendL_lower_root = getsecref(h.SThcell[0].dend0[2], dendLrefs)

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

	# Start iterative collapsing procedure
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

		############################################################################
		# 1. Merge a number of Y-sections
		logger.info("\n###############################################################"
					"\nFinding Y-sections to merge...\n")
		
		clusters = zip_fork_branches(allsecrefs, i_pass, zips_per_pass, Y_criterion='highest_level')

		############################################################################
		# 2. Create equivalent sections

		logger.info("\n###############################################################"
					"\nCreating equivalent sections...\n")
		# Mark Sections
		for secref in allsecrefs:
			secref.is_substituted = False
			secref.is_deleted = False

		eq_secs, newsecrefs = sub_equivalent_Y_sections(clusters, allsecrefs, path_props,
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

		# Print tree structure
		logger.info("Equivalent tree topology:")
		if logger.getEffectiveLevel() <= logging.DEBUG:
			h.topology()

	# Make sure all sections have identifiers
	assign_sth_indices(somaref, allsecrefs)
	return eq_secs, newsecrefs

def reduce_gillies_partial(delete_old_cells=True):
	"""
	Reduce Gillies & Willshaw STN neuron model by clustering and collapsing 
	its 3 identical trees independently.

	ALGORITHM
	- cluster each of 3 identical trees
	- merge/collapse segments in each cluster
	- create equivalent sections
		- conductances are interpolated according to distance L/lambda from soma
	"""
	############################################################################
	# 0. Load full model to be reduced (Gillies & Willshaw STN)
	for sec in h.allsec():
		if not sec.name().startswith('SThcell') and delete_old_cells:
			h.delete_section() # delete existing cells
	if not hasattr(h, 'SThcells'):
		h.xopen("gillies_cell_singleton.hoc")

	# Make sections accesible by name and index
	somaref = ExtSecRef(sec=h.SThcell[0].soma)
	dendLrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend0] # 0 is left tree
	dendRrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend1] # 1 is right tree
	allsecrefs = [somaref] + dendLrefs + dendRrefs

	# Get references to root sections of the 3 identical trees
	dendR_root = getsecref(h.SThcell[0].dend1[0], dendRrefs)
	dendL_juction = getsecref(h.SThcell[0].dend0[0], dendLrefs)
	dendL_upper_root = getsecref(h.SThcell[0].dend0[1], dendLrefs)
	dendL_lower_root = getsecref(h.SThcell[0].dend0[2], dendLrefs)

	# Assign indices used in Gillies code to read section properties from file
	somaref.tree_index = -1
	somaref.table_index = 0
	for j, dendlist in enumerate((dendLrefs, dendRrefs)):
		for i, secref in enumerate(dendlist):
			secref.tree_index = j # left tree is 0, right is 1
			secref.table_index = i+1 # same as in /sth-data/treeX-nom.dat

	############################################################################
	# 0. Pre-clustering: calculate properties

	# Assign lambda and electrotonic path length to each section and segment
	logger.info("Computing electrotonic path lengths...")
	f_lambda = 100. # frequency for electrotonic length constant lambda
	redtools.assign_electrotonic_length(somaref, allsecrefs, f_lambda, 
										gleak_name, allseg=True)

	############################################################################
	# 1. Cluster based on identified functional regions
	labels_by_order = ['soma', 'trunk', 'smooth', 'spiny']

	# Manually do soma cluster
	somaclu = Cluster('soma')
	somaclu.parent_label = 'soma'
	somaclu.parent_pos = 0.0
	somaclu.order = 0
	# Soma section
	somaref.cluster_labels = ['soma'] * somaref.sec.nseg
	# Junction sections
	dendL_juction.cluster_labels = ['soma'] * dendL_juction.sec.nseg
	soma_clu_refs = [somaref, dendL_juction]
	clusters = [somaclu]

	# Clustering parameters
	clu_fun = clutools.label_from_custom_regions
	clu_args = {'marker_mech': ('hh', 'gnabar')} # flag segments based on cluster label
	
	# Mark soma cluster manually
	for ref in soma_clu_refs:
		ref.sec.insert('hh')
		for seg in ref.sec:
			seg.gnabar_hh = 5

	# Cluster dendritic trees
	logger.info("Clustering left tree, upper...")
	clutools.clusterize_custom(dendL_upper_root, allsecrefs, clusters, '_0', clu_fun, clu_args)

	logger.info("Clustering left tree, lower...")
	clutools.clusterize_custom(dendL_lower_root, allsecrefs, clusters, '_1', clu_fun, clu_args)

	logger.info("Clustering right tree...")
	clutools.clusterize_custom(dendR_root, allsecrefs, clusters, '_2', clu_fun, clu_args)

	# Determine cluster relations/topology
	clutools.assign_topology(clusters, labels_by_order)
	# Manual edit
	clu_dendR_trunk = next((clu for clu in clusters if clu.label=='trunk_2'))
	clu_dendR_trunk.parent_pos = 0.0 # add to other end of soma

	############################################################################
	# 2. Create equivalent sections

	# Calculate equivalent properties by merging sections within each cluster
	logger.info("Merging within-cluster sections...")
	average_Ri = True
	for cluster in clusters:
		merge_seg_cluster(cluster, allsecrefs, average_Ri)

	# Create new sections
	logger.info("Creating equivalent sections...")
	eq_secs = redbush.equivalent_sections(clusters, allsecrefs, f_lambda, 
				gbar_scaling='area', interp_path=(1, (1,3,8)), interp_method='left_neighbor')

	############################################################################
	# 3. Finalize & Analyze
	
	# Delete original model sections & set ion styles
	for sec in h.allsec(): # makes each section the CAS
		if sec.name().startswith('SThcell') and delete_old_cells: # original model sections
			h.delete_section()
		else: # equivalent model sections
			h.ion_style("na_ion",1,2,1,0,1)
			h.ion_style("k_ion",1,2,1,0,1)
			h.ion_style("ca_ion",3,2,1,1,1)

	# Print tree structure
	logger.info("Equivalent tree topology:")
	if logger.getEffectiveLevel() <= logging.DEBUG:
		h.topology()

	return clusters, eq_secs

def reduce_gillies_pathRi(customclustering, average_Ri):
	""" Reduce Gillies & Willshaw STN neuron model

	To set active conductances, interpolates using axial path resistance
	values (sum(Ri)).

	@param customclustering		see param 'customclustering' in function
								cluster_sections()

	@param average_Ri			see param 'average_Ri' in function
								merge_sec_cluster
	"""

	# Initialize Gillies model
	for sec in h.allsec():
		if not sec.name().startswith('SThcell'): # original model sections
			h.delete_section()
	if not hasattr(h, 'SThcells'):
		h.xopen("gillies_cell_singleton.hoc")

	# Make sections accesible by both name and index + allow to add attributes
	somaref = ExtSecRef(sec=h.SThcell[0].soma)
	dendLrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend0]
	dendRrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend1]
	alldendrefs = dendLrefs + dendRrefs
	allsecrefs = [somaref] + alldendrefs
	or_refs = (somaref, dendLrefs, dendRrefs)

	############################################################################
	# 0. Pre-clustering: calculate properties

	# Assign Strahler numbers
	logger.info("Assingling Strahler's numbers...")
	clutools.assign_topology_attrs(dendLrefs[0], dendLrefs)
	clutools.assign_topology_attrs(dendRrefs[0], dendRrefs)
	somaref.order = 0 # distance from soma
	somaref.level = 0
	somaref.strahlernumber = dendLrefs[0].strahlernumber # same as root of left tree

	############################################################################
	# 1. Cluster based on identified functional regions

	# Cluster sections
	logger.info("Clustering sections...")
	dendLroot, dendRroot = dendLrefs[0], dendRrefs[0]

	# Cluster soma
	somaref.cluster_label = 'soma'
	somaclu = Cluster('soma')
	somaclu.parent_label = 'soma'
	somaclu.parent_pos = 0.0
	somaclu.order = 0
	clusters = [somaclu]

	# Cluster dendritic trees
	if customclustering:
		# Based on diameters: see Gilles & Willshaw (2006) fig. 1
		clu_fun = lambda (ref): 'spiny' if ref.sec.diam <= 1.0 else 'trunk'
		clu_args = {}
	else:
		clu_fun = clutools.label_from_strahler
		clu_args = {'thresholds':(1,2)}

	clutools.clusterize_custom(dendLroot, allsecrefs, clusters, '_0', clu_fun, clu_args)
	clutools.clusterize_custom(dendRroot, allsecrefs, clusters, '_1', clu_fun, clu_args)

	# Determine cluster relations/topology
	clutools.assign_topology(clusters, ['soma', 'trunk', 'smooth', 'spiny'])

	# Debug info
	for cluster in clusters:
		logger.debug("Cluster '{0}'' has parent cluster '{1}'".format(cluster.label, cluster.parent_label))

	############################################################################
	# 2. Create equivalent sections (compute ppties from cluster sections)

	# Calculate cluster properties
	logger.info("Merging within-cluster sections...")
	for cluster in clusters:
		# Merge sections within each cluster: 
		merge_sec_cluster(cluster, allsecrefs, average_Ri)

		# Map axial path resistance values to gbar values
		map_pathri_gbar(cluster, allsecrefs)

	# Create equivalent section for each cluster
	logger.info("Creating equivalent sections...")
	eq_secs, eq_secrefs = equivalent_sections(clusters, allsecrefs, gradients=True)

	# Sort equivalent sections
	eq_sorted = sorted(eq_secrefs, key=lambda ref: label_order(ref.cluster_label))
	eq_somaref = next(ref for ref in eq_secrefs if ref.cluster_label.startswith('soma'))
	eq_dendLrefs = [ref for ref in eq_secrefs if ref.cluster_label.endswith('0')]
	eq_dendRrefs = [ref for ref in eq_secrefs if ref.cluster_label.endswith('1')]
	eq_somasec = eq_somaref.sec
	eq_dendLsecs = [ref.sec for ref in eq_dendLrefs]
	eq_dendRsecs = [ref.sec for ref in eq_dendRrefs]

	############################################################################
	# 3. Finalize & Analyze

	# Compare full/reduced model
	eq_secs = (eq_somasec, eq_dendLsecs, eq_dendRsecs)
	eq_refs = (eq_somaref, eq_dendLrefs, eq_dendRrefs)
	analysis.compare_models(or_refs, eq_refs, [])

	# Delete original model sections
	for sec in h.allsec(): # makes each section the CAS
		if sec.name().startswith('SThcell'): # original model sections
			h.delete_section()
		else: # equivalent model sections
			h.ion_style("na_ion",1,2,1,0,1)
			h.ion_style("k_ion",1,2,1,0,1)
			h.ion_style("ca_ion",3,2,1,1,1)

	logger.info("Equivalent tree topology:")
	if logger.getEffectiveLevel() <= logging.DEBUG:
		h.topology() # prints topology
	
	# return data structures
	return clusters, eq_secs, eq_refs

if __name__ == '__main__':
	# clusters, eq_secs = reduce_gillies_partial(delete_old_cells=True)
	eq_secs, newsecrefs = reduce_gillies_incremental(n_passes=7, zips_per_pass=100)
	from neuron import gui # check in ModelView: conductance distribution, structure