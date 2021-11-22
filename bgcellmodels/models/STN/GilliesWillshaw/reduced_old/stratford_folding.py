"""
Reduce model by folding/collapsing branches according to algorithm described in following articles:

- Stratford, K., Mason, A., Larkman, A., Major, G., and Jack, J. J. B. (1989) - The modelling of pyramidal neurones in the visual cortex

- Fleshman, J. W., Segev, I., and Burke, R. E. (1988) - Electrotonic architecture of type-identified alpha-motoneurons in the cat spinal cord. J. Neurophysiol. 60: 60 85.

- Clements, J., and Redman, S. (1989) - Cable properties of cat spinal motoneurones measured by combining voltage clamp, current clamp and intracellular staining. J. Physiol. 409: 63 87.

- Algorithm overview in 'Methods in Neural Modeling', Chapter 3.4, equations 3.21 - 3.22

@author Lucas Koelman

@date	28-08-2017
"""

import re
import math

import redutils
import interpolation as interp

from redutils import get_sec_props_obj
from cluster import Cluster, assign_topology_attrs

from bgcellmodels.common import treeutils
from bgcellmodels.common.treeutils import getsecref, seg_index, interp_seg, next_segs
from bgcellmodels.common.electrotonic import calc_lambda, measure_transfer_impedance as measure_Zin, segs_at_dX

from neuron import h

# Global parameters
f_lambda = 0.0
alphabet_uppercase = [chr(i) for i in range(65,90+1)] # A-Z are ASCII 65-90

# logging of DEBUG/INFO/WARNING messages
import logging
logging.basicConfig(format='%(message)s %(levelname)s:@%(filename)s:%(lineno)s', level=logging.DEBUG)
logname = "stratford" # __name__
logger = logging.getLogger(logname) # create logger for this module

################################################################################
# Folding Algorithm
################################################################################

maxlvl = 0

def next_eq_diam(a_seg, a_i, root_X, dX, cluster, reduction, lvl=0):
	"""
	Compute equivalent diameter at next discretization step (a_i+1) * dX
	in the Stratford et al. (1989) reduction method.

	@param	a_seg		nrn.Segment : last segment encountered when ascending the tree.

	@param	a_i			discretization step index of last computed diameter
						(use 0 for first call)

	@param	root_X		electrotonic distance of the first segment (the subtree
						root segment) from soma or root of whole tree

	@param	dX			discretization step size
	"""
	allsecrefs = reduction.all_sec_refs
	gleak_name = reduction.gleak_name
	gbar_list = reduction.active_gbar_names

	# electrotonic distance to last segment where diam was calculated
	# last_X = root_X + a_i*dX
	# target_X = last_X + dX

	# electrotonic distance to last visited segment (not necessarily used)
	# a_ref		= getsecref(a_seg.sec, allsecrefs)
	# a_index		= seg_index(a_seg)
	# a_lambda	= a_ref.seg_lambda[a_index]
	# a_X0		= a_ref.seg_path_Lelec0[a_index]
	# a_X1		= a_ref.seg_path_Lelec1[a_index]
	# a_X			= interp_seg(a_seg, a_X0, a_X1)
	# assert target_X > a_X # should be guaranteed by last call
	
	pref = (lvl*"-") # prefix for printing
	global maxlvl
	if lvl > maxlvl:
		maxlvl += 1
		brk = ""
	else:
		maxlvl = lvl
		brk = "\n"

	# Get child segments
	child_segs = segs_at_dX(a_seg, dX, f_lambda, reduction.gleak_name) # segments/nodes connected to this one
	if not any(child_segs):
		if (a_seg.x == 1.0):
			logger.debug(pref+"Folding: reached end of branch @ {} (NOT USED)".format(a_seg))
			return # reached end of branch
		else:
			child_segs = [a_seg.sec(1.0)] # test end of branch
	
	# Get new segment at (a_i+1) * dX
	for iseg, b_seg in enumerate(child_segs):

		# NOTE: depth-first tree traversal: look for next point ta calculate diam
		#       (next discretization step) along current branch, then ascend further
		if len(child_segs) > 1:
			logger.debug(brk+pref+"Branched at {}".format(a_seg))
		logger.debug(pref+"Child [{}/{}] = {}".format(iseg+1, len(child_segs), b_seg))

		# Get electronic distance to following/adjacent segment
		# b_ref = getsecref(b_seg.sec, allsecrefs)
		# b_index = seg_index(b_seg)
		# b_X0 = b_ref.seg_path_Lelec0[b_index]
		# b_X1 = b_ref.seg_path_Lelec1[b_index]
		# b_X = interp_seg(b_seg, b_X0, b_X1)
		# last_dX = b_X - last_X

		# # If distance to last used segment is smaller than step size: ascend to next segment
		# if target_X > b_X: # (if dX > last_dX)
			
		# 	# Examples: (s=last_seg, a=a_seg, b=b_seg, b=target_seg)
		# 	#	[-----|-s-a-|b--t-] - [t----|-----|-----]
		# 	#	[-----|---sa|--b--] - [t----|-----|-----]
		# 	print(pref+"Skip child {}".format(b_seg))
		# 	# print((a_i*"-") + "Skipping {} since bX={} < tarX={}".format(b_seg, b_X, target_X))

		# 	# Ascend but don't increase index
		# 	next_eq_diam(b_seg, a_i, root_X, dX, cluster, reduction, lvl=lvl+1)

		# else: # target_X <= b_X
		# 	# target X value is between previous and next segment -> get new segment by interpolation
			
		# 	# Examples: (s=last_seg, a=a_seg, b=b_seg, b=target_seg)
		# 	#	[-----|-s-a-|t--b-] - [-----|-----|-----]
		# 	#	[-----|---s-|a-t--] - [b----|-----|-----]
		# 	print(pref+"Interp child {}".format(b_seg))
		# 	# print((a_i*"-") + "Interpolating between {} (X={}) and {} (X={})".format(a_seg, a_X, b_seg, b_X))
			
		# 	fX = (target_X - a_X) / (b_X - a_X)
		# 	new_seg = None
		
		# 	if b_seg.sec.same(a_seg.sec): # segment A and B are in same cylinder
				
		# 		# Examples: (s=last_seg, a=a_seg, b=b_seg, b=target_seg)
		# 		#	[-----|-s-a-|t--b-] - [t----|-----|-----]
				
		# 		# linear interpolation
		# 		assert b_seg.x > a_seg.x
		# 		d_x = b_seg.x - a_seg.x
		# 		new_x = a_seg.x + (fX * d_x)
		# 		new_seg = a_seg.sec(new_x)

		# 	else: # segment B is in next cylinder (Section)

		# 		# Examples: (s=last_seg, a=a_seg, b=b_seg, b=target_seg)
		# 		#	[-----|---s-|a-t--] - [b----|-----|-----]

		# 		if fX <= 0.5: # new segment is in cylinder A

		# 			d_x = 1.0 - a_seg.x
		# 			frac = fX/0.5
		# 			new_x = a_seg.x + (frac * d_x)
		# 			new_seg = a_seg.sec(new_x)

		# 		else: # new segment is in cylinder B

		# 			frac = (fX-0.5) / 0.5
		# 			new_x = frac*b_seg.x
		# 			new_seg = b_seg.sec(new_x)

		new_seg = b_seg

		# Save attributes at this electrotonic distance from root
		step_attrs = {}
		step_attrs['step_diams'] = [new_seg.diam]
		step_attrs['step_Rm'] = [1.0 / getattr(new_seg, gleak_name)]
		step_attrs['step_Ra'] = [new_seg.sec.Ra]
		step_attrs['step_cm'] = [new_seg.cm]
		step_attrs['step_gbar'] = [dict(((gname, getattr(new_seg, gname)) for gname in gbar_list))]

		# Save them on cluster object
		for attr_name, new_seg_vals in step_attrs.items():

			# Get cluster attribute
			step_attr_list = getattr(cluster, attr_name) # list(list()) : [[step1_d1...], [step2_d1...], ...]
			
			# Save values for new segment
			if len(step_attr_list) <= a_i:
				# First entry (first branch for this step) : start new list
				step_attr_list.append(new_seg_vals) # add new list to the [[step1...], [step2...], ...]
			else:
				# Subsequent branches: add to list for this step
				step_parallel_list = step_attr_list[a_i]
				step_parallel_list.extend(new_seg_vals)


		# Ascend and increase index
		next_eq_diam(new_seg, a_i+1, root_X, dX, cluster, reduction, lvl=lvl+1)

	return


def next_seg_diam(a_seg, a_i, cluster, reduction, lvl=0):
	"""
	Compute equivalent diameter at next parallel segments up the tree, assuming
	they have the same electrotonic distance from the current segment.
	
	This would be approximately true if the segment width is chosen as a fraction
	of lambda (but will unlikely be exactly true due to the discrete nseg).

	@param	a_seg		nrn.Segment : last segment encountered when ascending the tree.

	@param	a_i			discretization step index of last computed diameter
						(use 0 for first call)
	"""
	allsecrefs = reduction.all_sec_refs
	gleak_name = reduction.gleak_name
	gbar_list = reduction.active_gbar_names

	# electrotonic distance to last visited segment (not necessarily used)
	a_ref		= getsecref(a_seg.sec, reduction.all_sec_refs)
	a_index		= seg_index(a_seg)
	a_X0		= a_ref.seg_path_Lelec0[a_index]
	a_X1		= a_ref.seg_path_Lelec1[a_index]
	a_X			= interp_seg(a_seg, a_X0, a_X1)

	
	# Get child segments
	child_segs = next_segs(a_seg)
	if not any(child_segs):
		logger.debug("Folding: reached end of branch @ {} (NOT USED)".format(a_seg))
	
	# Keep track of electrotonic distance to children
	cur_dX = None
	tol = 0.01
	within_tol = lambda a,b: (1.0 - tol) <= (a / b) <= (1.0 + tol)
	
	# Get new segment at (a_i+1) * dX
	for iseg, b_seg in enumerate(child_segs):

		new_seg = b_seg

		# Get electronic distance to following/adjacent segment
		b_ref = getsecref(b_seg.sec, allsecrefs)
		b_index = seg_index(b_seg)
		b_X0 = b_ref.seg_path_Lelec0[b_index]
		b_X1 = b_ref.seg_path_Lelec1[b_index]
		b_X = interp_seg(b_seg, b_X0, b_X1)
		new_dX = b_X - a_X

		# Check if precondition for function is still true
		if (cur_dX is not None) and (not within_tol(new_dX, cur_dX)):
			# FIXME: this will only check at branch point, also test for deeper parallel branches: basically all step[i] must occur at the same X
			raise ValueError("Dendritic tree did not meet precondition that "
					"segments must be at equal electronic distance from branch points")

		# Save attributes at this electrotonic distance from root
		step_attrs = {}
		step_attrs['step_diams'] = [new_seg.diam]
		step_attrs['step_Rm'] = [1.0 / getattr(new_seg, gleak_name)]
		step_attrs['step_Ra'] = [new_seg.sec.Ra]
		step_attrs['step_cm'] = [new_seg.cm]
		step_attrs['step_gbar'] = [dict(((gname, getattr(new_seg, gname)) for gname in gbar_list))]

		# Save them on cluster object
		for attr_name, new_seg_vals in step_attrs.items():

			# Get cluster attribute
			step_attr_list = getattr(cluster, attr_name) # list(list()) : [[step1_d1...], [step2_d1...], ...]
			
			# Save values for new segment
			if len(step_attr_list) <= a_i:
				# First entry (first branch for this step) : start new list
				step_attr_list.append(new_seg_vals) # add new list to the [[step1...], [step2...], ...]
			else:
				# Subsequent branches: add to list for this step
				step_parallel_list = step_attr_list[a_i]
				step_parallel_list.extend(new_seg_vals)


		# Ascend and increase index
		next_seg_diam(new_seg, a_i+1, cluster, reduction, lvl=lvl+1)

	return


def calc_folds(target_Y_secs, i_pass, reduction, dX=0.1):
	"""
	Do folding operation

	@return			list of Cluster objects with properties of equivalent
					Section for each set of collapsed branches.
	"""

	gbar_list = reduction.active_gbar_names

	clusters = []
	for j_zip, root_ref in enumerate(target_Y_secs):

		# Fold name
		p = re.compile(r'(?P<cellname>\w+)\[(?P<cellid>\d+)\]\.(?P<seclist>\w+)\[(?P<secid>\d+)\]')
		pmatch = p.search(root_ref.sec.name())
		if pmatch:
			pdict = pmatch.groupdict()
			name_sanitized = "{}_{}".format(pdict['seclist'], pdict['secid'])
		else:
			name_sanitized = re.sub(r"[\[\]\.]", "", root_ref.sec.name())
		fold_label = "fold_{0}".format(name_sanitized)

		# Make Cluster object that represents collapsed segments
		cluster = Cluster(fold_label)

		# Branch parameters
		cluster.step_diams = []
		cluster.step_Rm = []
		cluster.step_Ra = []
		cluster.step_cm = []
		cluster.step_gbar = []

		# Topological info
		cluster.parent_seg = root_ref.sec(1.0)
		clu_secs = treeutils.subtree_secs(root_ref.sec)
		clu_segs = [seg for sec in clu_secs for seg in sec]

		# Geometrical and electrical info
		cluster.or_area = sum(seg.area() for seg in clu_segs)
		cluster.or_cmtot = sum(seg.cm*seg.area() for seg in clu_segs)
		cluster.or_cm = cluster.or_cmtot / cluster.or_area
		cluster.or_gtot = dict((gname, 0.0) for gname in reduction.gbar_names)
		
		for gname in reduction.gbar_names:
			cluster.or_gtot[gname] += sum(getattr(seg, gname)*seg.area() for seg in clu_segs)

		# Ascend subtree from root and compute diameters
		start_seg = root_ref.sec(1.0) # end of cylinder
		start_X = root_ref.pathLelec1
		logger.debug("Folding: start tree ascent @ {}".format(start_seg))
		next_eq_diam(start_seg, 0, start_X, dX, cluster, reduction)

		# Equivalent parameters
		cluster.num_sec = len(cluster.step_diams)
		cluster.eq_diam =	[0.0] * cluster.num_sec
		cluster.eq_L =		[0.0] * cluster.num_sec
		cluster.eq_Rm =		[0.0] * cluster.num_sec
		cluster.eq_Ra =		[0.0] * cluster.num_sec
		cluster.eq_cm =		[0.0] * cluster.num_sec
		cluster.eq_lambda =	[0.0] * cluster.num_sec
		cluster.eq_gbar =	[0.0] * cluster.num_sec

		# Finalize diam calculation
		for i_step, diams in enumerate(cluster.step_diams):

			# Number of parallel branches found at this step
			num_parallel = len(diams)
			logger.anal("Combining diameters of {} parallel branches in discretization step {}".format(num_parallel, i_step))
			assert len(cluster.step_diams[i_step]) == len(cluster.step_Rm[i_step]) == len(cluster.step_Ra[i_step]) == len(cluster.step_cm[i_step])

			# Equivalent properties at this step
			cluster.eq_diam[i_step] = sum((d**(3.0/2.0) for d in diams)) ** (2.0/3.0)

			# Passive electrical properties
			cluster.eq_Rm[i_step] = sum(cluster.step_Rm[i_step]) / num_parallel
			cluster.eq_Ra[i_step] = sum(cluster.step_Ra[i_step]) / num_parallel
			cluster.eq_cm[i_step] = sum(cluster.step_cm[i_step]) / num_parallel

			# Combine dict of gbar for each parallel section
			cluster.eq_gbar[i_step] = dict(((gname, 0.0) for gname in gbar_list))
			for gdict in cluster.step_gbar[i_step]:
				for gname, gval in gdict.items():
					cluster.eq_gbar[i_step][gname] += (gval / num_parallel) # average of parallel branches

			# Physical length 
			# NOTE: it is important to use the same lambda equation as the one
			#       used for electrotonic path lengths
			cluster.eq_lambda[i_step] = calc_lambda(f_lambda, cluster.eq_diam[i_step], cluster.eq_Ra[i_step], 
													1.0/cluster.eq_Rm[i_step], cluster.eq_cm[i_step])
			
			cluster.eq_L[i_step] = dX * cluster.eq_lambda[i_step]

		# Save cluster
		clusters.append(cluster)

	return clusters


def new_fold_section(cluster, i_step, j_sec, reduction):
	"""
	Start new Section for discretization step i of the cluster.
	"""
	sec, ref = treeutils.create_hoc_section("{}_{}".format(
								cluster.label, alphabet_uppercase[j_sec]))
	
	sec.L		= cluster.eq_L[i_step]
	sec.nseg	= 1
	sec.diam	= cluster.eq_diam[i_step]
	sec.Ra		= cluster.eq_Ra[i_step]
	sec.cm		= cluster.eq_cm[i_step]

	# Insert mechanisms
	for mech in reduction.mechs_gbars_dict:
		sec.insert(mech)

	# Set conductances
	setattr(sec(1.0), reduction.gleak_name, 1.0/cluster.eq_Rm[i_step])
	for gname in reduction.active_gbar_names:
		setattr(sec(1.0), gname, cluster.eq_gbar[i_step][gname])

	return sec, ref


def extend_fold_section(cluster, i_step, sec, reduction):
	"""
	Extend existing fold Section with a new segment for given discretization step.
	"""

	sec.L			= sec.L + cluster.eq_L[i_step]
	sec.nseg		= sec.nseg + 1
	sec(1.0).diam	= cluster.eq_diam[i_step]
	sec(1.0).cm		= cluster.eq_cm[i_step]
	# NOTE: Ra is not RANGE var -> only set once per Section

	# Set conductances
	setattr(sec(1.0), reduction.gleak_name, 1.0/cluster.eq_Rm[i_step])
	for gname in reduction.active_gbar_names:
		setattr(sec(1.0), gname, cluster.eq_gbar[i_step][gname])



def make_substitute_folds(clusters, reduction):
	"""
	Make equivalent Sections for each cluster, and substitute them into cell
	to replace the folded compartments.
	"""

	# Create one Section for each region with uniform diam
	# (uniform diam also means uniform length: see equation for physical length)
	diam_tolerance = 0.05 # 5 percent

	for i_clu, cluster in enumerate(clusters):

		# Save equivalent sections
		cluster.eq_refs = []

		# Create first Section for this fold
		cur_diam = cluster.eq_diam[0]
		cur_sec, cur_ref = new_fold_section(cluster, 0, 0, reduction)
		cur_sec.connect(cluster.parent_seg, 0.0) # for electrotonic properties
		
		num_fold_secs = 1
		cluster.eq_refs.append(cur_ref)

		# Create Section for each region with uniform diam
		for i_step, diam in enumerate(cluster.eq_diam):

			if (1.0 - diam_tolerance) <= (diam / cur_diam) <= (1.0 + diam_tolerance):
				
				# Uniform: extend current Section with new segment
				extend_fold_section(cluster, i_step, cur_sec, reduction)

			else:
				
				# Not uniform: start new Section
				new_sec, new_ref = new_fold_section(cluster, i_step, num_fold_secs, reduction)
				cluster.eq_refs.append(new_ref)

				# Connect to previous sections
				new_sec.connect(cur_sec(1.0), 0.0)

				# Update current value
				num_fold_secs += 1
				cur_diam = diam
				cur_sec, cur_ref = new_sec, new_ref

		# Calculate new total area
		cluster.eq_area = sum(seg.area() for ref in cluster.eq_refs for seg in ref.sec)
		logger.debug("Cluster @ {} : original area = {} , new area = {} [um^2]".format(
						cluster.parent_seg, cluster.or_area, cluster.eq_area))

		# Set nonuniform active conductances
		# TODO: plot response and test different dX
		for eq_ref in cluster.eq_refs:

			# first calculate electrotonic properties of equivalent sections
			redutils.sec_path_props(eq_ref, f_lambda, reduction.gleak_name)

			# interpolate conductance in each segment
			for j_seg, seg in enumerate(eq_ref.sec):
			
				# Get adjacent segments along path
				seg_X0, seg_X1 = eq_ref.seg_path_Lelec0[j_seg], eq_ref.seg_path_Lelec1[j_seg]
				seg_X = interp_seg(seg, seg_X0, seg_X1)
				bound_segs, bound_X = interp.find_adj_path_segs('path_L_elec', seg_X, reduction.interp_path)

				# Set conductances by interpolating neighbors
				for gname in reduction.active_gbar_names:
					gval = interp.interp_gbar_linear_neighbors(seg_X, gname, bound_segs, bound_X)
					seg.__setattr__(gname, gval)

			# Check if conductances need scaling
			gbar_scaling = reduction.get_reduction_param(reduction.active_method, 'gbar_scaling')
			if (gbar_scaling is None) or (gbar_scaling == 'none'):
				continue

			# Re-scale gbar distribution to yield same total gbar (sum(gbar*area))
			for gname in reduction.gbar_names:
				
				# Get original and current total conductance
				eq_gtot = sum(getattr(seg, gname)*seg.area() for seg in eq_ref.sec)
				if eq_gtot <= 0.:
					eq_gtot = 1.
				
				or_gtot = cluster.or_gtot[gname]
				
				# Scale each gbar
				for j_seg, seg in enumerate(eq_ref.sec):
					
					if gbar_scaling == 'area':
						# conserves ratio in each segment but not total original conductance
						scale = cluster.or_area/cluster.eq_area
					
					elif gbar_scaling == 'gbar_integral':
						# does not conserve ratio but conserves gtot_or since: sum(g_i*area_i * or_area/eq_area) = or_area/eq_area * sum(gi*area_i) ~= or_area/eq_area * g_avg*eq_area = or_area*g_avg
						scale = or_gtot/eq_gtot
					
					else:
						raise ValueError("Unknown gbar scaling method '{}'.".format(gbar_scaling))
					
					# Set gbar
					gval = getattr(seg, gname) * scale
					seg.__setattr__(gname, gval)


		# Disconnect substituted Sections
		clu_root_sec = cluster.eq_refs[0].sec
		redutils.sub_equivalent_Y_sec(clu_root_sec, cluster.parent_seg, [], 
								reduction.all_sec_refs, reduction.mechs_gbars_dict, 
								delete_substituted=True)

		# Set ion styles
		for ref in cluster.eq_refs:
			redutils.set_ion_styles(ref.sec, **reduction._ion_styles)



################################################################################
# Interface Functions (FoldReduction)
################################################################################


def preprocess_impl(reduction):
	"""
	Preprocess cell for Stratford reduction.

	(Implementation of interface declared in reduce_cell.CollapseReduction)

	@param	reduction		reduce_cell.CollapseReduction object
	"""
	
	dendL_secs = list(h.SThcell[0].dend0)
	dendR_secs = list(h.SThcell[0].dend1)

	allsecrefs = reduction.all_sec_refs
	dend_lists = [dendL_secs, dendR_secs]

	# Assign indices used in Gillies code (sth-data folder)
	for somaref in reduction.soma_refs:
		somaref.tree_index = -1
		somaref.table_index = 0

	for secref in reduction.dend_refs:
		for i_dend, dendlist in enumerate(dend_lists):
			if any([sec.same(secref.sec) for sec in dendlist]):
				secref.tree_index = i_dend
				secref.table_index= next((i+1 for i,sec in enumerate(dendlist) if sec.same(secref.sec)))

	# Choose representative path for interpolation
	interp_tree_id = 1
	interp_table_ids = (1,3,8) # from soma to end of SThcell[0].dend1[7]
	reduction.interp_path = []

	# Save properties along this path
	for ref in allsecrefs:
		if (ref.tree_index==interp_tree_id) and (ref.table_index in interp_table_ids):
			
			# calculate properties
			redutils.sec_path_props(ref, f_lambda, reduction.gleak_name)

			# Save properties
			reduction.interp_path.append(get_sec_props_obj(ref, reduction.mechs_gbars_dict,
										['pathL_elec'], ['pathLelec0', 'pathLelec1']))
	
	# Compute input resistance at soma
	method				= reduction.active_method
	somaref				= reduction.soma_refs[0]
	Z_freq				= reduction.get_reduction_param(method, 'Z_freq')
	init_func			= reduction.get_reduction_param(method, 'Z_init_func')
	linearize_gating	= reduction.get_reduction_param(method, 'Z_linearize_gating')
	init_func()
	Zin_DC = measure_Zin(somaref.sec(0.5), freq=0.0, linearize_gating=linearize_gating)
	Zin_AC = measure_Zin(somaref.sec(0.5), freq=Z_freq, linearize_gating=linearize_gating)
	somaref.orig_Zin_DC = Zin_DC
	somaref.orig_Zin_AC = Zin_AC
		


def prepare_folds_impl(reduction):
	"""
	Prepare next collapse operation: assign topology information
	to each Section.

	(Implementation of interface declared in reduce_cell.CollapseReduction)

	NOTE: topology is ONLY used for when using multiple folding passes, 
	      for determining the folding branch point
	"""
	all_refs = reduction.all_sec_refs
	root_ref = reduction._soma_refs[0]

	# Assign topology info (order, level, strahler number)
	root_ref.order = 0
	root_ref.level = 0
	assign_topology_attrs(root_ref, all_refs, root_order=0)
	root_ref.strahlernumber = max((ref.strahlernumber for ref in all_refs))
	
	# For each segment: compute path length, path resistance, electrotonic path length
	for secref in reduction.all_sec_refs:
		redutils.sec_path_props(secref, f_lambda, reduction.gleak_name)
		secref.max_passes = 100


def calc_folds_impl(reduction, i_pass):
	"""
	Collapse branches at branch points identified by given criterion.
	"""
	# Get sections
	allsecrefs = reduction.all_sec_refs

	# Find collapsable branch points
	# target_Y_secs = folding.find_collapsable(allsecrefs, i_pass, 'highest_level')
	target_Y_secs = reduction._fold_root_refs

	# Do collapse operation at each branch points
	dX = reduction.get_reduction_param(reduction.active_method, 'fold_dX')
	clusters = calc_folds(target_Y_secs, i_pass, reduction, dX=dX)

	# Save results
	reduction.clusters = clusters


def make_folds_impl(reduction):
	"""
	Make equivalent Sections for branches that have been folded.
	"""
	clusters = reduction.clusters

	# Make new sections
	make_substitute_folds(clusters, reduction) # SectionRef stored on each Cluster

	# Update Sections
	new_dend_refs = sum((cluster.eq_refs for cluster in clusters), []) # concatenate lists
	reduction.update_refs(dend_refs=new_dend_refs) # prepare for next iteration


def post_process_impl(reduction):
	"""
	Post-process cell for Stratford reduction.

	(interface declared in reduce_cell.CollapseReduction)

	@pre		somatic section must consist of single compartment (segment)

	@effect		Measure DC and AC input impedance at soma, and compensate discrepancy
				with original cell model by adding a somatic resistive shunt 
				(decreasing membrane resistance, i.e. increasing leak conductance)
				and a somatic capacitive shunt (increasing membrane capacitance).

	@post		segment.gleak and segment.cm of middle somatic segment will be
				modified
	"""
	somaref = reduction.soma_refs[0]
	somasec = somaref.sec
	somaseg = somasec(0.5)
	assert somasec.nseg == 1, "Method only works for soma sections consisting of one segment."

	# Get parameters
	red_method = reduction.active_method
	Z_freq = reduction.get_reduction_param(red_method, 'Z_freq')
	init_cell = reduction.get_reduction_param(red_method, 'Z_init_func')
	linearize_gating = reduction.get_reduction_param(red_method, 'Z_linearize_gating')
	
	# Compute transfer impedances
	init_cell()
	imp = h.Impedance() # imp = h.zz # previously created
	imp.loc(somaseg.x, sec=somasec) # injection site

	# DC Measurement ###########################################################
	
	# Compute input impedance at soma
	imp.compute(0.0, int(linearize_gating)) 
	red_Zin_DC = imp.input(somaseg.x, sec=somasec)

	# Get discrepancy between Zin (Gin), and compensate using somatic shunt 
	# (gradient descent until Gin matches?)
	old_Zin_DC = somaref.orig_Zin_DC
	diff_Gin_DC = (1.0 / old_Zin_DC) - (1.0 / red_Zin_DC) # red_Gin should be lower -> diff_Gin > 0

	# Gtot = Glocal + Gend
	# we can change Glocal = jwC + 1/Rm
	# first correct Rm using DC measurement
	# put an extra gm in parallel to make up the difference
	gm_old = getattr(somaseg, reduction.gleak_name)
	gm_shunt = diff_Gin_DC / somaseg.area()
	gm_new = gm_old + gm_shunt
	setattr(somaseg, reduction.gleak_name, gm_new)

	logger.debug('[DC:0Hz] Original Zin = {}, Reduced model has Zin = {}'.format(old_Zin_DC, red_Zin_DC))
	logger.debug('Added shunt conductance of {} [S/cm2] to segment {}'.format(gm_shunt, somaseg))
	logger.debug('Changed gleak from {} -> {} [S/cm2]'.format(gm_old, gm_new))

	# Check result
	imp.compute(0.0, int(linearize_gating)) 
	red_Zin_DC_new = imp.input(somaseg.x, sec=somasec)
	tolerance = 0.05
	
	if not (1.0-tolerance <= (red_Zin_DC_new / old_Zin_DC) <= 1.0+tolerance):
		logger.warning("Zin_DC not fixed within tolerance: "
				"Zin_old = {}, Zin_new = {}".format(old_Zin_DC, red_Zin_DC_new))
	
	logger.debug('After resistive shunt placement: Zin = {}'.format(red_Zin_DC_new))

	# AC Measurement ###########################################################

	if Z_freq == 0.0:
		return # same as previous calculation

	# Compute AC input impedance at soma
	imp.compute(Z_freq, int(linearize_gating)) # compute A(x->loc) for all x where A is Vratio/Zin/Ztransfer
	red_Zin_AC = imp.input(somaseg.x, sec=somasec)
	
	# Get discrepancy between Zin (Gin), and compensate using capactive shunt
	old_Zin_AC = somaref.orig_Zin_AC
	ratio_Zin_AC = red_Zin_AC / old_Zin_AC
	diff_Gin_AC = (1.0/old_Zin_AC) - (1.0/red_Zin_AC)

	# Test if shunt is needed
	logger.debug('[AC:{}Hz] Original Zin = {}, Reduced model has Zin = {}'.format(Z_freq, old_Zin_AC, red_Zin_AC))
	if (1.0-tolerance <= ratio_Zin_AC <= 1.0+tolerance):
		logger.debug("Zin_AC is already within tolerance. Skip scaling Cm.")
		return

	# Gtot = Glocal + Gend
	# we will only change jwC part of Glocal = jwC + 1/Rm
	cm_old = somaseg.cm
	cm_shunt = diff_Gin_AC / (somaseg.area() * 2.0 * math.pi * Z_freq)
	cm_new = cm_old + cm_shunt
	somaseg.cm = cm_new

	logger.debug('Added shunt capacitance of {} [uf/cm2] to segment {}'.format(gm_shunt, somaseg))
	logger.debug('Changed cm from {} -> {} [uf/cm2]'.format(cm_old, cm_new))

	# Check result
	imp.compute(Z_freq, int(linearize_gating)) 
	red_Zin_AC_new = imp.input(somaseg.x, sec=somasec)
	
	if not (1.0-tolerance <= (red_Zin_AC_new / old_Zin_AC) <= 1.0+tolerance):
		logger.warning("Zin_AC not fixed within tolerance: "
				"Zin_old = {}, Zin_new = {}".format(old_Zin_AC, red_Zin_AC_new))

	logger.debug('After capacitive shunt placement: Zin = {}'.format(red_Zin_AC_new))
