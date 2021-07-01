"""
Marasco (2012) reduction, cluster-based approach.

Cluster-based means: sections are partitioned into clusters based on some cirterion,
then the equivalent section for each cluster is computed. Sections that are clustered
together are not necessarily adjoining / contiguous.


@author		Lucas Koelman
"""

# Python modules
import math
PI = math.pi

# NEURON modules
import neuron
h = neuron.h
h.load_file("stdlib.hoc") # Load the standard library

import cluster as clutools
from cluster import Cluster
import reduce_bush_sejnowski as redbush
from marasco_merging import merge_seg_subtree, merge_sec_sequential
import analyze_reduction as analysis

# Gillies STN model
from gillies_model import gillies_gdict, gillies_mechs, gillies_glist
mechs_chans = gillies_gdict
glist = gillies_glist
gleak_name = 'gpas_STh'
f_lambda = 100.0

# logging of DEBUG/INFO/WARNING messages
import logging
logging.basicConfig(format='%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logname = "reduction" # __name__
logger = logging.getLogger(logname) # create logger for this module
fmtr = logging.Formatter('%(levelname)s:%(message)s @%(filename)s:%(lineno)s')


################################################################################
# Merging Clusters
################################################################################

def merge_sec_cluster(cluster, allsecrefs, average_Ri):
	"""
	Calculate cluster equivalent section properties using section-based merging.

	@param average_Ri	If True, Ra is calculated so that Ri (absolute axial
						resistance) of the equivalent section for each cluster is 
						the average Ri of all disconnected subtrees merged into that
						section. This does NOT conserve input resistance of the tree.

							If False, Ra is the average Ra in the cluster and the diameter 
						is calculated so that the absolute axial resistance Ri is 
						equivalent to the parallel circuit of all the unconnected 
						subtrees. This preserves input resistance

	ALGORITHM
	
	- find next root section of a within-cluster connected subtree
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
		L_eq, diam_eq, Ra_eq, ri_eq, cmtot_eq, gtot_eq = merge_sec_sequential(rootref, allsecrefs)

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


def merge_seg_cluster(cluster, allsecrefs, average_Ri):
	"""
	Calculate cluster equivalent section properties using segment-based merging.

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
					not ref.visited[i] and not cluster_has_parentseg(seg, cluster, clu_secs)))

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
		eq_props, eq_br = merge_seg_subtree(rootseg, allsecrefs)

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


################################################################################
# Utility functions
################################################################################


def cluster_has_parentseg(aseg, cluster, clu_secs):
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


def label_to_order(label):
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


################################################################################
# Equivalent Section creation
################################################################################


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



def make_equivalent_secs(clusters, allsecrefs, gradients):
	"""
	Create the reduced/equivalent cell by creating a section for each cluster 

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


################################################################################
# Reduction Experiments
################################################################################


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
	"""
	Reduce Gillies & Willshaw STN neuron model

	ALGORITHM

	- cluster left and right dendritic tree using either Strahler number criterion
	  or diameter criterion
	
	- use section-based merging within each cluster
	
	- create equivalent sections
		- conductances are interpolated according to distance L/lambda from soma

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
	eq_secs, eq_secrefs = make_equivalent_secs(clusters, allsecrefs, gradients=True)

	# Sort equivalent sections
	eq_sorted = sorted(eq_secrefs, key=lambda ref: label_to_order(ref.cluster_label))
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