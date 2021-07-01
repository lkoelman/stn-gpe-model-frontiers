"""
Reduce Gillies & Willshaw (2006) STN neuron model using the method
described in Bush & Sejnowski (1993)


@author Lucas Koelman
@date	02-02-2017
"""

# Python modules
import math
import pickle
import re

# Enable logging
import logging
logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__) # create logger for this module

# Load NEURON
import neuron
h = neuron.h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# Our own modules
import redutils as redtools
from redutils import ExtSecRef, getsecref # for convenience
from cluster import Cluster
import interpolation as interp
import analyze_reduction as analysis
from gillies_model import gillies_gdict, gillies_mechs, gillies_glist

mechs_chans = gillies_gdict
glist = gillies_glist
gleak_name = 'gpas_STh'

def cluster_sec_properties(clusters, orsecrefs):
	"""
	Calculate cluster properties based on Sections

	@post		properties are stored on each Cluster
	"""
	for cluster in clusters: # SECTION-BASED CLUSTERING
		# Gather member sections
		clu_secs = [secref for secref in orsecrefs if secref.cluster_label==cluster.label]
		logger.debug("Cluster %s contains %i sections" % (cluster.label, len(clu_secs)))

		# Compute equivalent properties
		cluster.eqL = sum(secref.sec.L for secref in clu_secs) / len(clu_secs)
		cluster.eqdiam = math.sqrt(sum(secref.sec.diam**2 for secref in clu_secs))
		cluster.eqRa = sum(secref.sec.Ra for secref in clu_secs) / len(clu_secs)
		cluster.or_area = sum(sum(seg.area() for seg in secref.sec) for secref in clu_secs)
		cluster.or_cmtot = sum(sum(seg.cm*seg.area() for seg in secref.sec) for secref in clu_secs)
		cluster.or_cm = cluster.or_cmtot / cluster.or_area
		cluster.or_gtot = dict((gname, 0.0) for gname in glist) # sum of gbar multiplied by area
		for gname in glist:
			cluster.or_gtot[gname] += sum(sum(getattr(seg, gname)*seg.area() for seg in secref.sec) for secref in clu_secs)

def cluster_seg_properties(clusters, orsecrefs):
	"""
	Calculate cluster properties based on individual segments

	@post		properties are stored on each Cluster
	"""
	for cluster in clusters: # SEGMENT-BASED CLUSTERING
		# Gather member segments
		clu_filter = lambda secref, iseg: secref.cluster_labels[iseg] == cluster.label
		clu_segs = [seg for ref in orsecrefs for i, seg in enumerate(ref.sec) if clu_filter(ref, i)]

		# Compute equivalent properties
		cluster.eqL = sum(seg.sec.L/seg.sec.nseg for seg in clu_segs) / len(clu_segs)
		cluster.eqdiam = math.sqrt(sum(seg.diam**2 for seg in clu_segs))
		cluster.eqRa = sum(seg.sec.Ra for seg in clu_segs) / len(clu_segs)
		cluster.or_area = sum(seg.area() for seg in clu_segs)
		cluster.or_cmtot = sum(seg.cm*seg.area() for seg in clu_segs)
		cluster.or_cm = cluster.or_cmtot / cluster.or_area
		cluster.or_gtot = dict((gname, 0.0) for gname in glist) # sum of gbar multiplied by area
		for gname in glist:
			cluster.or_gtot[gname] += sum(getattr(seg, gname)*seg.area() for seg in clu_segs)

def equivalent_sections(clusters, orsecrefs, f_lambda, 
						gbar_scaling='area', interp_method='linear_neighbors',
						interp_path=None, conserve_gbar_ratios=True):
	"""
	Create equivalent Section for each cluster.

	@type	clusters	list(Cluster)
	@param	clusters	list of Cluster objects (containing cluster label 
						and properties)

	@type	orsecrefs	list(SectionRef)
	@param	orsecrefs	references to all sections in the cell

	@type	gbar_scaling	string
	@param	gbar_scaling	
							- 'area': use the ratio of original area to new area
							to scale Cm and all conductances (including gpas) in
							each section so their total value from the full model
							is conserved. This method is used in the Bush &
							Sejnowski algorithm.

							- 'integral': use ratio or original to new total conductance
							(integral of gbar*area) to scale gbar in each section.
							This conserves the total conductance but not the ration
							of conductances in each segment

							- None: don't use scaling

	@type	interp_method	string
	@param	interp_method	
							'linear_neighbors' 
								= linear interpolation of 'adjacent' segments in full model 
								(i.e. next lower and higher electrotonic path length). 
							
							'linear_dist' 
								= estimate linear distribution and interpolate it
							
							'left_neighbor', 'right_neighbor', 'nearest_neighbor'
								= extrapolation of 'neighoring' segments in terms of L/lambda

	@post				Equivalent section properties are available as
						attributed on given Cluster object

	@return	eq_refs		list(SectionRef) with references to equivalent section
						for each cluster (in same order as param clusters)
	"""
	# List with equivalent section for each cluster
	eq_secs = [None] * len(clusters)

	# Create equivalent sections and passive electric structure
	logger.debug("Building passive section topology...")
	for i, cluster in enumerate(clusters):
		# Create equivalent section
		if cluster.label in [sec.name() for sec in h.allsec()]:
			raise Exception('Section named {} already exists'.format(cluster.label))
		h("create %s" % cluster.label)
		sec = getattr(h, cluster.label)

		# Set geometry 
		sec.L = cluster.eqL
		sec.diam = cluster.eqdiam
		cluster.eq_area = sum(seg.area() for seg in sec) # should be same as cluster eqSurf

		# Passive electrical properties (except Rm/gleak)
		sec.Ra = cluster.eqRa

		# Append to list of equivalent sections
		eq_secs[i] = sec

		logger.debug("Summary for cluster '%s' : L=%f \tdiam=%f \tRa=%f" % (cluster.label,
							cluster.eqL, cluster.eqdiam, cluster.eqRa))

	# Connect equivalent sections
	for i, clu_i in enumerate(clusters):
		for j, clu_j in enumerate(clusters):
			if clu_j is not clu_i and clu_j.parent_label == clu_i.label:
				eq_secs[j].connect(eq_secs[i], clu_j.parent_pos, 0)

	# Set active properties and finetune
	for i, cluster in enumerate(clusters):
		logger.debug("Scaling properties of cluster %s ..." % clusters[i].label)
		sec = eq_secs[i]

		# Insert all mechanisms
		for mech in mechs_chans.keys():
			sec.insert(mech)

		# Scale passive electrical properties
		area_ratio = cluster.or_area / cluster.eq_area
		logger.debug("Ratio of areas is %f" % area_ratio)

		# Scale Cm
		eq_cm1 = cluster.or_cm * area_ratio
		eq_cm2 = cluster.or_cmtot / cluster.eq_area # more accurate than cm * or_area/eq_area
		eq_cm = eq_cm2
		logger.debug("Cm scaled by ratio is %f (equivalently, cmtot/eq_area=%f)" % (eq_cm1, eq_cm2))

		# Scale Rm
		or_gleak = cluster.or_gtot[gleak_name] / cluster.or_area
		eq_gleak = or_gleak * area_ratio # same as reducing Rm by area_new/area_old
		logger.debug("gleak scaled by ratio is %f (old gleak is %f)" % (eq_gleak, or_gleak))

		# Set number of segments based on rule of thumb electrotonic length
		sec.nseg = redtools.calc_min_nseg_hines(100., sec.L, sec.diam, sec.Ra, eq_cm)

		# Save Cm and conductances for each section for reconstruction
		cluster.nseg = sec.nseg # save for reconstruction
		cluster.eq_gbar = dict((gname, [float('NaN')]*cluster.nseg) for gname in glist)
		cluster.eq_cm = [float('NaN')]*cluster.nseg

		# Set Cm and gleak (Rm) for each segment
		if gbar_scaling is not None:
			for j, seg in enumerate(sec):
				setattr(seg, 'cm', eq_cm)
				setattr(seg, gleak_name, eq_gleak)
				cluster.eq_cm[j] = eq_cm
				cluster.eq_gbar[gleak_name][j] = eq_gleak

		# Get active conductances
		active_glist = list(glist)
		active_glist.remove(gleak_name) # get list of active conductances

		# Set initial conductances by interpolation
		if interp_path is None:
			interp_path = (1, (1,3,8)) # dendritic path (default: right tree (1,2,4,6,8))

		if interp_method == 'linear_dist': # estimate distribution and use that
			# First calculate parameters of distribution
			gdist_params = {}
			for gname in active_glist:
				gdist_params[gname] = interp.calc_gdist_params(gname, 
							orsecrefs[0], orsecrefs, interp_path[0], interp_path[1])

			# Then set segment conductances by interpolating/sampling distribution
			for j, seg in enumerate(sec):
				L_elec = redtools.seg_path_L_elec(seg, f_lambda, gleak_name)
				for gname in active_glist:
					if cluster.label.startswith('soma'):
						ibounds, pbounds, gbounds = interp.calc_gdist_params(gname, 
							orsecrefs[0], orsecrefs, -1, (0,))
						logger.debug("Calculated parameters of {} conductace distribution for soma".format(gname))
					else:
						ibounds, pbounds, gbounds = gdist_params[gname]
					gval = interp.interp_gbar_linear_dist(L_elec, ibounds, pbounds, gbounds)
					seg.__setattr__(gname, gval) # Set segment conductance
					cluster.eq_gbar[gname][j] = gval

		else: # linear/nearest/left/right neighbor
			for j, seg in enumerate(sec):
				# First calculate electrotonic path length to segment
				L_elec = redtools.seg_path_L_elec(seg, f_lambda, gleak_name)
				if cluster.label.startswith('soma'):
					tree_id, path_ids = -1, (0,)
				else:
					tree_id, path_ids = interp_path
				path_secs = [secref for secref in orsecrefs if (secref.tree_index == interp_tree_id and 
													secref.table_index in interp_table_ids)]

				# Then get 'neighbor segments'
				bound_segs, bound_L = find_adj_path_segs('path_L_elec', L_elec, path_secs)

				# Set conductances by interpolating neighbors
				for gname in active_glist:
					if interp_method == 'linear_neighbors':
						gval = interp.interp_gbar_linear_neighbors(L_elec, gname, bound_segs, bound_L)
					else:
						match_method = re.search(r'^[a-z]+', interp_method)
						method = match_method.group() # should be nearest, left, or right
						assert len(bound_segs)==1, "Found more than two boundary segments along path"
						gval = interp.interp_gbar_pick_neighbor(L_elec, gname, 
											bound_segs[0], bound_L[0], method)
					seg.__setattr__(gname, gval)
					cluster.eq_gbar[gname][j] = gval


		# Re-scale gbar distribution to yield same total gbar (sum(gbar*area))
		if gbar_scaling is not None:
			for gname in active_glist:
				eq_gtot = sum(getattr(seg, gname)*seg.area() for seg in sec)
				if eq_gtot <= 0.:
					eq_gtot = 1.
				or_gtot = cluster.or_gtot[gname]
				for j, seg in enumerate(sec):
					if gbar_scaling == 'area':
						# conserves ratio in each segment but not total original conductance
						gval = getattr(seg, gname)*clusters[i].or_area/clusters[i].eq_area
					elif gbar_scaling == 'integral':
						# does not conserve ratio but conserves gtot_or since: sum(g_i*area_i * or_area/eq_area) = or_area/eq_area * sum(gi*area_i) ~= or_area/eq_area * g_avg*eq_area = or_area*g_avg
						gval = getattr(seg, gname) * or_gtot/eq_gtot
					else:
						raise Exception("Unknown gbar scaling method'{}'.".format(gbar_scaling))
					seg.__setattr__(gname, gval)
					cluster.eq_gbar[gname][j] = gval # save for reconstruction

		# Check gbar calculation
		for gname in active_glist:
			gtot_or = cluster.or_gtot[gname]
			gtot_eq_scaled = sum(getattr(seg, gname)*seg.area() for seg in sec)
			logger.debug("conductance %s : gtot_or = %.3f ; gtot_eq = %.3f",
							gname, gtot_or, gtot_eq_scaled)

		# Debugging info:
		logger.debug("Created equivalent section '%s' with \n\tL\tdiam\tcm\tRa\tnseg\
		\n\t%.3f\t%.3f\t%.3f\t%.3f\t%d\n\n", clusters[i].label, sec.L, sec.diam, sec.cm, sec.Ra, sec.nseg)
	
	return eq_secs

def rebuild_sections(clusters, eq_secs=None):
	""" Build the reduced model from a previous reduction where the equivalent
		section properties are stored in each Cluster object.

	@param	eq_secs		Section list (if already created)
	"""
	if eq_secs is None:
		existing_sections = False
		eq_secs = [None] * len(clusters)
	else:
		# order according to clusters
		clu_sec = lambda clu: next(sec for sec in eq_secs if sec.name().startswith(clu.label))
		eq_secs = [clu_sec(cluster) for cluster in clusters]
		existing_sections = True

	# Create equivalent section for each cluster
	for i, cluster in enumerate(clusters):
		sec = eq_secs[i]
		if sec is None:
			# Create equivalent section
			if cluster.label in [sec.name() for sec in h.allsec()]:
				raise Exception('Section named {} already exists'.format(cluster.label))
			h("create %s" % cluster.label)
			sec = getattr(h, cluster.label)
			eq_secs[i] = sec

		# Set geometry & global passive properties
		sec.L = cluster.eqL
		sec.diam = cluster.eqdiam
		sec.Ra = cluster.eqRa
		sec.nseg = cluster.nseg

		# Insert all mechanisms
		for mech in mechs_chans.keys():
			sec.insert(mech)
		# Get active conductances
		active_glist = list(glist)
		active_glist.remove(gleak_name) # get list of active conductances

		# Set ion styles
		if not existing_sections:
			sec.push()
			h.ion_style("na_ion",1,2,1,0,1)
			h.ion_style("k_ion",1,2,1,0,1)
			h.ion_style("ca_ion",3,2,1,1,1)
			h.pop_section()

		# Set Cm and conductances for each section (RANGE variables)
		for j, seg in enumerate(sec):
			# passive properties
			setattr(seg, 'cm', cluster.eq_cm[j])
			setattr(seg, gleak_name, cluster.eq_gbar[gleak_name][j])
			# active conductances
			for gname in active_glist:
				if len(cluster.eq_gbar[gname]) != cluster.nseg:
					raise Exception("Number of gbar values does not match number of segments")
				if any(math.isnan(g) for g in cluster.eq_gbar[gname]):
					raise Exception("Conductance vector contains NaN values")
				seg.__setattr__(gname, cluster.eq_gbar[gname][j])

	# Connect equivalent sections
	if not existing_sections:
		for i, clu_i in enumerate(clusters):
			for j, clu_j in enumerate(clusters):
				if clu_j is not clu_i and clu_j.parent_label == clu_i.label:
					eq_secs[j].connect(eq_secs[i], clu_j.parent_pos, 0)
	
	return eq_secs

def save_clusters(clusters, filepath):
	""" Save list of Cluster objects to file """
	try:
		clu_file = open(filepath, "wb")
		pickle.dump(clusters, clu_file)
	finally:
		clu_file.close()
		print("Successfully saved file as " + filepath)

def load_clusters(filepath):
	""" Load list of Cluster objects from file """
	try:
		clu_file = open(filepath, "rb")
		clusters = pickle.load(clu_file)
		return clusters
	finally:
		clu_file.close()

def reduce_bush_sejnowski(delete_old_cells=True):
	""" Reduce STN cell according to Bush & Sejnowski (2016) method """

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
	somaref.tree_index = -1
	somaref.table_index = 0
	for j, dendlist in enumerate((dendLrefs, dendRrefs)):
		for i, secref in enumerate(dendlist):
			secref.tree_index = j
			secref.table_index = i+1

	############################################################################
	# 1. Cluster into functional regions/according to electrotonic length

	# Compute properties used for clustering
	logger.info("Computing electrotonic path lengths...")
	f_lambda = 100. # frequency for electrotonic length
	redtools.assign_electrotonic_length(somaref, allsecrefs, f_lambda, 
										gleak_name, allseg=True)

	# Cluster soma
	somaclu = Cluster('soma')
	somaref.cluster_label = 'soma' # used if section-based clustering
	somaref.cluster_labels = ['soma'] * somaref.sec.nseg # used if segment-based clsutering
	clusters = [somaclu]

	# Cluster segments in each dendritic tree
	thresholds = (0.4, 1.0) # (0.4, 1.0) determine empirically for f=100Hz
	logger.info("Clustering left tree (dend0)...")
	redtools.clusterize_sec_electrotonic(dendLrefs[0], allsecrefs, thresholds, clusters)

	logger.info("Clustering right tree (dend1)...")
	redtools.clusterize_sec_electrotonic(dendRrefs[0], allsecrefs, thresholds, clusters)

	# Determine cluster topology
	redtools.assign_topology(clusters, ['soma', 'trunk', 'smooth', 'spiny'])

	############################################################################
	# 2. Create equivalent sections (compute ppties from cluster sections)
	use_segments = False

	# Calculate cluster properties
	logger.debug("Calculating cluster properties...")
	if use_segments:
		cluster_seg_properties(clusters, orsecrefs)
	else:
		cluster_sec_properties(clusters, orsecrefs)

	# Create new sections
	logger.info("Creating equivalent sections...")
	eq_secs = equivalent_sections(clusters, allsecrefs, f_lambda, 
				use_segments=use_segments, gbar_scaling='area', 
				interp_path=(1, (1,3,8)), interp_method='left_neighbor')
	
	# Delete original model sections & set ion styles
	for sec in h.allsec(): # makes each section the CAS
		if sec.name().startswith('SThcell') and delete_old_cells: # original model sections
			h.delete_section()
		else: # equivalent model sections
			h.ion_style("na_ion",1,2,1,0,1)
			h.ion_style("k_ion",1,2,1,0,1)
			h.ion_style("ca_ion",3,2,1,1,1)

	return clusters, eq_secs

if __name__ == '__main__':
	clusters, eq_secs = reduce_bush_sejnowski(delete_old_cells=True)
	filepath = "C:\\Users\\lkoelman\\cloudstore_m\\simdata\\bush_sejnowski\\stn_reduced_bush.p"
	save_clusters(clusters, filepath)
