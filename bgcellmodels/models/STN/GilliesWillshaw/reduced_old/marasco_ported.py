"""
This is an attempt to a literal port of the code & algorithms supplied
with the Marasco & Miglieri (2013) paper to Python. It is incomplete
and untested, and was solely used to gain insight into the algorithms.

DEPRECATED: all functionality is in reduce_marasco and redutils

@author Lucas Koelman
@date	28-11-2016
"""
# Python modules
from collections import OrderedDict
from operator import attrgetter
import math
PI = math.pi
import logging
logger = logging.getLogger(__name__) # create logger for this module
logger.setLevel(logging.DEBUG)

# NEURON modules
import neuron
h = neuron.h

# Own modules
from redutils import * # since general functions were moved to that module

################################################################################
# Merging
################################################################################

def calc_mrgRiSurf(secref):
	""" Calculate axial resistance and surface of section
		based on mutable its 'mrg' properties """
	# textbook definition of axial resistance
	mrgri = (secref.sec.Ra*secref.mrgL)/(PI*secref.mrgdiam**2/4*100)
	# cylinder surface based on merged L and diam
	mrgsurf = secref.mrgL*secref.mrgdiam*PI
	return mrgri, mrgsurf

def prep_merge(secrefs):
	""" Prepare sections for merging procedure by computing
		and storing metrics needed in the merging operation
	"""
	for secref in secrefs: # mark all as unvisited and unmerged
		sec = secref.sec

		# geometrical properties
		secref.secSurf = 0 # total area (diam/area is range var)
		secref.secri = 0 # total axial resistance from 0 to 1 end (seg.ri is between segments)
		for seg in sec:
			secref.secSurf += seg.area()
			secref.secri += seg.ri()

		# Axial path resistance
		sec_path_ri(secref)
		
		# mutable properties for merging
		secref.mrgL = abs(sec.L)
		secref.mrgdiamSer = sec.diam
		secref.mrgdiamPar = sec.diam**2
		secref.mrgdiam = math.sqrt(sec.Ra*sec.L*4./secref.secri/100/PI)
		secref.mrgri = secref.secri
		secref.mrgri2 = secref.secri
		secref.mrgSurf = secref.secSurf
		
		# properties for merging iteration
		secref.visited = False
		secref.doMerge = False

def find_mergeable(secrefs, clusterlabel):
	""" Find next mergeable sections

	@param secrefs			list of mutable SectionRef containing all sections
	@param clusterlabel		label of the cluster to merge
	@effect					- find next available section A and its sibling B
							- set them to visited (unavailable for search)
							- find their parent P if available and within cluster
							- if P not found: set A and B available for merge
	@return 				tuple of Section Ref (secA, secB, secP)
	"""

	# Get sections in current cluster and sort by order
	clustersecs = [sec for sec in secrefs if sec.cluster_label == clusterlabel]
	
	# Sort by order (nseg from soma), descending
	secsbyorder = clustersecs.sort(key=attrgetter('order'), reverse=True)

	# Keep looping until all sections visited
	while any(not sec.visited for sec in clustersecs):
		# get next unvisited section furthest from soma
		secA = next(sec for sec in secsbyorder if not sec.visited)
		secA.visited = True

		# get its brother section (same parent) if available
		secB = next((sec for sec in secsbyorder if (sec is not secA 
					and not sec.visited and sameparent(sec, secA))), None)
		if secB is not None: secB.visited = True

		# get their parent section if in same cluster
		secP = secA.parent
		secP = getsecref(secP, clustersecs) # None if not in cluster
		if secP is None:
			if secA: secA.doMerge = True
			if secB: secB.doMerge = True

		# Return sec refs from generator
		yield secA, secB, secP

def mergeChildWithParent(refA, refP):
	"""
	Merge child into parent section

	SOURCE: mergingParentSec() -> mergingSerialMethod()
	"""
	
	# Combine properties (call to shortnms)
	Amri, Amsurf = calc_mrgRiSurf(refA)
	Pmri, Pmsurf = calc_mrgRiSurf(refP)

	# Combine (call to mergingSerialMethod())
	newL = refA.mrgL + refP.mrgL
	newRa = (refA.sec.Ra+refP.sec.Ra)/2
	newri = Amri + Pmri
	newSurf = max(Amsurf, Pmsurf)
	newdiam = math.sqrt(newRa*newL*4/newri/PI/100.)

	# Update parent properties
	refP.mrgL = newL
	refP.mrgri = newri
	refP.mrgri2 = newri
	refP.mrgSurf = newSurf
	refP.mrgdiam = newdiam
	# FIXME: unassigned/unused
	# refP.mrgdiamSer = refA.mrgdiamSer + refP.mrgdiamSer
	# refP.mrgdiamPar = refA.mrgdiamPar + refP.mrgdiamPar

def mergeYWithParent():
	"""
	Merge Y (two branched children) into parent

	SOURCE: mergingYSec() -> mergingYMethod() -> (mergingParallelMethod() + mergingSerialMethod())
	"""
	# Current equivalent properties (call to shortnms)
	Amri, Amsurf = calc_mrgRiSurf(refA)
	Bmri, Bmsurf = calc_mrgRiSurf(refB)
	Pmri, Pmsurf = calc_mrgRiSurf(refP)

	# Combine properties of parallel sections (Call to mergingParallelMethod())
	L12 = (refA.mrgL*Amsurf+refB.mrgL*Bmsurf)/(Amsurf+Bmsurf) # lengths weighted by surfaces
	Ra12 = (refA.sec.Ra+refB.sec.Ra)/2
	diam12 = math.sqrt(refA.mrgdiam**2 + refB.mrgdiam**2)
	cross12 = PI*diam12**2/4 # cross-section area

	# Equivalent propties of merged parallel sections
	ri12 = Ra12*L12/cross12/100 # equivalent axial resistance
	surf12 = L12*diam12*PI # equivalent surface

	# Combine properties of serial sections (Call to mergingSerialMethod())
	newL = L12 + refP.mrgL
	newRa = (Ra12+refP.sec.Ra)/2
	newri = ri12 + Pmri # NOTE: not uses as r_a in paper (see newri2 which corresponds to r_a)
	newdiam = math.sqrt(newRa*newL*4/newri/PI/100.)
	# FIXME: unassigned/unused
	# newSurf = max(surf12, Pmsurf) # overwritten by newSurf in mergingYMethod()
	# newdiamSer = refA.mrgdiamSer + refP.mrgdiamSer
	# newdiamPar = refA.mrgdiamPar + refP.mrgdiamPar

	# Equivalent properties of merged serial sections
	newri2 = (Amri*Bmri/(Amri+Bmri))+Pmri # NOTE: this one is used as ri (r_a in paper)
	newSurf = max(Amsurf+Bmsurf,Pmsurf)

	# Update parent properties
	refP.mrgL = newL
	refP.mrgri = newri
	refP.mrgri2 = newri2
	refP.mrgSurf = newSurf
	refP.mrgdiam = newdiam
	# FIXME: unassigned/unused
	# refP.mrgdiamSer = newdiamSer + refP.mrgdiamSer
	# refP.mrgdiamPar = newdiamPar + refP.mrgdiamPar

def calc_eq_ppties_mrg(cluster, allsecrefs):
	"""
	Calculate properties of equivalent section for cluster that are based
	on the sections marked for merging

	SOURCE: see mergingAllMethod() in mergingMethods.hoc()
	"""
	# Initialize cluster properties
	cluster.maxL = 0 # L of longest sec in cluster
	cluster.totsurf = 0 # sum of surf calc from mrg dimensions for each merge available sec
	cluster.totdiam = 0 # sum of mrgdiam for each merge available sec
	cluster.eqri2Prod = 1
	cluster.eqri2Sum = 0
	cluster.surfSum = 0 # sum of mrgSum property for each merge available sec

	# Gather sections for merging
	mrg_secs = [secref for secref in allsecrefs if (secref.doMerge and 
									secref.cluster_label==cluster.label)]
	cluster.weightsL = [1.0]*len(mrg_secs)

	# Update cluster properties based on mergeable sections
	for i, secref in enumerate(mrg_secs):

		if cluster.maxL < secref.mrgL:
			cluster.maxL = secref.mrgL
		cluster.totsurf += secref.mrgL*PI*secref.mrgdiam
		cluster.totdiam += secref.mrgdiam
		cluster.eqri2Prod *= secref.mrgri2
		cluster.eqri2Sum += secref.mrgri2
		cluster.surfSum += secref.mrgsurf

		# Equivalent dimensions
		cluster.normalFact = totsurf
		cluster.weightsL[i] = secref.mrgL*PI*secref.mrgdiam
		cluster.eqL += secref.mrgL*cluster.weightsL[i]
		cluster.eqdiam += (secref.sec.Ra*cluster.maxL*4.)/(secref.mrgri*PI*100.) # rho squared

	cluster.eqdiam = math.sqrt(cluster.eqdiam)
	cluster.eqri2 = cluster.eqri2Sum/len(mrg_secs)
	cluster.eqSurf = cluster.eqL*PI*cluster.eqdiam
	cluster.orSurf1 = cluster.surfSum

	return cluster

def calc_eq_ppties_all(cluster, allsecrefs):
	"""
	Calculate properties of equivalent section for cluster that are based
	on all of its sections

	@param allsecrefs	list of SectionRef (mutable) with first element
						a ref to root/soma section

	SOURCE: see first part of createPurkEqCell() in useful&InitProc.hoc
	"""

	somaref = allsecrefs[0]

	# Compute cluster properties based on all its sections
	clu_secs = [secref for secref in secrefs if secref.cluster_label==cluster.label]
	cluster.eqRa = sum(secref.sec.Ra for secref in clu_secs)/len(clu_secs)
	cluster.orMaxpathri = max(secref.pathri1 for secref in clu_secs)
	cluster.orMinpathri = min(secref.pathri0 for secref in clu_secs)

	# Clusters that contain only one section (trunk/soma)
	if len(clu_secs) == 1:
		secref = clu_secs[0]
		if secref.has_parent() and secref.parent() is not somaref.sec:
			parref = getsecref(secref.parent(), allsecrefs)
			orMinpathri = parref.pathri1
		else:
			orMinpathri = 0

	# Correct if min > max
	if orMinpathri > orMaxpathri:
		orMinpathri, orMaxpathri = orMaxpathri, orMinpathri # swap


################################################################################
# Equivalent/reduced cell
################################################################################

def create_equivalent_sections(eq_clusters, allsecrefs):
	""" Create the reduced/equivalent cell by creating 
		a section for each cluster 

	@param eq_clusters	list of Cluster objects containing data
						for each cluster
	@param allsecrefs	list of SectionRef (mutable) with first element
						a ref to root/soma section

	SOURCE: in RedPurk.hoc/init(): see call to createPurkEqCell() and everything that follows it
	"""
	# These are scaling factors for each cluster
	factCm = 1.0 # fixed/not set in Marasco 2012/2013 method
	factRm = 1.0 # fixed/not set in Marasco 2012/2013 method
	factRa = 1.0 # fixed/not set in Marasco 2012/2013 method
	factEqRa = 1.0 # fixed/not set in Marasco 2012/2013 method

	# Scaling factors based on surface
	SURFFACT_PARTIAL_TOTAL = 1
	SURFFACT_PARTIAL_METHOD = 2
	# Vector containing sum of original section surface for each cluster
	orSurfVect = [sum(sec.secSurf for sec in allsecrefs if sec.clusterlabel == clu.label) for clu in eq_clusters]
	eqSurfVect = [clu.eqdiam*PI*clu.eqL for clu in eq_clusters]
	# Total original surface and total surface of equivalent sections
	TotOrSurf = sum(sec.secSurf for sec in allsecrefs)
	TotEqSurf = sum(clu.eqdiam*PI*clu.eqL for clu in eq_clusters)
	surfFact = [1.0]*len(eq_clusters)
	surfFactRmCm = [1.0]*len(eq_clusters)
	if SURFFACT_PARTIAL_METHOD == 0:
		surfFact = [orSurf/eqSurf for orSurf,eqSurf in zip(orSurfVect,eqSurfVect)]
		surfFactRmCm = [orSurf/eqSurf for orSurf,eqSurf in zip(orSurfVect,eqSurfVect)]
	if SURFFACT_PARTIAL_TOTAL == 0:
		surfFact = [TotOrSurf/TotEqSurf]*len(eq_clusters)
		surfFactRmCm = [TotOrSurf/TotEqSurf]*len(eq_clusters)



	# Create equivalent section for each clusters
	eq_sections = [h.Section() for clu in eq_clusters]
	for i, clu_i in enumerate(eq_clusters):
		for j, clu_j in enumerate(eq_clusters):
			if clu_j is not clu_i and clu_j.parent_label == clu_i.label:
				eq_sections[j].connect(eq_sections[i], clu_j.parent_pos, 0)

	# Set L/diam/Ra/cm/g_leak for each equivalent section
	
	# Calculate equivalent path resistance for each cluster/equivalent section
	# FIXME: this was wrong, calculation happend on equivalent sections, after created
	for i, secref in enumerate(eq_secs):

		# Cluster corresponding to current section
		cluster = next(clu for clu in clusterList if clu.label == secref.cluster_label)
		cluster.eqsecpathri1 = 0.0
		cluster.eqsecpathri0 = 0.0

		# Get path from root node to this sections
		rootsec = treeroot(secref)
		calc_path = h.RangeVarPlot('v')
		rootsec.push()
		calc_path.begin(.5)
		sec.push()
		calc_path.end(.5)
		root_path = h.SectionList()
		calc_path.list(root_path) # store path in list
		h.pop_section()
		h.pop_section()

		# Compute axial path resistances
		path_secs = list(root_path)
		path_len = len(root_path)
		for i, psec in enumerate(path_secs):
			for seg in psec:
				cluster.eqsecpathri1 += seg.ri()
				if i < path_len-1:
					cluster.eqsecpathri0 += seg.ri()

################################################################################
# Main routine
################################################################################

def reduce_marasco():
	""" Implementation of Marasco (2013) CellPurkAnalysis() & PurkReduction() """

	# Initialize Gillies model
	h.xopen("gillies_cell_singleton.hoc")
	# Make sections accesible by both name and index + allow to add attributes
	somaref = ExtSecRef(sec=h.SThcell[0].soma)
	dendLrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend0]
	dendRrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend1]
	alldendrefs = dendLrefs + dendRrefs
	allsecrefs = [somaref] + alldendrefs

	# Cluster Soma
	somaref.cluster_label = 'soma'
	somaref.parent_label = 'soma'
	somaref.parent_pos = 0.0
	
	# Cluster subtree of each trunk section
	clu_fun = clutools.label_from_strahler
	clu_args = {'thresholds':(3,5)}
	clutools.clusterize_custom(dendLrefs[0], allsecrefs, clusters, '_0', clu_fun, clu_args)
	clutools.clusterize_custom(dendRrefs[0], allsecrefs, clusters, '_1', clu_fun, clu_args)
	cluster_labels = list(set(secref.cluster_label for secref in allsecrefs)) # unique labels
	eq_clusters = [Cluster(label) for label in cluster_labels]

	# Merge within-cluster branches until islands remain
	for cluster in eq_clusters:

		# First merge children into parents (iteratively update parent properties)
		for secA, secB, secP in find_mergeable(alldendrefs, cluster.label):
			# Now merge these sections
			if secP is not None and secB is not None:
				mergeYWithParent()
			elif secP is not None:
				mergeChildWithParent()
				

		# Calculate cluster/equivalent section properties
		calc_eq_ppties_mrg(cluster, allsecrefs)
		calc_eq_ppties_all(cluster, allsecrefs)

	# Initialize the equivalent model (one equivalent sec per cluster)
	create_equivalent_sections(eq_clusters, allsecrefs)

if __name__ == '__main__':
	stn_tree = redutils.combinedtree()
	assign_strahler(stn_tree)