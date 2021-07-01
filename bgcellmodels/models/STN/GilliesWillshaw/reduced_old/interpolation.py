"""
Manual reduction of Gillies & Willshaw (2006) STN neuron model


@author Lucas Koelman
@date	09-02-2017
"""

# Math modules
import numpy as np

# Enable logging
import logging
logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logger = logging.getLogger('redops') # create logger for this module

# Load NEURON
import neuron
h = neuron.h

# Our own modules
from bgcellmodels.common.nrnutil import ExtSecRef, getsecref # for convenience
from bgcellmodels.models.STN.GilliesWillshaw.gillies_model import (
	gillies_gdict, gillies_glist)

mechs_chans = gillies_gdict
glist = gillies_glist
gleak_name = 'gpas_STh'


def loadgstruct(gext):
	""" Return structured array with conductance values for given
		channel.

	@param gext	channel name/file extension in sth-data folder
	"""
	gfile = os.path.normpath(os.path.join(scriptdir, 'sth-data', 'cell-'+gext))
	gstruct = np.loadtxt(gfile, dtype={'names': ('dendidx', 'branchidx', 'x', 'g'),
									   'formats': ('i4', 'i4', 'f4', 'f4')})
	return np.unique(gstruct) # unique and sorted rows


def gillies_gstructs():
	""" Load all structured arrays with conductance values """
	gfromfile = ["gk_KDR", "gk_Kv31", "gk_Ih", "gk_sKCa", "gcaT_CaT", "gcaN_HVA", "gcaL_HVA"]
	gmats = {gname: loadgstruct(gname) for gname in gfromfile}
	return gmats


def calc_gdist_params(gname, secref, orsecrefs, tree_index, path_indices, xgvals=None):
	"""
	Calculate parameters of the linear conductance distribution
	as defined in the Gillies & Willshaw paper (A/B/C/D)

	@type	orsecrefs	list(h.SectionRef)
	@param	orsecrefs	References to sections in original model. Eac SectionRef
						must have properties pathLelec0, pathLelec1, and pathL_elec
						containing the electrotonic path length up to the 0 and
						1-end, and up to each segment

	@return			tuple (L0, g0), (L1, g1), (Lstart, gstart), (Lend, gend), (gmin, gmax)
					gmin: 	min gbar along path
					gmax:	max gbar along path
					L0: 	electrotonic length of first segment with gbar > gmin
					L1:		electrotonic length of last segment with gbar > gmin
					Lstart:	electrotonic length of first segment on path
					Lend:	electrotonic length of last segment on path
	"""
	if secref is None:
		return
	first_call = xgvals is None
	if first_call:
		xgvals = []
	secfilter = lambda secref: ((secref.tree_index==tree_index) 
								and (secref.table_index in path_indices))

	# Measure current section
	if secfilter(secref):
		for i, seg in enumerate(secref.sec): # for each segment: save electrotonic path length & gbar
			xgvals.append((secref.pathL_elec[i], getattr(seg, gname)))
	
	# Collect x and gbar values from chidlren
	for childsec in secref.sec.children():
		childref = getsecref(childsec, orsecrefs)
		calc_gdist_params(gname, childref, orsecrefs, tree_index, path_indices, xgvals=xgvals) # this updates xgvals

	# Calculate parameters from gbar distribution
	if first_call:
		xgvals = sorted(xgvals, key=lambda xg: xg[0]) # sort by ascending x (L_elec)
		xvals, gvals = zip(*xgvals)
		gmin = min(gvals)
		gmax = max(gvals)
		xg0 = next((xg for xg in xgvals if xg[1] > gmin), xgvals[0]) # first g > gmin
		xg1 = next((xg for xg in reversed(xgvals) if xg[1] > gmin), xgvals[-1]) # last g > gmin
		return (xg0, xg1), (xgvals[0], xgvals[-1]), (gmin, gmax)


def find_adj_path_segs(interp_prop, interp_L, path_secs):
	"""
	Find segments that are immediately before or after given path length value.

	@param interp_prop	property used for calculation of path length, one of
						'path_L', 'path_ri', or 'path_L_elec'

	@param interp_L		path length value

	@param path_secs	list of SectionRef or EqProp objects containing stored properties
						for each Section on the path

	@pre				path lengths must be assigned along searched path
	"""
	# attribute names corresponding to path property (see redutils.sec_path_props())
	if interp_prop == 'path_L':
		sec_prop0 = 'pathL0'
		sec_prop1 = 'pathL1'
		seg_prop = 'pathL_seg'
	
	elif interp_prop == 'path_ri':
		sec_prop0 = 'pathri0'
		sec_prop1 = 'pathri1'
		seg_prop = 'pathri_seg'
	
	elif interp_prop == 'path_L_elec':
		sec_prop0 = 'pathLelec0'
		sec_prop1 = 'pathLelec1'
		seg_prop = 'pathL_elec'
	
	else:
		raise ValueError("Unknown path property '{}'".format(interp_prop))

	# 1. Find sections where interp_L falls between path length to start and end
	# filter_path = lambda secref: ((secref.tree_index==tree_index) and (secref.table_index in path_indices))
	# path_secs = [secref for secref in orsecrefs if filter_path(secref)]
	filter_L = lambda secref: (getattr(secref, sec_prop0) <= interp_L <= getattr(secref, sec_prop1))
	map_secs = [secref for secref in path_secs if filter_L(secref)]

	# If none found: extrapolate
	if len(map_secs) == 0:
		# find section where L(1.0) is closest to interp_L
		L_closest = min([abs(interp_L-getattr(secref, sec_prop1)) for secref in path_secs])
		# all sections where interp_L-L(1.0) is within 5% of this
		map_secs = [secref for secref in path_secs if (0.95*L_closest <= abs(interp_L-getattr(secref, sec_prop1)) <= 1.05*L_closest)]
		logger.debug("Path length did not map onto any original section: " 
					 "extrapolating from {} sections".format(len(map_secs)))

	# 2. In mapped sections: find segments that bound interp_L
	bound_segs = [] # bounding segments (lower, upper)
	bound_L = []    # path lengths to (lower, upper)
	for secref in map_secs:
		
		# in each section find the two segments with L(seg_a) <= interp_L <= L(seg_b)
		if isinstance(path_secs[0], ExtSecRef):
			segs_internal = [seg for seg in secref.sec]
			vals_internal = getattr(secref, seg_prop)
		else:
			segs_internal = secref.seg # array of dict() representing each segment
			vals_internal = [sprops[seg_prop] for sprops in secref.seg]

		if len(segs_internal) == 1: # single segment: just use midpoint
			
			midseg = segs_internal[0]
			midL = vals_internal[0]
			bound_segs.append((midseg, midseg))
			bound_L.append((midL, midL))

		else:
			# Get last/furthest segment with L <= interp_L
			lower = [(i, pval) for i, pval in enumerate(vals_internal) if interp_L >= pval]
			i_a, L_a = next(reversed(lower), (0, vals_internal[0])) # default is first seg

			# Get first/closest segment with L >= interp_L
			higher = ((i, pval) for i, pval in enumerate(vals_internal) if interp_L <= pval)
			i_b, L_b = next(higher, (-1, vals_internal[-1])) # default is last seg

			# Append bounds
			bound_segs.append((segs_internal[i_a], segs_internal[i_b]))
			bound_L.append((L_a, L_b))

	# Return pairs of boundary segments and boundary path lengths
	return bound_segs, bound_L


def interp_gbar_linear_dist(L_elec, bounds, path_bounds, g_bounds):
	""" Linear interpolation of gbar according to the given
		electrotonic path length and boundaries.

	@param	L_elec		electrotonic path length where gbar is desired

	@param	bounds		((L0, g0), (L1, g1)): electrotonic path length of
						first and last segment with gbar > gmin and their
						gbar values

	@param	path_bounds	((L_start, g_start), (L_end, g_end)): electrotonic
						path length of first and last segment on path and
						their gbar values

	@param	g_bounds	min and max gbar value
	"""
	L_a, g_a = bounds[0]
	L_b, g_b = bounds[1]
	L_min, g_start = path_bounds[0]
	L_max, g_end = path_bounds[1]
	gmin, gmax = g_bounds

	if L_a <= L_elec <= L_b:
		if L_b == L_a:
			alpha = 0.5
		else:
			alpha = (L_elec - L_a)/(L_b - L_a)
		if alpha > 1.0:
			alpha = 1.0
		if alpha < 0.0:
			alpha = 0.0
		return g_a + alpha * (g_b - g_a)
	elif L_elec < L_a:
		if L_a == L_min:
			return g_a # proximal distribution, before bounds
		else:
			return gmin # distal distribution, before bounds
	else: #L_elec > L_b:
		if L_b == L_max:
			return g_b # distal distribution, after bounds
		else:
			return gmin # proximal distribution, after bounds


def interp_gbar_linear_neighbors(L_elec, gname, bound_segs, bound_L):
	""" For each pair of boundary segments (and corresponding electrotonic
		length), do a linear interpolation of gbar in the segments according
		to the given electrotonic length. Return the average of these
		interpolated values.

	@type	gname		str
	@param	gname		full conductance name (including mechanism suffix)

	@type	bound_segs	list(tuple(Segment,Segment)) OR list(tuple(dict,dict)
	@param	bound_segs	pairs of boundary segments

	@type	bound_L		list(tuple(float, float))
	@param	bound_segs	electrotonic lengths of boundary segments

	@return		gbar_interp: the average interpolated gbar over all boundary
				pairs
	"""
	# Allow to give a dict instead of segment
	if isinstance(bound_segs[0][0], dict):
		getfunc = dict.__getitem__
	else:
		getfunc = getattr

	# Interpolate
	gbar_interp = 0.0
	for i, segs in enumerate(bound_segs):
		seg_a, seg_b = segs
		L_a, L_b = bound_L[i]
		assert L_b >= L_a, "Path lengths of bounding segments must be in ascending order."

		# Linear interpolation of gbar in seg_a and seg_b according to electrotonic length
		if L_elec <= L_a:
			gbar_interp += getfunc(seg_a, gname)
			continue
		if L_elec >= L_b:
			gbar_interp += getfunc(seg_b, gname)
			continue
		if L_b == L_a:
			alpha == 0.5
		else:
			alpha = (L_elec - L_a)/(L_b - L_a)
		if alpha > 1.0:
			alpha = 1.0 # if too close to eachother
		gbar_a = getfunc(seg_a, gname)
		gbar_b = getfunc(seg_b, gname)
		gbar_interp += gbar_a + alpha * (gbar_b - gbar_a)

	gbar_interp /= len(bound_segs) # take average
	return gbar_interp


def interp_gbar_pick_neighbor(L_elec, gname, bound_segs, bound_L, method='nearest'):
	""" Nearest/left/right neighbor interpolation 

	@param	method	interpolation method: 'nearest'/'left'/'right'

	@see	interp_gbar_linear_neighbors() for other parameters
	"""
	seg_a, seg_b = bound_segs
	L_a, L_b = bound_L
	if isinstance(seg_a, dict):
		getfunc = type(seg_a).__getitem__
	else:
		getfunc = getattr

	if method == 'left':
		return getfunc(seg_a, gname)
	elif method == 'right':
		return getfunc(seg_b, gname)
	elif method == 'nearest':
		if abs(L_elec-L_a) <= abs(L_elec-L_b):
			return getfunc(seg_a, gname)
		else:
			return getfunc(seg_b, gname)
	else:
		raise Exception("Unknown interpolation method '{}'".format(method))


def interpconductances(sec, tree_index, path_indices, glist=None):
	""" Interpolate conductances along given path of branches using
		tabular data provided with Gillies & Willshaw paper

	@type	sec			Hoc.Section()
	@param	sec			section to set conductances for

	@type	tree_index		int
	@param	tree_index		index of the dendritic tree

	@type	path_indices		sequence of int
	@param	path_indices		indices of branches in full model, specifying path
						along dendritic tree from soma outward

	@type	glist		dict<str, str>
	@param	glist		list of conductances to set (including mechanism suffix)
	"""

	# Load channel conductances from file
	allgmats = redutils.loadgstructs()
	if glist is None:
		glist = list(gillies_glist)

	# Na & NaL are not from file
	h("default_gNa_soma = 1.483419823e-02") 
	h("default_gNa_dend = 1.0e-7")
	h("default_gNaL_soma = 1.108670852e-05")
	h("default_gNaL_dend = 0.81e-5")

	# branch indices along longest path
	geostruct = redutils.loadgeotopostruct(tree_index)
	pathL = np.array([geostruct[i-1]['L'] for i in path_indices]) # length of each branch along path

	# Distributed conductances: interpolate each conductance along longest path
	for seg in sec:
		lnode = seg.x*sum(pathL) # equivalent length along longest path

		# first determine on which branch we are and how far on it
		nodebranch = np.NaN # invalid branch
		xonbranch = 0
		for i, branchidx in enumerate(path_indices): # map node to branch and location
			if lnode <= pathL[0:i+1].sum(): # location maps to this branch
				begL = pathL[0:i].sum()
				endL = pathL[0:i+1].sum()
				nodebranch = branchidx
				xonbranch = (lnode-begL)/(endL-begL) # how far along this branch are we
				break
		if np.isnan(nodebranch):
			raise Exception('could not map to branch')

		# now interpolate all conductances from file
		for gname, gmat in allgmats.iteritems():
			if gname not in glist:
				print('Skipping conductance: '+gname)
				continue
			branchrows = (gmat['dendidx']==tree_index) & (gmat['branchidx']==nodebranch-1)
			gnode = np.interp(xonbranch, gmat[branchrows]['x'], gmat[branchrows]['g'])
			sec(xnode).__setattr__(gname, gnode)

		# Conductances with constant value (vals: see tools.hoc/washTTX)
		gNa = 1.483419823e-02 if tree_index==-1 else 1.0e-7 # see h.default_gNa_soma/dend in .hoc file
		gNaL = 1.108670852e-05 if tree_index==-1 else 0.81e-5 # see h.default_gNaL_soma/dend in .hoc file
		gNarsg = 0.016 # same as in .mod file and Akeman papeer
		g_fixed = {'gna_Na':gNa, 'gna_NaL':gNaL, 'gbar_Narsg':gNarsg} # NOTE: Narsg is NOT in Gillies model
		for gname, gval in g_fixed.iteritems():
			if gname in glist:
				setattr(seg, gname, gval)


def setconductances(sec, dendidx, fixbranch=None, fixloc=None, glist=None):
	""" Set conductances at the node/midpoint of each segment
		by interpolating values along longest path
		(e.g. along branch 1-2-5 in dend1)

	@param dendidx		index of the denritic tree where longest path should
						be followed (1/0/-1)

	@param fixbranch	if you want to map the section to a fixed branch
						instead of following the longest path, provide its index

	@param fixloc		if you want to map all segments/nodes
						to a fixed location on the mapped branch,
						provide a location (0<=x<=1)

	@type	glist		dict<str, str>
	@param	glist		list of conductances to set (including mechanism suffix)
	"""

	# Load channel conductances from file
	allgmats = redutils.loadgstructs()
	if glist is None:
		glist = list(gillies_glist)

	# Na & NaL are not from file
	h("default_gNa_soma = 1.483419823e-02")
	h("default_gNa_dend = 1.0e-7")
	h("default_gNaL_soma = 1.108670852e-05")
	h("default_gNaL_dend = 0.81e-5")

	# branch indices along longest path
	if dendidx == 1:
		geostruct = redutils.loadgeotopostruct(dendidx)
		longestpath = np.array([1,2,5])
		pathL = np.array([geostruct[i-1]['L'] for i in longestpath])
	elif dendidx == 0:
		geostruct = redutils.loadgeotopostruct(dendidx)
		longestpath = np.array([1,2,4,7])
		pathL = np.array([geostruct[i-1]['L'] for i in longestpath])
	else: # -1: soma
		# dimensions not in treeX-nom.dat file
		longestpath = np.array([1]) # soma is dendidx=-1, branchidx=0 in file
		pathL = np.array([18.8])

	# Distributed conductances: interpolate each conductance along longest path
	for iseg in range(1, sec.nseg+1):
		xnode = (2.*iseg-1.)/(2.*sec.nseg) # arclength of current node (segment midpoint)
		lnode = xnode*sum(pathL) # equivalent length along longest path

		# first determine on which branch we are and how far on it
		if (fixbranch is not None) and (fixloc is not None):
			nodebranch = fixbranch
			xonbranch = fixloc
		else:
			nodebranch = np.NaN # invalid branch
			xonbranch = 0
			for i, branchidx in enumerate(longestpath): # map node to branch and location
				if lnode <= pathL[0:i+1].sum(): # location maps to this branch
					begL = pathL[0:i].sum()
					endL = pathL[0:i+1].sum()
					nodebranch = branchidx
					xonbranch = (lnode-begL)/(endL-begL) # how far along this branch are we
					break
			if np.isnan(nodebranch):
				raise Exception('could not map to branch')

		# now interpolate all conductances from file
		for gname, gmat in allgmats.iteritems():
			if gname not in glist:
				print('Skipping conductance: '+gname)
				continue
			branchrows = (gmat['dendidx']==dendidx) & (gmat['branchidx']==nodebranch-1)
			gnode = np.interp(xonbranch, gmat[branchrows]['x'], gmat[branchrows]['g'])
			sec(xnode).__setattr__(gname, gnode)

		# Conductances with constant value (vals: see tools.hoc/washTTX)
		gNa = 1.483419823e-02 if dendidx==-1 else 1.0e-7 # see h.default_gNa_soma/dend in .hoc file
		gNaL = 1.108670852e-05 if dendidx==-1 else 0.81e-5 # see h.default_gNaL_soma/dend in .hoc file
		gNarsg = 0.016 # same as in .mod file and Akeman papeer
		g_fixed = {'gna_Na':gNa, 'gna_NaL':gNaL, 'gbar_Narsg':gNarsg} # NOTE: Narsg is NOT in Gillies model
		for gname, gval in g_fixed.iteritems():
			if gname not in glist: continue
			sec(xnode).__setattr__(gname, gval)