"""
Manual reduction of Gillies & Willshaw (2006) STN neuron model


@author Lucas Koelman
@date   09-02-2017
"""

# Load NEURON
import neuron
h = neuron.h

# Our own modules
from bgcellmodels.common import logutils
from bgcellmodels.common.nrnutil import ExtSecRef, get_range_var
from . import redutils

# Enable logging
logger = logutils.getLogger('redops') # create logger for this module


def find_adj_path_segs(interp_prop, interp_L, path_secs):
    """
    Find segments that are immediately before or after given path length value.

    @param interp_prop  property used for calculation of path length, one of
                        'path_L', 'path_ri', or 'path_L_elec'

    @param interp_L     path length value

    @param path_secs    list of SectionRef or EqProp objects containing stored properties
                        for each Section on the path

    @pre                path lengths must be assigned along searched path
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


def get_interpolation_path_sections(secref, tree_properties):
    """
    Return Sections forming a path from soma to dendritic terminal/endpoint.
    The path is used for interpolating spatially non-uniform properties.

    @param  secref      <SectionRef> section for which a path is needed

    @return             <list(Section)>
    """

    # Return each Section/SecProps in secref.merged_sec_gids
    absorbed = redutils.find_secprops(tree_properties,
                            lambda sec: sec.gid in secref.merged_sec_gids)

    # Find the farthest one (highest path length from soma)
    max_dist = max((sec.pathL1 for sec in absorbed))
    farthest_sec = next((sec for sec in absorbed if sec.pathL1==max_dist))

    # Descend farther, following child branch with largest mean diameter
    def mean_diam(sec):
        return sum((sec.seg[i]['diam'] for i in range(sec.nseg))) / sec.nseg
    
    terminal_sec = farthest_sec
    while any(terminal_sec.children):
        best_diam = 0.0
        for child in terminal_sec.children:
            child_diam = mean_diam(child)
            if child_diam > best_diam:
                best_diam = child_diam
                terminal_sec = child
        assert best_diam > 0.0

    # Construct path to terminal section
    path_secs = [terminal_sec]
    while path_secs[-1].parent is not None:
        path_secs.append(path_secs[-1].parent)

    path_secs = list(reversed(path_secs)) # soma to terminal order
    
    if path_secs is None: raise Exception('No path sections found')
    return path_secs


def interp_gbar_linear_dist(L_elec, bounds, path_bounds, g_bounds):
    """ Linear interpolation of gbar according to the given
        electrotonic path length and boundaries.

    @param  L_elec      electrotonic path length where gbar is desired

    @param  bounds      ((L0, g0), (L1, g1)): electrotonic path length of
                        first and last segment with gbar > gmin and their
                        gbar values

    @param  path_bounds ((L_start, g_start), (L_end, g_end)): electrotonic
                        path length of first and last segment on path and
                        their gbar values

    @param  g_bounds    min and max gbar value
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

    @type   gname       str
    @param  gname       full conductance name (including mechanism suffix)

    @type   bound_segs  list(tuple(Segment,Segment)) OR list(tuple(dict,dict)
    @param  bound_segs  pairs of boundary segments

    @type   bound_L     list(tuple(float, float))
    @param  bound_segs  electrotonic lengths of boundary segments

    @return     gbar_interp: the average interpolated gbar over all boundary
                pairs
    """
    # Allow to give a dict instead of segment
    if isinstance(bound_segs[0][0], dict):
        get_gval = lambda seg: seg.get(gname, 0.0)
    else:
        get_gval = lambda seg: get_range_var(seg, gname, 0.0)

    # Interpolate
    gbar_interp = 0.0
    for i, segs in enumerate(bound_segs):
        seg_a, seg_b = segs

        # path length of bounding segments
        L_a, L_b = bound_L[i]
        assert L_b >= L_a, "Path lengths of bounding segments must be in ascending order."

        # Linear interpolation of gbar in seg_a and seg_b according to electrotonic length
        if L_elec <= L_a:
            gbar_interp += get_gval(seg_a)
            continue
        
        if L_elec >= L_b:
            gbar_interp += get_gval(seg_b)
            continue
        
        if L_b == L_a:
            alpha = 0.5
        else:
            alpha = (L_elec - L_a)/(L_b - L_a)
        if alpha > 1.0:
            alpha = 1.0 # if too close to eachother
        
        gbar_a = get_gval(seg_a)
        gbar_b = get_gval(seg_b)
        gbar_interp += gbar_a + alpha * (gbar_b - gbar_a)

    gbar_interp /= len(bound_segs) # take average
    return gbar_interp


def interp_gbar_pick_neighbor(L_elec, gname, bound_segs, bound_L, method='nearest'):
    """ Nearest/left/right neighbor interpolation 

    @param  method  interpolation method: 'nearest'/'left'/'right'

    @see    interp_gbar_linear_neighbors() for other parameters
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
