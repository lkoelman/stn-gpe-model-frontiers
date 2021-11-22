"""
Utilities for dealing with NEURON cell models

@author Lucas Koelman
"""

from neuron import h
import numpy as np

from .nrnutil import seg_index, seg_xmin, seg_xmax, seg_xmid, seg_at_index, getsecref
from io import StringIO

# aliases to avoid repeatedly doing multiple hash-table lookups
_h_section_ref = h.SectionRef
_h_allsec = h.allsec
_h_parent_connection = h.parent_connection
_h_section_orientation = h.section_orientation

def parent(sec):
    """
    Return the parent of sec or None if sec is a root

    Copied from NEURON RxD module

    @author R.A McDougal
    """
    sref = _h_section_ref(sec=sec)
    if sref.has_trueparent():
        return sref.trueparent().sec
    elif sref.has_parent():
        temp = sref.parent().sec
        # check if temp owns the connection point
        if _h_section_ref(sec=temp).has_parent() and _h_parent_connection(sec=temp) == _h_section_orientation(sec=temp):
            # connection point belongs to temp's ancestor
            return parent(temp)
        return temp
    else:
        return None


def parent_loc(sec, trueparent):
    """
    Return the position on the (true) parent where sec is connected

    Note that _h_section_orientation(sec=sec) is which end of the section is
    connected.

    Copied from NEURON RxD module

    @author R.A McDougal
    """
    # TODO: would I ever have a parent but not a trueparent (or vice versa)
    sref = _h_section_ref(sec=sec)
    parent = sref.parent().sec
    while parent != trueparent:
        sec, parent = parent, _h_section_ref(sec=sec).parent().sec
    return _h_parent_connection(sec=sec)


def prev_seg(curseg, x_loc='mid'):
    """
    Get segment preceding seg: this can be on same or parent Section
    """
    # prevseg = [seg for seg in curseg.sec if seg_index(seg) < seg_index(curseg)][-1]
    i_seg = seg_index(curseg)
    if i_seg > 0:
        prevseg = seg_at_index(curseg.sec, i_seg-1)
    else:
        parseg = curseg.sec.parentseg()
        if parseg is None:
            return None
        parsec = parseg.sec
        # Do not return the 0-area nodes at 0-end and 1-end
        prevseg = parsec(seg_xmid(parseg))

    # Adjust x-loc if required
    if x_loc == 'mid':
        return prevseg
    elif x_loc == 'min':
        prevseg = prevseg.sec(seg_xmin(prevseg, side='inside'))
    elif x_loc == 'max':
        prevseg = prevseg.sec(seg_xmax(prevseg, side='inside'))
    else:
        raise ValueError("Invalid value {} for argument 'x-loc'".format(x_loc))

    return prevseg

def next_segs(curseg, x_loc='mid'):
    """
    Get child segments of given segment

    @param  x_loc: str

            Which x-value to associate with each returned segment:
            - 'mid': center of the segment
            - 'min': minimum x-value inside the segment
            - 'max': maximum x-value inside the segment
    """
    cursec = curseg.sec
    i_rootseg = seg_index(curseg)
    child_segs = []

    # find next segments
    if i_rootseg < cursec.nseg-1:
        # not end segment -> return next segment of same Section
        child_segs.append(next((seg for seg in cursec if seg_index(seg)>i_rootseg)))
    else:
        # end segment -> return first segment of each child section
        child_segs = [next((seg for seg in sec)) for sec in cursec.children()]

    # Adjust x-loc if required
    if x_loc == 'min':
        child_segs = [seg.sec(seg_xmin(seg, side='inside')) for seg in child_segs]
    elif x_loc == 'max':
        child_segs = [seg.sec(seg_xmax(seg, side='inside')) for seg in child_segs]
    elif x_loc != 'mid':
        raise ValueError("Invalid value {} for argument 'x-loc'".format(x_loc))

    return child_segs


def next_segs_dx(curseg, dx):
    """
    Get next segments (in direction 0 to 1 end) with discretization step x
    """
    x = curseg.x

    if x + dx <= 1.0:
        return [curseg.sec(x + dx)]
    else:
        xb = x + dx - 1.0
        return [sec(xb) for sec in curseg.sec.children()]


def next_segs_dL(curseg, dL):
    """
    Get next segments (in direction 0 to 1 end) with step size dL [um]
    """
    aL = curseg.x * curseg.sec.L
    bL = aL + dL

    if bL <= curseg.sec.L:
        return [curseg.sec(bL/curseg.sec.L)]
    else:
        sL = bL - curseg.sec.L
        return [sec(sL/sec.L) for sec in curseg.sec.children()]


def interp_seg(seg, a, b):
    """
    Interpolate values [a,b] at segment boundaries using segment's x-loc
    """
    seg_dx = 1.0/seg.sec.nseg
    iseg = min(int(seg.x/seg_dx), seg.sec.nseg-1)
    x_a = iseg * seg_dx
    return a + (seg.x - x_a)/seg_dx * (b-a)



def sameparent(secrefA, secrefB):
    """ Check if sections have same parent section """
    if not (secrefA.has_parent() and secrefB.has_parent()):
        return False
    apar = secrefA.parent # changes CAS
    bpar = secrefB.parent # changes CAS
    h.pop_section()
    h.pop_section()
    return apar.same(bpar)


def subtreeroot(secref):
    """
    Find the root section of the tree that given sections belongs to.
    I.e. the first section after the root of the entire cell.
    """
    # Get root section of tree
    orig = secref.root # changes the cas
    h.pop_section()

    for root in orig.children():
        # Get subtree of the current root
        roottree = h.SectionList()

        # Fill SectionList with subtree of CAS
        roottree.subtree(sec=root)

        # Check if given section in in subtree
        if secref.sec in roottree:
            return root

    return orig


def wholetreeroot(secref, allsecrefs):
    """
    Find absolute root of tree.

    @return     nrn.Section
                absolute root of tree
    """
    root_sec = allsecrefs[0].root; h.pop_section() # pushes CAS
    return root_sec


def subtree_secs(rootsec):
    """
    Get all Sections in subtree of but not including rootsec.
    """
    tree_secs = h.SectionList()
    tree_secs.subtree(sec=rootsec) # includes rootsec itself

    return [sec for sec in tree_secs if not rootsec.same(sec)]


def wholetree_secs(sec):
    """
    Get all Sections in the same cell (i.e. that have a path to the given section)
    """
    tree_secs = h.SectionList()
    tree_secs.wholetree(sec=sec)

    return list(tree_secs)


def root_sections():
    """
    Returns a list of all sections that have no parent.

    @author     Alex H Williams
    @citation   DOI:10.5281/zenodo.12576
    """
    roots = []
    for section in h.allsec():
        sref = h.SectionRef(sec=section)
        if not sref.has_parent():
            roots.append(section)
    return roots


def root_section(tree_sec):
    """
    Return the root section for tree that given Section is part of.
    """
    ref = h.SectionRef(sec=tree_sec)
    root = ref.root
    h.pop_section()
    return root


def leaf_sections(root_sec=None, subtree=False):
    """
    Returns a list of all sections that have no children.

    @param  tree_sec : Hoc.Section
            Any section that is part of the tree of which you want to get
            the leaves.
    """
    leaves = []
    if root_sec is None:
        all_sec = h.allsec() # iterator
    elif subtree:
        all_sec = h.SectionList()
        all_sec.subtree(sec=root_sec) # includes rootsec itself
    else:
        all_sec = wholetree_secs(root_sec)

    for sec in all_sec:
        ref = h.SectionRef(sec=sec)
        if ref.nchild() < 0.9:
            leaves.append(sec)

    return leaves


def path_sections(target, source=None):
    """
    Get all sections on path from source to target Section.

    @param  source : nrn.Section
            If source is None, it will be set to root section of tree

    @return path : list(nrn.Section)
            List of Section on path between source and target
    """
    if source is None:
        ref_target = h.SectionRef(sec=target)
        source = ref_target.root; h.pop_section()

    # Get path from soma (not including) up to and including this section
    calc_path = h.RangeVarPlot('v')
    source.push()
    calc_path.begin(0.5)
    target.push()
    calc_path.end(0.5)

    path_secs = h.SectionList()
    calc_path.list(path_secs) # copy path sections to SectionList
    h.pop_section()
    h.pop_section()

    return list(path_secs)


def path_segments(target, source=None):
    """
    Get all segments on path from source to target segment.

    @param  source : nrn.Segment
            If source is None, it will be set to the center of the root
            Section of the tree

    @return path : list(nrn.Segment)
            List of segments on path between source and target
    """

    ref_target = h.SectionRef(sec=target)
    root = ref_target.root(0.5); h.pop_section()

    if source is None:
        source = root

    # get all Sections between source and target segment
    calc_path = h.RangeVarPlot('v')
    source.push()
    calc_path.begin(0.5)
    target.push()
    calc_path.end(0.5)
    path_secs = h.SectionList()
    calc_path.list(path_secs) # copy path sections to SectionList
    h.pop_section()
    h.pop_section()
    path_secs = list(path_secs)

    # All calls to distance() are relative to target segment
    h.distance(0, target.x, sec=target.sec)

    # Check second section distance to see if we're walking toward or away from soma
    path_segments = []
    last_dist = h.distance(source.x, sec=source.sec)
    for i_sec, sec in enumerate(path_secs):

        # if in different section than target -> can add any segment with dist < start, and add them in order of decreasing distance

        # TODO: if in same section as target -> find end (0/1) that is closest to last segment (of previous section), then add all segments with x between that end_segment and target_segment

        # Distance must always become smaller
        current_segments = []
        current_distances = []
        for seg in sec: # iterate in 0-to-1 order
            distance = h.distance(seg.x, sec=sec)
            if distance <= last_dist:
                current_segments.append(seg)
                current_distances.append(distance)
                last_dist = distance

        # Segments must be added in order of decreasing distance to target
        if current_distances[0] > current_distances[-1]:
            path_segments.extend(current_segments)
        else:
            path_segments.extend(current_segments[::-1])


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


def dfs_iter_tree_recursive(node):
    """
    Return generator that does depth-first tree traversal starting at given node.

    @note   makes new generator for each child, not as elegant
    """
    yield node

    for child_node in node.children():
        for cn in dfs_iter_tree_recursive(child_node):
            yield cn


def ascend_sectionwise_dfs(start_node):
    """
    Return generator that does depth-first tree traversal starting at given node.

    @note   non-recursive, avoids making a new generator per descent
            and avoids blowing up the stack trace
    """
    stack = [start_node]
    while stack:
        node = stack.pop() # LIFO stack
        yield node
        for child in node.children():
            stack.append(child)


def ascend_sectionwise_bfs(start_node):
    """
    Return generator that does breadth-first tree traversal starting at given node.

    NOTE: this generates sections in the same order as NEURON functions
          h.SectionList().wholetree(sec=root) and 
          h.SectionList().wholetree(sec=root)

    @note   non-recursive, avoids making a new generator per descent
            and avoids blowing up the stack trace
    """
    queue = [start_node]
    while queue:
        node = queue.pop(0) # FIFO queue
        yield node
        for child in node.children():
            queue.append(child)


def ascend_segmentwise_dfs(start_section):
    """
    Return generator that does depth-first tree traversal of segments
    starting at 0-end of given section.

    @note   non-recursive, avoids making a new generator per descent
            and avoids blowing up the stack trace
    """
    stack = [start_section]
    while stack:
        section = stack.pop()
        for segment in section:
            yield segment
        for child_section in section.children():
            stack.append(child_section)


def ascend_with_fixed_spacing(start_segment, dL):
    """
    Generator that does depth-first tree traversal of segments sampled
    at a fixed spacing from each other.

    @param      start_segment : nrn.Segment
                Segment to start the descent.

    @effect     If start_segment is in the root Section of the tree,
                all subtrees will be ascended no matter whether they
                are attached to the 0-end or 1-end. If start_segment
                is not in the root_section, ascent will proceed
                in the 1-direction only.
    """
    assert dL > 0, "dL must be strictly positive"

    stack = [start_segment]
    root_sec = root_section(start_segment.sec)
    start_seg_x = start_segment.x

    while stack:
        curseg = stack.pop()
        yield curseg

        stepping_backward = curseg.sec.same(root_sec) and curseg.x < start_seg_x
        both_directions =  curseg.sec.same(root_sec) and curseg.x == start_seg_x

        La = curseg.x * curseg.sec.L    # distance from 0-end in micron

        if (not stepping_backward) or both_directions:
            Lb = La + dL                # distance from 0-end of next sample

            if Lb <= curseg.sec.L:      # next sample in same Section
                stack.append(curseg.sec(Lb/curseg.sec.L))
            else:
                Lrem = Lb - curseg.sec.L
                stack.extend([sec(Lrem/sec.L) for sec in curseg.sec.children() if sec.parentseg().x!=0])

        if stepping_backward or both_directions:
            # Case where we are walking in reverse direction (1-end to 0-end)
            Lb = La - dL
            if Lb >= 0:
                stack.append(curseg.sec(Lb/curseg.sec.L))
            else:
                Lrem = abs(Lb)
                stack.extend([sec(Lrem/sec.L) for sec in curseg.sec.children() if sec.parentseg().x==0])


def sample_tree_uniformly(root_sec, num_seg, spacing,
                          filter_func=None, rng=None, replace=False):
    """
    Sample dendritic tree with uniform spacing (in micron).

    @param      root_sec : nrn.Section
                Root section of the tree

    @param      num_seg : int
                Number of sections to pick

    @param      spacing : float
                Spacing between segments in micron

    @param      filter_func : callable(nrn.Segment) -> bool (optional)
                Filter function to mark segments as eligible to be picked

    @param      replace : bool
                Sample with replacement

    @param      rng : numpy.Random (optional)
                Random number generator.
    """
    # Get random number generator
    if rng is None:
        rng = np.random

    # Generate sample segments with uniform spacing
    segment_sampler = ascend_with_fixed_spacing(root_sec(0.5), spacing)
    if filter_func is not None:
        eligible_segs = [seg for seg in segment_sampler if filter_func(seg)]
    else:
        eligible_segs = [seg for seg in segment_sampler]
    sample_indices = rng.choice(len(eligible_segs), num_seg, replace=replace)
    return [eligible_segs[i] for i in sample_indices]


def subtree_topology(sub_root, max_depth=1e9):
    """
    Like h.topology() but for subtree of section.

    @see    nrnhome/nrn/src/nrnoc/solve.c::nrnhoc_topology()
    """

    buff = StringIO.StringIO()

    def dashes(sec, offset, lead_char, dist=0):
        """
        @param  dist    distance from first section (subtree root)
        """

        orient = int(h.section_orientation(sec=sec))
        direc = "({}-{})".format(orient, 1-orient)

        # Print section in format -----| with one dash per segment
        buff.write(" " * offset)
        buff.write(lead_char)
        if dist == max_depth+1:
            # truncate end sections
            num_child = len(subtree_secs(sec))
            buff.write("-..       %s%s + %i children\n" % (sec.name(), direc, num_child))
            return # end of recursion
        else:
            buff.write("-" * sec.nseg)

        # Print termination symbol and section description
        buff.write("|       %s%s\n" % (sec.name(), direc))

        for child_sec in reversed(sec.children()): # reversed since NEURON uses stack + pop
            # get index of segment where child connects to parent
            con_seg_idx = seg_index(child_sec.parentseg())
            # buff.write(" ")
            dashes(child_sec, con_seg_idx+offset+3, "`", dist+1)


    dashes(sub_root, 0, ".|")

    buff_string = buff.getvalue()
    buff.close()
    return buff_string


def check_tree_constraints(sections, raise_any=False):
    """
    Check unbranched cable assumption and orientation constrained.

    @param      raise_any : bool
                Raise exception if any inconsistency is encountered

    @return     a tuple (unbranched, oriented, branched, misoriented) of type
                tuple(bool, bool, list(Section), list(Section))
                with following entries:

                    bool: all Sections are unbranched

                    bool: all sections are correctly oriented, i.e. the 0-end is
                          connected to the 1-end of the parent, except if the parent
                          is the root section in which case a connection to the 0-end
                          is permitted.

                    list: all branched sections

                    list: all misoriented sections
    """
    # check both connect(child(x), parent(y))
    # parent_y: for all sections in whole tree: sec.parentseg().x must be 1.0
    #           EXCEPT for sections connected to root
    # child_x:  for all sections in whole tree: sec.orientation() (h.section_orientation(sec=sec)) must 0.0

    is_unbranched = True
    is_oriented = True
    branched = set()
    misoriented = set()

    first_ref = h.SectionRef(sec=sections[0])
    tree_root = first_ref.root; h.pop_section() # pushes CAS

    for sec in sections:

        parent_sec = sec.parentseg().sec
        parent_y = sec.parentseg().x
        orient_parent_ok = parent_y==1.0 or (parent_y==0.0 and parent_sec.same(tree_root))
        branch_parent_ok = parent_y==1.0 or parent_y==0.0

        self_x = sec.orientation() # see h.section_orientation()
        orient_self_ok = self_x==0.0
        branch_self_ok = self_x==0.0 or self_x==1.0

        is_unbranched = is_unbranched and branch_parent_ok and branch_self_ok
        is_oriented = is_oriented and orient_parent_ok and orient_self_ok

        if raise_any and not (orient_parent_ok and orient_self_ok and
                              branch_parent_ok and branch_self_ok):
            raise Exception("Branching topology constraints not met at "
                            "section {}".format(sec.name()))


        if not orient_parent_ok:
            misoriented.update(parent_sec)
        if not orient_self_ok:
            misoriented.update(sec)

        if not branch_parent_ok:
            branched.update(parent_sec)
        if not branch_self_ok:
            branched.update(sec)

    return is_unbranched, is_oriented, list(branched), list(misoriented)


def gather_rangevar_along_paths(
        root,
        leaves=None, 
        measure_funcs=None,
        rangevar_names=None):
    """
    Measure electrotonic properties along paths.
    
    @param  root : nrn.Section
            Root section for measurement

    @param  leaves : iterable[nrn.Section]
            Endpoints for measurements

    @param  measure_funcs : function(nrn.Segment) -> any

    @return path_measurements : list[dict[measure, list[float]]]
            For each leaf node, a dictionary containing all the measurements
            as a function of distance from the root node to the leaf node
    """
    # Get leaf sections
    if leaves is None:
        leaves = leaf_sections(root, subtree=False) # all leaves of tree
    if measure_funcs is None:
        measure_funcs = {}
    if rangevar_names is None:
        rangevar_names = []

    # make paths to leaves and measure each measure
    leaf_path_measurements = [] # one dict for each leaf

    # Get path lengths
    h.distance(0, 0.5, sec=root)
    path_dist_um = lambda seg: h.distance(seg.x, sec=seg.sec)

    # Get electrotonic measurements
    for leaf_sec in leaves:

        path_measurements = {}

        path_segments = [seg for sec in path_sections(leaf_sec) for seg in sec]

        # Get path lengths
        path_measurements['path_dist_um'] = map(path_dist_um, path_segments)

        # Collect rangevar names
        for rangevar_name in rangevar_names:
            # has_mechanism = lambda seg, mech: h.ismembrane(rangevar_name.split('_')[-1], sec=seg.sec) if '_' in mech else True
            # measure_func = lambda seg: getattr(seg, rangevar_name) if h.ismembrane(rangevar_name.split('_')[-1], sec=seg.sec) else 0.0
            path_measurements[rangevar_name] = map(
                lambda seg: getattr(seg, rangevar_name),
                path_segments)

        # Measure using custom functions
        for measure_name, measure_func in measure_funcs.items():
            path_measurements[measure_name] = map(measure_func, path_segments)

        leaf_path_measurements.append(path_measurements)
        
    return leaf_path_measurements


def plot_rangevar_per_section(
        root_sec,
        propname,
        secfilter=None,
        labelfunc=None,
        y_range=None,
        fig=None,
        ax=None,
        plt=None):
    """
    Plot neuron RANGE variable in each section, with one figure axis
    per secton.


    @param  root_sec : nrn.Section
            Root section to start tree ascent

    @param  propname : str
            RANGE var name

    @param  secfilter : function(nrn.Section) -> bool
            filter function applied to each Section
    """
    if root_sec is None:
        return

    # Default arguments
    if secfilter is None:
        secfilter = lambda sec: True

    # Set up plotting
    first_call = fig is None
    if first_call:
        fig = plt.figure()
    if y_range is None:
        y_range = [float('inf'), float('-inf')]

    # Plot current node
    if secfilter(root_sec):
        if not ax:
            nax = len(fig.axes)
            for i, ax in enumerate(fig.axes):
                ax.change_geometry(nax+1, 1, i+1) # change grid and position in grid
            ax = fig.add_subplot(nax+1, 1, nax+1)

        # get data to plot
        sec = root_sec
        xg = [(seg.x, getattr(seg, propname)) for seg in sec]
        xvals, gvals = zip(*xg)
        if min(gvals) < y_range[0]:
            y_range[0] = min(gvals)
        if max(gvals) > y_range[1]:
            y_range[1] = max(gvals)

        # plot it
        ax.plot(xvals, gvals, 'o')
        ax.plot(xvals, gvals, '-')
        # ax.set_xlabel('x')
        ax.set_ylabel(labelfunc(sec))

        if y_range is not None:
            ax.set_ylim(y_range)

    # plot children
    for child_sec in sec.children():
        plot_rangevar_per_section(child_sec, propname, secfilter, labelfunc, y_range, fig)

    if first_call:
        plt.suptitle(propname)
        y_span = y_range[1]-y_range[0]
        for ax in fig.axes:
            ax.set_ylim((y_range[0]-0.1*y_span, y_range[1]+0.1*y_span))
        plt.show(block=False)
    
    return fig


# seg_builtin_attrs = ['area', 'cm', 'diam', 'hh', 'na_ion', 'k_ion', 'ca_ion', 'next', 'node_index', 'point_processes', 'ri', 'sec', 'v', 'x']

# pp_builtin_attrs = ['allsec', 'baseattr', 'cas', 'get_loc', 'has_loc', 'loc', 'hname', 'hocobjptr', 'next', 'ref', 'same', 'setpointer', 'Section']


# def print_pp_info(cell_sec=None, mechs_params=None):
#   """
#   Print information about all POINT_PROCESS mechanisms.

#   WARNING: this doesn't seem to work in Python ?!
#   """
#   if cell_sec:
#       all_cell_secs = wholetree_secs(cell_sec)
#   else:
#       all_cell_secs = [sec for sec in h.allsec()]

#   for sec in all_cell_secs:
#       for seg in sec:
#           for pp in seg.point_processes():

#               print '\nInfo for point process {} @ {}'.format(pp, seg)

#               mech_name = get_mod_name(pp)
#               if mechs_params and mech_name in mechs_params:
#                   param_names = mechs_params[mech_name]
#               else:
#                   param_names = [attr for attr in dir(pp) if (not attr.startswith('__') and (attr not in pp_builtin_attrs))]

#               for param_name in param_names:
#                   print '{} : {}'.format(param_name, getattr(pp, param_name, 'not found'))


