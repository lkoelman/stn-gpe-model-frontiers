"""
Cell reduction helper functions.


@author Lucas Koelman
@date   03-11-2016
@note   must be run from script directory or .hoc files not found

"""

import numpy as np
from neuron import h

from bgcellmodels.common import logutils
from bgcellmodels.common.nrnutil import (
    getsecref, seg_index, seg_xmin, seg_xmax,
    copy_ion_styles, get_ion_styles, set_ion_styles
)
from bgcellmodels.common.treeutils import subtreeroot
from bgcellmodels.common.electrotonic import seg_lambda

# Load NEURON function libraries
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# logging of DEBUG/INFO/WARNING messages
logger = logutils.getLogger('redops')


################################################################################
# Electrotonic structure
################################################################################


def sec_path_L_elec(secref, f, gleak_name):
    """ Calculate electrotonic path length up to but not including soma
        section (the topmost root section).

    ALGORITHM
    - walk each segment from root section (child of top root) to the given
      section and sum L/lambda for each segment

    @return     tuple pathLelec0, pathLelec1
    @post       pathLelec0 and pathLelec1 are available as attributes on secref

    FIXME: in root node, start walking segments only from midpoint
    """
    rootsec = subtreeroot(secref)
    rootparent = rootsec.parentseg()
    if rootparent is None:
        secref.pathLelec0 = 0.0
        secref.pathLelec1 = 0.0
        return 0.0, 0.0 # if we are soma/topmost root: path length is zero

    # Get path from soma (not including) up to and including this section
    calc_path = h.RangeVarPlot('v')
    calc_path.begin(0.5, sec=rootsec)
    calc_path.end(0.5, sec=secref.sec)
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList

    # Compute electrotonic path length
    secref.pathLelec1 = 0.0 # path length from root sec to 1 end of this sec
    secref.pathLelec0 = 0.0 # path length from root sec to 0 end of this sec
    path_secs = list(root_path)
    path_len = len(path_secs)
    for i, psec in enumerate(path_secs):
        L_seg = psec.L/psec.nseg # segment length
        for seg in psec:
            lamb_seg = seg_lambda(seg, gleak_name, f)
            L_elec = L_seg/lamb_seg
            secref.pathLelec1 += L_elec
            if i < path_len-1:
                secref.pathLelec0 += L_elec

    return secref.pathLelec0, secref.pathLelec1


def sec_path_props(secref, f, gleak_name, linearize_gating=False, init_cell=None):
    """
    Assign path properties to start and end of section, and to all internal segments

    @post   the given SectionRef will have the following attributes:

            - pathL0, pathL1    path length to start and end of section
            - pathL_seg         path length to END of each segment

            - pathri0, pathri1  axial path resistance to start and end of section
            - pathri_seg        axial path resistance to START of each segment

            - pathLelec0            electrotonic path length to START of SECTION
            - pathLelec1            electrotonic path length to END of SECTION
            
            - seg_path_Lelec0       electrotonic path length to START of each SEGMENT
            - seg_path_Lelec1       electrotonic path length to END of each SEGMENT
    """
    rootsec = subtreeroot(secref) # first child section of absolute root

    # Create attribute for storing path resistance
    secref.pathL_seg        = [0.0] * secref.sec.nseg
    secref.pathri_seg       = [0.0] * secref.sec.nseg
    secref.seg_path_Lelec0  = [0.0] * secref.sec.nseg
    secref.seg_path_Lelec1  = [0.0] * secref.sec.nseg
    secref.seg_lambda       = [0.0] * secref.sec.nseg
    # Aliases
    secref.pathL_elec = secref.seg_path_Lelec0
    
    # Initialize path length calculation
    ## Get path from root section to end of given section
    calc_path = h.RangeVarPlot('v')
    calc_path.begin(0.0, sec=rootsec) # x doesn't matter since we only use sections
    calc_path.end(1.0, sec=secref.sec) # x doesn't matter (idem)
    
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList

    # Initialize variables
    path_L = 0.0
    path_L_elec = 0.0
    path_ri = 0.0

    # Compute path length
    path_secs = list(root_path)
    path_len = len(path_secs)
    
    for isec, psec in enumerate(path_secs):
        arrived_sec = (isec==path_len-1) # alternative: use sec.same()
        
        # Start at 0-end of section
        for j_seg, seg in enumerate(psec):
            
            # store path length to start of segment
            if arrived_sec:
                secref.pathri_seg[j_seg] = path_ri
                secref.seg_path_Lelec0[j_seg] = path_L_elec
                if j_seg==0:
                    secref.pathL0 = path_L
                    secref.pathLelec0 = path_L_elec
                    secref.pathri0 = path_ri

            # Update path lengths to end
            seg_L = psec.L/psec.nseg
            seg_lamb = seg_lambda(seg, gleak_name, f)
            seg_L_elec = seg_L/seg_lamb
            path_L += seg_L
            path_L_elec += seg_L_elec
            path_ri += seg.ri()

            # store path length to end of segment
            if arrived_sec:
                secref.pathL_seg[j_seg] = path_L
                secref.seg_path_Lelec1[j_seg] = path_L_elec
                secref.seg_lambda[j_seg] = seg_lamb
                if (j_seg==psec.nseg-1):
                    secref.pathL1 = path_L
                    secref.pathLelec1 = path_L_elec
                    secref.pathri1 = path_ri


def seg_path_L_elec(endseg, f, gleak_name, endpoint):
    """ 
    Calculate electrotonic path length from after root section up to
    start of given segment.

    ALGORITHM
    - walk each segment from root section (child of top root) to the given
      segment and sum L/lambda for each segment

    @param  endpoint: str
            Path endpoint on given segment. Valid values are:
                'segment_start': measure distance to start of segment
                'segment_end': measure distance to end of segment
                'segment_x': measure distance to x-loc of segment

    @return     electrotonic path length
    """
    secref = h.SectionRef(sec=endseg.sec)
    rootsec = subtreeroot(secref)
    j_endseg = seg_index(endseg)
    rootparent = rootsec.parentseg()
    if rootparent is None:
        return 0.0 # if we are soma/topmost root: path length is zero

    # Get path from soma (not including) up to and including this section
    calc_path = h.RangeVarPlot('v')
    calc_path.begin(0.5, sec=rootsec)
    calc_path.end(0.5, sec=secref.sec)
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList

    # Initialize path walk
    path_L_elec = 0.0
    path_secs = list(root_path)
    assert(endseg.sec in path_secs)

    # Compute electrotonic path length
    for psec in path_secs:
        L_seg = psec.L/psec.nseg # segment length

        for j_seg, seg in enumerate(psec):
            # Compute L/lambda for this segment
            lamb_seg = seg_lambda(seg, gleak_name, f)
            seg_L_elec = L_seg/lamb_seg

            # Check if we have reached target segment
            if seg.sec.same(endseg.sec) and j_seg==j_endseg:

                if endpoint == 'segment_start':
                    return path_L_elec
                
                elif endpoint == 'segment_end':
                    return path_L_elec + seg_L_elec
                
                elif endpoint == 'segment_x':
                    partial_L_elec = np.interp(endseg.x,
                        [seg_xmin(endseg, side=None), seg_xmax(endseg, side=None)]
                        [0.0, seg_L_elec])
                    return path_L_elec + partial_L_elec
                else:
                    raise ValueError(endpoint)
                
                return path_L_elec
            
            else:
                path_L_elec += seg_L_elec


    raise Exception('End segment not reached')


def seg_path_L(endseg, endpoint):
    """
    Calculate path length from center of root section to given segment

    @param  endpoint: str
            Path endpoint on given segment. Valid values are:
                'segment_start': measure distance to start of segment
                'segment_end': measure distance to end of segment
                'segment_x': measure distance to x-loc of segment
    """
    secref = h.SectionRef(sec=endseg.sec)
    rootsec = subtreeroot(secref)
    rootparent = rootsec.parentseg()
    if rootparent is None:
        return 0.0 # if we are soma/topmost root: path length is zero
    
    # Get path from root section to endseg
    calc_path = h.RangeVarPlot('v')
    calc_path.begin(0.0, sec=rootsec) # x doesn't matter since we only use path sections
    calc_path.end(endseg.x, sec=endseg.sec)
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList

    # Compute path length
    path_secs = list(root_path)
    path_L = 0.0
    for isec, psec in enumerate(path_secs):
        if not psec.same(endseg.sec):
            path_L += psec.L # full length of Section
        else:
            assert isec == len(path_secs)-1
            if endpoint == 'segment_start':
                return path_L + seg_xmin(endseg, side=None)*psec.L
            elif endpoint == 'segment_end':
                return path_L + seg_xmax(endseg, side=None)*psec.L
            elif endpoint == 'segment_x':
                return path_L + endseg.x*psec.L
            else:
                raise ValueError(endpoint)


class SecProps(object):
    """
    Equivalent properties of merged sections

    NOTE: this is the 'Bunch' recipe from the python cookbook
          al alternative would be `myobj = type('Bunch', (object,), {})()`

    ATTRIBUTES
    ----------

        L           <float> section length
        Ra          <float> axial resistance
        nseg        <int> number of segments
        seg         <list(dict)> RANGE properties for each segment, including
                    'diam' and 'cm'
        children    <list(EqProps)> children attached to 1-end of section
    
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

EqProps = SecProps # alias


def dfs_traversal_iter(start_node):
    """
    Return generator that does depth-first tree traversal.

    @param  start_node : object
            Object with iterable attribute 'children'

    @note   non-recursive, avoids making a new generator per descent
            and avoids blowing up the stack trace
    """
    stack = [start_node]
    while stack:
        node = stack.pop() # LIFO stack
        yield node
        for child in node.children:
            stack.append(child)


def get_sec_range_props(src_sec, mechs_pars):
    """
    Get RANGE properties for each segment in Section.

    @return     list(dict()) of length nseg containing property names and 
                values for each segment in Section
    """
    pnames = [
        par+'_'+mech if mech != '' else par
            for mech, pars in mechs_pars.items() for par in pars
    ]

    seg_props = [{} for i in range(src_sec.nseg)]
    
    # Store segment RANGE properties
    for j_seg, seg in enumerate(src_sec):

        # built-in RANGE properties
        seg_props[j_seg]['diam'] = seg.diam
        seg_props[j_seg]['cm'] = seg.cm
        
        # mechanism RANGE properties
        for pname in pnames:
            if hasattr(seg, pname):
                seg_props[j_seg][pname] = getattr(seg, pname)
    
    return seg_props


def get_sec_props_ref(
        secref, 
        mechs_pars, 
        seg_assigned=None, 
        sec_assigned=None):
    """
    Get both RANGE properties and assigned properties for each segment.

    The properties are stored on an struct-like object named EqProps.


    @pre        all listed properties in seg_assigned and sec_assigned must
                be already computed and stored on the SectionRef object


    @param  secref      SectionRef object to Section for which you want to
                        save the properties

    @param  seg_assigned    list of SectionRef attributes you wish to save
                            (properties per segment)

    @param  sec_assigned    list of SectionRef attributes you wish to save
                            (properties of entire Section)


    @return     EqProps object with requested properties stored as attributes
    """
    if seg_assigned is None:
        seg_assigned = []
    if sec_assigned is None:
        sec_assigned = []

    # Store section properties (non-Range)
    sec_props = EqProps(
                    L=secref.sec.L, 
                    Ra=secref.sec.Ra, 
                    nseg=secref.sec.nseg)
    
    # custom properties
    for prop in sec_assigned:
        setattr(sec_props, prop, getattr(secref, prop))

    # Save membrane mechanisms (density, NOT POINT_PROCESS)
    sec_props.mechanisms = [mech for mech in mechs_pars.iterkeys() if secref.sec.has_membrane(mech)]

    # Initialize segment RANGE properties
    sec_props.seg = [{} for i in range(secref.sec.nseg)]
    bprops = [
        par+'_'+mech if mech != '' else par
            for mech, pars in mechs_pars.items() for par in pars
    ]
    
    # Store segment RANGE properties
    for j_seg, seg in enumerate(secref.sec):
        
        # Store NEURON properties
        for prop in bprops:
            if hasattr(seg, prop):
                sec_props.seg[j_seg][prop] = getattr(seg, prop)
        
        # Store CUSTOM properties (stored on SectionRef)
        for prop in seg_assigned:
            sec_props.seg[j_seg][prop] = getattr(secref, prop)[j_seg]
    
    return sec_props

# Aliases
get_sec_props_obj = get_sec_props_ref


def get_sec_props(sec, mechs_pars):
    """
    Get Section properties and save in new SecProps object.

    @param  sec             <nrn.Section> NEURON section

    @param  mechs_pars      dict mechanism_name -> [parameter_names] with segment 
                            properties to save. To include diam and cm, use a
                            key-value pair {'' : ['diam', 'cm']}

    @return                 <SecProps> object
    """
    # Store section properties (non-Range)
    sec_props = SecProps(
                    L=sec.L, 
                    Ra=sec.Ra,
                    nseg=sec.nseg)

    # Initialize dicts with RANGE properties
    sec_props.seg = [dict() for i in range(sec.nseg)]
    parnames = [
        par+'_'+mech if mech != '' else par
            for mech, pars in mechs_pars.items() for par in pars
    ]
    
    # Store segment RANGE properties
    for j_seg, seg in enumerate(sec):
        for pname in parnames:
            sec_props.seg[j_seg][pname] = getattr(seg, pname)

    return sec_props


def merge_sec_properties(
        src_props,
        tar_sec,
        mechs_pars,
        ion_styles=True,
        check_uniform=True
    ):
    """
    Merge section properties from multiple SecProps objects into target Section.

    WARNING: properties are queried at center segment and assigned to the
    target section like sec.ppty = val (i.e. not on a per-segment basis)

    @param  src_props   list(SecProps), each with attribute 'seg' containing one dict
                        of segment attributes per segment

    @param  mechs_pars  dict mechanism_name -> [parameter_names] with segment 
                        properties (RANGE properties) to copy
    """

    # keep track of assigned parameters
    assigned_params = {
        par+'_'+mech if mech != '' else par: None
            for mech,pars in mechs_pars.items() for par in pars
    }

    # merge all mechanism RANGE properties
    for src_sec in src_props:
        nseg = src_sec.nseg
        segs = src_sec.seg
        for mechname, parnames in mechs_pars.items():
            
            # Check if any segment in source has the mechanism
            if mechname in src_sec.mechanisms:
                # if any((p.endswith(mechname) for i in range(nseg) for p in segs[i].keys())):
                tar_sec.insert(mechname)
            else:
                continue
            
            # Copy parameter values
            for pname in parnames:
                fullpname = pname+'_'+mechname if mechname != '' else pname

                # Check that parameter values are uniform in section
                mid_val = segs[int(nseg)/2].get(fullpname, None)
                
                if check_uniform and any(
                    (mid_val!=segs[i].get(fullpname, None) for i in range(nseg))):
                    raise Exception("Parameter {} is not uniform in source section {}."
                                    "Cannot assign non-uniform value to Section.".format(
                                        fullpname, src_sec))

                # Check that parameter value is same as in other source sections
                prev_val = assigned_params[fullpname]
                
                if check_uniform and (prev_val is not None) and mid_val!=prev_val:
                    raise Exception("Parameter {} is not uniform between source sections."
                                    "Cannot assign non-uniform value to Section.".format(
                                        fullpname))

                # Copy value to entire section
                try:
                    setattr(tar_sec, fullpname, mid_val)
                    assigned_params[fullpname] = mid_val
                except AttributeError:
                    pass # setattr fails if property is a GLOBAL PARAM (default granularity for PARAM)

    if ion_styles:
        merge_ion_styles(src_props, tar_sec, check_uniform)


def merge_ion_styles(src_props, tar_sec, check_uniform=True):
    """
    Set correct ion styles for each Section.

    @param          src_props : iterable(SecProps)
                    Object containing Section properties

    @effect         By default, ion styles of merged Sections are looked up
                    in the saved original tree properties and copied to the
                    sections (if they are all the same). This method may be
                    overridden.
    """
    # Look up ion styles info of sections that have been merged into this one
    ionstyles_dicts = [p.ion_styles for p in src_props]

    # They must all be the same or we have a problem with mechanisms
    final_styles = ionstyles_dicts[0]
    for styles_dict in ionstyles_dicts[1:]:
        for ion, style_flags in styles_dict.items():
            
            if ion not in final_styles:
                final_styles[ion] = style_flags
            
            elif final_styles[ion] != style_flags:
                raise ValueError("Cannot merge Sections with distinct ion styles!")

    # Copy styles to target Section
    set_ion_styles(tar_sec, **final_styles)


def copy_sec_properties(src_sec, tar_sec, mechs_pars):
    """
    Copy section properties from source to target Section
    """
    # Number of segments and mechanisms
    tar_sec.nseg = src_sec.nseg
    for mech in mechs_pars.keys():
        if hasattr(src_sec(0.5), mech):
            tar_sec.insert(mech)

    # Geometry and passive properties
    tar_sec.L = src_sec.L
    tar_sec.Ra = src_sec.Ra

    # copy RANGE properties
    for seg in src_sec:
        tar_sec(seg.x).diam = seg.diam # diameter
        tar_sec(seg.x).cm = seg.cm # capacitance
        for mech in mechs_pars.keys():
            for par in mechs_pars[mech]:
                prop = par+'_'+mech
                setattr(tar_sec(seg.x), prop, getattr(seg, prop))

    # ion styles
    copy_ion_styles(src_sec, tar_sec)


def save_tree_properties(node_sec, mechs_pars):
    """
    Save properties of all Sections in subtree of given Section

    @param      mechs_pars      dict mechanism_name -> [parameter_names] with segment 
                                properties to save

    @param      save_ion_styles list(str) containing ion names: ion styles you want to save
    
    @return     EqProps tree with requested properties stored as attributes
    """
    # Create SecProps object for current node
    sec_props = get_sec_props(node_sec, mechs_pars)

    # Call for each child and add to children
    sec_props.children = []
    for child_sec in node_sec.children(): # Depth-first tree traversal
        sec_props.children.append(save_tree_properties(child_sec, mechs_pars))

    return sec_props


def save_tree_properties_ref(
        node_ref,
        all_refs,
        mechs_pars,
        sec_assigned_props=None,
        seg_assigned_props=None,
        save_ion_styles=True
    ):
    """
    Save properties of all Sections in subtree of given SectionRef

    @param      mechs_pars      dict mechanism_name -> [parameter_names] with segment 
                                properties to save

    @param      save_ion_styles: bool or list(str)
                ion styles you want to save
    
    @return     EqProps tree with requested properties stored as attributes
    """
    # Process current node: create SecProps object
    sec_props = get_sec_props_ref(
                    node_ref,
                    mechs_pars,
                    sec_assigned=sec_assigned_props,
                    seg_assigned=seg_assigned_props)

    if save_ion_styles is not None:
        sec_props.ion_styles = get_ion_styles(node_ref.sec, ions=save_ion_styles)

    # Process child nodes
    sec_props.parent = None
    sec_props.children = []
    for child_sec in node_ref.sec.children(): # Depth-first tree traversal
        child_ref = getsecref(child_sec, all_refs)
        child_props = save_tree_properties_ref(
                        child_ref, all_refs, mechs_pars, 
                        sec_assigned_props=sec_assigned_props,
                        seg_assigned_props=seg_assigned_props,
                        save_ion_styles=save_ion_styles)
        child_props.parent = sec_props
        sec_props.children.append(child_props)

    return sec_props


def subtree_assign_attributes(noderef, allsecrefs, attr_dict):
    """
    Assign attributes to all SectionRefs in subtree of given section

    @param attr_dict    dictionary of key-value pairs (attribute_name, attribute_value)
    """
    # Process self
    for aname, aval in attr_dict.items():
        setattr(noderef, aname, aval)

    # Process children
    childsecs = noderef.sec.children()
    childrefs = [getsecref(sec, allsecrefs) for sec in childsecs]
    for childref in childrefs:
        subtree_assign_attributes(childref, allsecrefs, attr_dict)


def subtree_assign_gids_dfs(node_ref, all_refs, parent_id=0):
    """
    Label nodes with int identifier by doing depth-first tree traversal
    and increment the id upon each node visit.
    """
    if node_ref is None:
        return

    highest = node_ref.gid = parent_id + 1

    childsecs = node_ref.sec.children()
    childrefs = [getsecref(sec, all_refs) for sec in childsecs]
    for childref in childrefs:
        highest = subtree_assign_gids_dfs(childref, all_refs, highest)

    return highest


def find_secprops(node, filter_fun, find_all=True):
    """
    Find SecProps object in tree satisfying given filter function

    @param  find_all    True if all nodes satisfying filter function should
                        be returned

    @return             list(SecProps) that match the filter function (may be
                        empty)
    """
    nodes_gen = dfs_traversal_iter(node)

    if find_all:
        return list(filter(filter_fun, nodes_gen))
    else:
        try:
            match = next((n for n in nodes_gen if filter_fun(n)))
            return [match]
        except StopIteration:
            return []


def find_roots_at_level(level, all_refs):
    """
    Find all root Sections art a particular level.

    @param  level       (int) level of branch points to return

    @pre                ref.level and ref.end_branchpoint attributes have been set, 
                        e.g. via a call to cluster.assign_topology_attrs().
    """
    return [ref for ref in all_refs if (ref.level==level and ref.end_branchpoint)]


def find_outer_branchpoints(node_ref, inc_refs):
    """
    @param      node_ref: SectionRef
                Section to start tree traversal

    @param      inc_refs: list(SectionRef)
                Sections included in search, other sections are skipped
    
    @return     list(SectionRef)
    """

    # Check if current cylinder ends in branch point (fork)
    child_refs = [getsecref(sec, inc_refs) for sec in node_ref.sec.children()]
    child_refs = [ch for ch in child_refs if ch is not None] # getsecref returns None if sec not in all_refs
    if len(child_refs) == 0:
        return []
    
    # Else, collect child branch points
    child_branchpoints = sum((find_outer_branchpoints(ch, inc_refs) for ch in child_refs), [])

    if len(child_branchpoints) > 0:
        return child_branchpoints
    elif node_ref.end_branchpoint:
        return [node_ref]
    else:
        return []


if __name__ == '__main__':
    print("__main__ @ redutils.py")
