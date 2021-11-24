"""
Cell reduction helper functions.


@author Lucas Koelman
@date   03-11-2016
@note   must be run from script directory or .hoc files not found

"""

import numpy as np
import matplotlib.pyplot as plt
import math

import logging
logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logname = "redops" # __name__
logger = logging.getLogger(logname) # create logger for this module

import neuron
from neuron import h

from bgcellmodels.common.treeutils import *
from bgcellmodels.common.electrotonic import *

# Load NEURON function libraries
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

################################################################################
# Electrotonic structure
################################################################################

def sec_path_ri(secref, store_seg_ri=False):
    """ Calculate axial path resistance from root to 0 and 1 end of each section

    @effect     calculate axial path resistance from root to 0/1 end of sections
                and set as properties pathri0/pathri1 on secref

    @return     tuple pathri0, pathri1
    """
    # Create attribute for storing path resistance
    if store_seg_ri and not(hasattr(secref, 'pathri_seg')):
        secref.pathri_seg = [0.0] * secref.sec.nseg

    # Get subtree root
    rootsec = subtreeroot(secref)
    rootparent = rootsec.parentseg()

    # Calculate path Ri
    if rootparent is None: # if we are soma/top root: path length is zero
        secref.pathri0 = 0.
        secref.pathri1 = 0.
        return 0., 0.
    
    # Get path from root node to this sections
    calc_path = h.RangeVarPlot('v')
    rootsec.push()
    calc_path.begin(0.5)
    secref.sec.push()
    calc_path.end(0.5)
    root_path = h.SectionList()
    calc_path.list(root_path) # store path in list
    h.pop_section()
    h.pop_section()

    # Compute axial path resistances
    pathri0 = 0. # axial path resistance from root to start of target Section
    pathri1 = 0. # axial path resistance from root to end of target Section
    path_secs = list(root_path)
    path_len = len(path_secs)
    for isec, psec in enumerate(path_secs):
        arrived = bool(psec.same(secref.sec))
        for jseg, seg in enumerate(psec):
            # Axial path resistance to start of current segment
            if store_seg_ri and arrived:
                secref.pathri_seg[jseg] = pathri1
            # Axial path resistance to end of current segment
            pathri1 += seg.ri()
            # Axial path resistance to start of target section
            if isec < path_len-1:
                pathri0 = pathri1

    secref.pathri0 = pathri0
    secref.pathri1 = pathri1
    return pathri0, pathri1

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
    rootsec.push()
    calc_path.begin(0.5)
    secref.sec.push()
    calc_path.end(0.5)
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList
    h.pop_section()
    h.pop_section()

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

def sec_path_L(secref):
    """
    Assign path length to end of all segments and to 0 and 1-end of Section

    @param to_end   if True, return distance to end of segment, else
                    return distance to start of segment
    """
    rootsec = subtreeroot(secref)
    j_endseg = seg_index(endseg)
    rootparent = rootsec.parentseg()
    if rootparent is None:
        return 0.0 # if we are soma/topmost root: path length is zero

    # Create attribute for storing path resistance
    secref.pathL_seg = [0.0] * secref.sec.nseg
    
    # Get path from root section to endseg
    calc_path = h.RangeVarPlot('v')
    rootsec.push()
    calc_path.begin(0.0) # x doesn't matter since we only use path sections
    endseg.sec.push()
    calc_path.end(endseg.x)
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList
    h.pop_section()
    h.pop_section()

    # Compute path length
    path_secs = list(root_path)
    path_L = 0.0
    for isec, psec in enumerate(path_secs):
        arrived_sec = bool(psec.same(endseg.sec))
        seg_L = psec.L/psec.nseg
        for j_seg, seg in enumerate(psec):
            if arrived_sec:
                if j_seg==0:
                    secref.pathL0 = path_L
                elif (j_seg==j_endseg):
                    secref.pathL1 = path_L + seg_L
                    return secref.pathL1
            path_L += seg_L
            secref.pathL_seg[j_seg] = path_L


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
    rootsec = subtreeroot(secref)
    rootparent = rootsec.parentseg()
    if rootparent is None:
        return # if we are soma/topmost root: path length is zero

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
    
    rootsec.push()
    calc_path.begin(0.0) # x doesn't matter since we only use sections
    
    secref.sec.push()
    calc_path.end(1.0) # x doesn't matter (idem)
    
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList
    
    h.pop_section()
    h.pop_section()


    # Initialize variables
    path_L = 0.0
    path_L_elec = 0.0
    path_ri = 0.0

    # Compute path length
    path_secs = list(root_path)
    path_len = len(path_secs)
    
    for isec, psec in enumerate(path_secs):
        arrived_sec = (isec==path_len-1) # alternative: use sec.same()
        
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

def seg_path_ri(endseg, f, gleak_name):
    """ 
    Calculate axial path resistance from start of segment up to but 
    not including soma section (the topmost root section).

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
    rootsec.push()
    calc_path.begin(0.5)
    secref.sec.push()
    calc_path.end(0.5)
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList
    h.pop_section()
    h.pop_section()

    # Compute electrotonic path length
    pathri0 = 0. # axial path resistance from root sec to 0 end of this sec
    path_secs = list(root_path)
    assert(endseg.sec in path_secs)
    for i, psec in enumerate(path_secs): # walk sections
        for j_seg, seg in enumerate(psec): # walk segments
            if seg.sec.same(endseg.sec) and j_seg==j_endseg: # reached end segment
                assert same_seg(seg, endseg)
                pathri1 = pathri0 + seg.ri()
                return pathri0, pathri1
            pathri0 += seg.ri()

    raise Exception('End segment not reached')

def seg_path_L_elec(endseg, f, gleak_name):
    """ 
    Calculate electrotonic path length from after root section up to
    start of given segment.

    ALGORITHM
    - walk each segment from root section (child of top root) to the given
      segment and sum L/lambda for each segment

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
    rootsec.push()
    calc_path.begin(0.5)
    secref.sec.push()
    calc_path.end(0.5)
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList
    h.pop_section()
    h.pop_section()

    # Compute electrotonic path length
    path_L_elec = 0.0
    path_secs = list(root_path)
    assert(endseg.sec in path_secs)
    for i, psec in enumerate(path_secs):
        L_seg = psec.L/psec.nseg # segment length
        for j_seg, seg in enumerate(psec):
            # Check if we have reached our target segment
            if seg.sec.same(endseg.sec) and j_seg==j_endseg:
                assert same_seg(seg, endseg)
                return path_L_elec
            lamb_seg = seg_lambda(seg, gleak_name, f)
            L_elec = L_seg/lamb_seg
            path_L_elec += L_elec

    raise Exception('End segment not reached')

def seg_path_L(endseg, to_end):
    """
    Calculate path length from center of root section to given segment

    @param to_end   if True, return distance to end of segment, else
                    return distance to start of segment
    """
    secref = h.SectionRef(sec=endseg.sec)
    rootsec = subtreeroot(secref)
    j_endseg = seg_index(endseg)
    rootparent = rootsec.parentseg()
    if rootparent is None:
        return 0.0 # if we are soma/topmost root: path length is zero
    
    # Get path from root section to endseg
    calc_path = h.RangeVarPlot('v')
    rootsec.push()
    calc_path.begin(0.0) # x doesn't matter since we only use path sections
    endseg.sec.push()
    calc_path.end(endseg.x)
    root_path = h.SectionList() # SectionList structure to store path
    calc_path.list(root_path) # copy path sections to SectionList
    h.pop_section()
    h.pop_section()

    # Compute path length
    path_secs = list(root_path)
    path_L = 0.0
    for isec, psec in enumerate(path_secs):
        arrived = bool(psec.same(endseg.sec))
        seg_L = psec.L/psec.nseg 
        for j_seg, seg in enumerate(psec):
            if arrived and j_seg==j_endseg:
                # assert same_seg(seg, endseg)
                if to_end:
                    return path_L + seg_L
                else:
                    return path_L
            path_L += seg_L

def assign_electrotonic_length(rootref, allsecrefs, f, gleak_name, allseg=False):
    """ 
    Assign length constant (lambda) and electrotonic path length (L/lambda)
    to 0- and 1-end for each section in subtree of given section.

    @type   rootref     SectionRef
    @param  rootref     Section reference to current node

    @type   allsecrefs  list(SectionRef)
    @param  allsecrefs  references to all sections in the cell

    @post               all section references have the following attributes:
                        - 'f_lambda': frequency at which length constant is computed
                        - 'lambda_f': section's length constant at given frequency
                        - 'pathLelec0': electrotonic path length to 0-end
                        - 'pathLelec1': Electrotonic path length to 1-end
    """
    if rootref is None:
        return

    # Compute length constant
    gleak = sum([getattr(seg, gleak_name) for seg in rootref.sec])
    gleak /= rootref.sec.nseg # average gleak of segments in section
    lambda_f = sec_lambda(rootref.sec, gleak, f)
    rootref.lambda_f = lambda_f
    rootref.f_lambda = f

    # Compute electrotonic path length
    if allseg:
        rootref.pathL_elec = [0.0]*rootref.sec.nseg
        for i, seg in enumerate(rootref.sec):
            # visit every segment expcent 0/1 end zero-area segments
            pathL = seg_path_L_elec(seg, f, gleak_name)
            rootref.pathL_elec[i] = pathL
    pathLelec0, pathLelec1 = sec_path_L_elec(rootref, f, gleak_name)
    rootref.pathLelec0 = pathLelec0
    rootref.pathLelec1 = pathLelec1

    # Assign to children
    for childsec in rootref.sec.children():
        childref = getsecref(childsec, allsecrefs)
        assign_electrotonic_length(childref, allsecrefs, f, gleak_name, allseg=allseg)


################################################################################
# Clustering & Topology
################################################################################

class EqProps(object):
    """
    Equivalent properties of merged sections

    NOTE: this is the 'Bunch' recipe from the python cookbook
          al alternative would be `myobj = type('Bunch', (object,), {})()`
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

################################################################################
# Editing cells/dendritic trees
################################################################################

# def find_first_cut(noderef):
#   """ Utility function to find first segment in cluster """
#   for jseg, seg in enumerate(noderef.sec):
#       if noderef.cluster_labels[jseg] == cluster.label:
#           return seg
#
#   child_secs = noderef.sec.children()
#   if not any(child_secs):
#       return None
#
#   next_cut = None
#   for childsec in child_secs:
#       childref = getsecref(childsec, allsecrefs)
#       next_cut = find_cut(childref)
#       if next_cut is not None:
#           break
#   return next_cut

def find_prox_cuts(nodeseg, cuts, cluster, allsecrefs):
    """
    Find proximal cuts for given cluster.

    i.e. cluster boundaries given by pairs of segments where the first segment
    is proximal to the root section and in a different cluster, and the second
    segment is in the given cluster and distal to the first one.
    """
    noderef = getsecref(nodeseg.sec, allsecrefs)
    i_node = seg_index(nodeseg)
    node_label = noderef.cluster_labels[i_node]

    if node_label==cluster.label:
        # return nodeseg
        parseg = prev_seg(nodeseg)
        parref = getsecref(parseg.sec, allsecrefs)
        parlabel = parref.cluster_labels[seg_index(parseg)]
        if parlabel != cluster.label:
            cuts.add((parseg, nodeseg))
        else:
            return # don't descend any further
    else:
        chi_segs = next_segs(nodeseg)
        for chiseg in chi_segs:
            find_prox_cuts(chiseg, cuts, cluster, allsecrefs)

def find_dist_cuts(nodeseg, cuts, cluster, allsecrefs):
    """
    Find distal cuts for given cluster.

    @post   argument 'cuts' is filled with pairs (seg_a, seg_b) where
            seg_a is last segment in cluster along its path, and seg_b
            is first segment in next cluster
    """
    noderef = getsecref(nodeseg.sec, allsecrefs)
    i_node = seg_index(nodeseg)
    node_label = noderef.cluster_labels[i_node]

    if node_label != cluster.label:
        return
    else:
        chi_segs = next_segs(aseg)
        # if no child segs: this is a cut
        if not any(chi_segs):
            cuts.add((nodeseg,None))
            return
        # if child segs: check their labels
        for chiseg in chi_segs:
            chiref = getsecref(chiseg.sec, allsecrefs)
            chilabel = chiref.cluster_labels[seg_index(chiseg)]
            # if same label: continue searching
            if chilabel == cluster.label:
                find_dist_cuts(chiseg, cuts, cluster, allsecrefs)
            else: # if different label: this is a cut
                cuts.add((nodeseg,chiseg))

def substitute_cluster(rootref, cluster, clu_eqsec, allsecrefs, mechs_pars, 
                        delete_substituted=False):
    """
    Substitute equivalent section of cluster in dendritic tree.

    @param rootref      SectionRef to root section where algorithm
                        should start to descend tree.

    @param cluster      Cluster object containing information about
                        the cluster

    @param clu_eqsec    equivalent Section for cluster

    ALGORITHM
    - descend tree from root:
    - find first segment assigned to cluster
    - descend tree along all children and find last segments assigned to cluster
    - cut section of first segment (segments after first)
    - cut sections of last segments (segments before last)
    - fix the cut parts: re-set nseg, and active conductances
    - substitution: connect everything
    """

    # Mark all sections as uncut
    for secref in allsecrefs:
        if not hasattr(secref, 'is_cut'):
            setattr(secref, 'is_cut', False)

    # Find proximal cuts
    prox_cuts = set()
    find_prox_cuts(rootref.sec(0.0), prox_cuts, cluster, allsecrefs)
    if not any(prox_cuts):
        raise Exception("Could not find segment assigned to cluster in subtree of {}".format(rootref))
    elif len(prox_cuts) > 1:
        raise Exception("Inserting cluster at location where it has > 1 parent Section not supported")
    first_cut_seg = prox_cuts[0][1]
    
    # Find distal cuts
    dist_cuts = set() # pairs (seg_a, seg_b) where a is in same cluster, and b in next cluster
    find_dist_cuts(first_cut_seg, dist_cuts, cluster, allsecrefs)

    # Sever tree at proximal cuts
    cut_seg_index = seg_index(first_cut_seg)
    if cut_seg_index > 0:
        cut_sec = first_cut_seg.sec
        
        # Calculate new properties
        pre_cut_L = cut_sec.L/cut_sec.nseg * cut_seg_index
        cut_props = get_sec_properties(cut_sec, mechs_pars)
        
        # Set new properties
        cut_sec.nseg = cut_seg_index
        cut_sec.L = pre_cut_L
        for jseg, seg in enumerate(cut_sec):
            pseg = cut_props[jseg]
            for pname, pval in pseg.items():
                seg.__setattr__(pname, pval)

        # Mark as being cut
        cut_ref = getsecref(cut_sec, allsecrefs)
        cut_ref.is_cut = True

        # Set parent for equivalent section
        eqsec_parseg = cut_sec(1.0)
    else:
        eqsec_parseg = cut_sec.parentseg()

    # Sever tree at distal cuts
    eqsec_child_segs = []
    for dist_cut in dist_cuts:
        cut_seg = dist_cut[1] # Section will be cut before this segment
        cut_sec = cut_seg.sec
        cut_seg_index = seg_index(cut_seg)

        # If Section cut in middle: need to resize
        if cut_seg_index > 0:
            # Check if not already cut
            cut_ref = getsecref(cut_sec, allsecrefs)
            if cut_ref.is_cut:
                raise Exception('Section has already been cut before')

            # Calculate new properties
            new_nseg = cut_sec.nseg-(cut_seg_index+1)
            post_cut_L = cut_sec.L/cut_sec.nseg * new_nseg
            cut_props = get_sec_properties(cut_sec, mechs_pars)
            
            # Set new properties
            cut_sec.nseg = new_nseg
            cut_sec.L = post_cut_L
            for jseg, seg in enumerate(cut_sec):
                pseg = cut_props[cut_seg_index + jseg]
                for pname, pval in pseg.items():
                    seg.__setattr__(pname, pval)

        # Set child segment for equivalent section
        eqsec_child_segs.append(cut_sec(0.0))

    # Connect the equivalent section
    clu_eqsec.connect(eqsec_parseg, 0.0) # see help(sec.connect)
    for chiseg in eqsec_child_segs:
        chiseg.sec.connect(clu_eqsec, 1.0, chiseg.x) # this disconnects them

    # Disconnect the substituted parts
    for subbed_sec in eqsec_parseg.sec.children():
        chiref = getsecref(subbed_sec, allsecrefs)
        if chiref.cluster_labels[0] == cluster.label:
            # NOTE: only need to disconnect proximal ends of the subtree: the distal ends have already been disconnected by connecting them to the equivalent section for the cluster
            h.disconnect(sec=subbed_sec)

            # Delete disconnected subtree if requested
            if delete_substituted:
                child_tree = h.SectionList()
                subbed_sec.push()
                child_tree.subtree()
                h.pop_section()
                for sec in child_tree: # iterates CAS
                    h.delete_section()

def sub_equivalent_Y_sec(eqsec, parent_seg, bound_segs, allsecrefs, mechs_pars, 
                            delete_substituted):
    """
    Substitute equivalent section of a collapsed Y into the tree.

    @param eqsec    equivalent Section for cluster

    @param parent_seg   parent segment (where equivalent section will be attached)

    @param bound_segs   list of distal boundary segments of the cluster
                        of segments that the equivalent section is substituting for

    """

    # Mark all sections as uncut
    for secref in allsecrefs:
        if not hasattr(secref, 'is_cut'):
            setattr(secref, 'is_cut', False)

    # Parent of equivalent section
    parent_seg = parent_seg

    # Sever tree at distal cuts
    eqsec_child_segs = []
    for bound_seg in bound_segs: # last absorbed/collapsed segment
        cut_segs = next_segs(bound_seg) # first segment after cut (not collapsed)
        for i in range(len(cut_segs)):
            cut_seg = cut_segs[i] # cut_seg may be destroyed by resizing so don't make loop variable
            cut_sec = cut_seg.sec
            cut_seg_index = seg_index(cut_seg)

            # If Section cut in middle: need to resize
            if cut_seg_index > 0:
                # Check if not already cut
                cut_ref = getsecref(cut_sec, allsecrefs)
                if cut_ref.is_cut:
                    raise Exception('Section has already been cut before')

                # Calculate new properties
                new_nseg = cut_sec.nseg-(cut_seg_index)
                post_cut_L = cut_sec.L/cut_sec.nseg * new_nseg
                cut_props = get_sec_properties(cut_sec, mechs_pars)
                
                # Set new properties
                cut_sec.nseg = new_nseg
                cut_sec.L = post_cut_L
                for jseg, seg in enumerate(cut_sec):
                    pseg = cut_props[cut_seg_index + jseg]
                    for pname, pval in pseg.items():
                        seg.__setattr__(pname, pval)
                logger.debug("Cut section {0} and re-assigned segment properties.".format(cut_sec.name()))

            # Set child segment for equivalent section
            eqsec_child_segs.append(cut_sec(0.0))

    # Connect the equivalent section
    eqsec.connect(parent_seg, 0.0) # see help(sec.connect)
    for chi_seg in eqsec_child_segs:
        chi_seg.sec.connect(eqsec, 1.0, chi_seg.x) # this disconnects them

    # Disconnect the substituted parts
    for subbed_sec in parent_seg.sec.children():
        # Skip the equivalent/substituting section
        if subbed_sec.same(eqsec):
            continue

        # Only need to disconnect proximal ends of the subtree: the distal ends have already been disconnected by connecting them to the equivalent section for the cluster
        logger.debug("Disconnecting substituted section '{}'...".format(subbed_sec.name()))
        h.disconnect(sec=subbed_sec)

        # Delete disconnected subtree if requested
        child_tree = h.SectionList()
        subbed_sec.push()
        child_tree.subtree()
        h.pop_section()
        for sec in child_tree: # iterates CAS
            subbed_ref = getsecref(sec, allsecrefs)
            subbed_ref.is_substituted = True
            if delete_substituted:
                logger.debug("Deleting substituted section '{}'...".format(sec.name()))
                subbed_ref.is_deleted = True # can also be tested with NEURON ref.exists()
                h.delete_section()
                

def split_section(src_sec, mechs_pars, delete_src=False):
    """
    Split section by deleting it and adding two sections in series

    @param mechs_pars   dictionary of mechanism name -> [parameter names]
                        that need to be copied
    """
    # Create two halves
    halfA_name = src_sec.name() + "_A"
    h("create %s" % halfA_name)
    secA = getattr(h, halfA_name)

    halfB_name = src_sec.name() + "_B"
    h("create %s" % halfB_name)
    secB = getattr(h, halfB_name)

    # Copy properties
    for tar_sec in [secA, secB]:
        # Copy passive properties
        tar_sec.L = src_sec.L / 2.
        tar_sec.Ra = src_sec.Ra
        if src_sec.nseg % 2 == 0:
            tar_sec.nseg = src_sec.nseg / 2
        else:
            logger.warning("Splitting section with uneven number of segments")
            tar_sec.nseg = int(src_sec.nseg) / 2 + 1 # don't lose accuracy

        # Copy mechanisms
        for mechname in mechs_pars.keys():
            if hasattr(src_sec(0.5), mechname):
                tar_sec.insert(mechname)

        # Copy range variables
        for tar_seg in tar_sec:
            src_seg = src_sec(tar_seg.x/2.) # segment that it maps to
            tar_seg.diam = src_seg.diam
            tar_seg.cm = src_seg.cm
            for mech in mechs_pars.keys():
                for par in mechs_pars[mech]:
                    prop = par+'_'+mech
                    setattr(tar_seg, prop, getattr(src_seg, prop))

        # Copy ion styles
        copy_ion_styles(src_sec, tar_sec)

    # Connect A to B
    secB.connect(secA, 1, 0)

    # Connect A to parent of src_sec
    secA.connect(src_sec.parentseg().sec, src_sec.parentseg().x, 0)

    # Connect children of src_sec to B
    for childsec in src_sec.children():
        xloc = childsec.parentseg().x
        # NOTE: connecting section disconnects it from previous parent
        if xloc >= 0.5:
            childsec.connect(secB, 2*(xloc-0.5), 0)
        else:
            childsec.connect(secA, 2*xloc, 0)

    # Disconnect and delete src_sec
    h.disconnect(sec=src_sec)
    if delete_src:
        h.delete_section(sec=src_sec)

def join_unbranched_sections(nodesec, max_joins=1e9):
    """
    Join unbranched sections into a single section.

    This operation preserves the number of segments and distribution of 
    diameter and electrical properties.

    @pre    the method assumes Ra is uniform in all the sections
    """
    pass


def get_sec_properties(src_sec, mechs_pars):
    """
    Get RANGE properties for each segment in Section.

    @return     list(dict()) of length nseg containing property names and 
                values for each segment in Section
    """
    props = []
    for seg in src_sec:
        pseg = {}
        pseg['diam'] = seg.diam
        pseg['cm'] = seg.cm
        for mech in mechs_pars.keys():
            for par in mechs_pars[mech]:
                prop = par+'_'+mech
                pseg[prop] = getattr(seg, prop)
        props.append(pseg)
    return props


def get_range_props(secref, prop_names):
    """
    Get Section's requested RANGE properties for each segment.
    """

    seg_prop_dicts = [dict() for i in range(secref.sec.nseg)]

    for j_seg, seg in enumerate(secref.sec):
        for prop in prop_names:
            seg_prop_dicts[j_seg][prop] = getattr(seg, prop)

    return seg_prop_dicts


def set_range_props(secref, seg_prop_dicts):
    """
    Set Section's RANGE properties stored in given dictionaries.
    """
    if secref.sec.nseg != len(seg_prop_dicts):
        raise ValueError("Must provide one dictionary for each segment.")
    
    for j_seg, seg in enumerate(secref.sec):
        propdict = seg_prop_dicts[j_seg]
        for pname, pval in propdict.items():
            setattr(seg, pname, pval)


def get_sec_props_obj(secref, mechs_pars, seg_assigned, sec_assigned):
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


    @return     object EqProps with the desired properties stored as attributes
    """

    # Store section properties (non-Range)
    sec_props = EqProps(L=secref.sec.L, Ra=secref.sec.Ra, nseg=secref.sec.nseg)
    for prop in sec_assigned:
        setattr(sec_props, prop, getattr(secref, prop))

    # Initialize segment RANGE properties
    sec_props.seg = [dict() for i in range(secref.sec.nseg)]
    bprops = [par+'_'+mech for mech,pars in mechs_pars.items() for par in pars] # NEURON properties
    
    # Store segment RANGE properties
    for j_seg, seg in enumerate(secref.sec):
        
        # Store built-in properties
        for prop in bprops:
            sec_props.seg[j_seg][prop] = getattr(seg, prop)
        
        # Store self-assigned properties (stored on SectionRef)
        for prop in seg_assigned:
            sec_props.seg[j_seg][prop] = getattr(secref, prop)[j_seg]
    
    return sec_props


def store_seg_props(secref, mechs_pars, attr_name='or_seg_props', assigned_props=None):
    """
    Store each segment's properties (RANGE variables) in a dictionary
    on the given SectionRef.

    @param mechs_pars       dict mechanism_name -> [paramerter_names] with segment properties
                            to store

    @param assigned_props   list of assigned properties (assigned attributed of SectionRef)
                            that should be saved

    @param attr_name        the attribute name to use for saving the properties, in format:
                            list(dict(param : value))

    @effect     SectionRef will have an attribute 'attr_name' that is
                a list with a dict for each segment containing all its properties.
    """

    # Initialize segment RANGE properties
    secref.__setattr__(attr_name, [dict() for i in range(secref.sec.nseg)])
    nrn_props = [par+'_'+mech for mech,pars in mechs_pars.items() for par in pars] # NEURON properties
    ref_props = [] if assigned_props is None else assigned_props

    # Store segment RANGE properties
    for j_seg, seg in enumerate(secref.sec):
        # Store built-in properties
        for prop in nrn_props:
            secref.__getattribute__(attr_name)[j_seg][prop] = getattr(seg, prop)
        # Store self-assigned properties (stored on SectionRef)
        for prop in ref_props:
            secref.__getattribute__(attr_name)[j_seg][prop] = getattr(secref, prop)[j_seg]




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


def copy_ion_styles(src_sec, tar_sec, ions=None):
    """
    Copy ion styles from source to target section

    NOTE:

    oldstyle = ion_style("name_ion")

    oldstyle = int:
        int( 1*c_style + 4*cinit + 8*e_style + 32*einit + 64*eadvance )
        c_style:    0, 1, 2, 3  (2 bits)
        e_style:    0, 1, 2, 3  (2 bits)
        einit:      0, 1        (1 bits)
        eadvance:   0, 1        (1 bits)
        cinit:      0, 1        (1 bits)

    ion_style("name_ion", c_style, e_style, einit, eadvance, cinit)

    """
    if ions is None:
        ions = ['na', 'k', 'ca']

    # Get ion style for each ion species
    src_sec.push()
    styles = dict(((ion, h.ion_style(ion+'_ion')) for ion in ions))
    
    # Copy to target Section
    set_ion_styles(tar_sec, **styles)

    h.pop_section()


def get_ion_styles(src_sec, ions=None):
    """
    Get ion styles as integer for each ion.
    """
    if ions is None:
        ions = ['na', 'k', 'ca']

    # Get ion style for each ion species
    src_sec.push()
    styles = dict(((ion, h.ion_style(ion+'_ion')) for ion in ions))
    h.pop_section()

    return styles


def set_ion_styles(tar_sec, **kwargs):
    """
    Set ion styles from integer containing bit flags.

    @param  kwargs      keyword arguments ion_name: style_int
    """
    # Copy to target Section
    tar_sec.push()
    for ion, style in kwargs.items():

        # Decompose int into bit flags
        c_style = int(style) & (1+2)
        cinit = (int(style) & 4) >> 2
        e_style = (int(style) & (8+16)) >> 3
        einit = (int(style) & 32) >> 5
        eadvance = (int(style) & 64) >> 6

        # Copy to new section
        h.ion_style(ion+'_ion', c_style, e_style, einit, eadvance, cinit)
    
    h.pop_section()


def duplicate_subtree(rootsec, mechs_pars, tree_copy):
    """ Duplicate tree of given section
    @param rootsec      root section of the subtree
    @param mechs_pars   dictionary mechanism_name -> parameter_name
    @param tree_copy    out argument: list to be filled
    """
    # Copy current root node
    copyname = 'copyof_' + rootsec.name()
    i = 0
    while h.issection(copyname):
        if i > 1000:
            raise Exception('Too many copies of this section!')
        i += 1
        copyname = ('copy%iof' % i) + rootsec.name()
    h("create %s" % copyname)
    root_copy = getattr(h, copyname)
    copy_sec_properties(rootsec, root_copy, mechs_pars)
    tree_copy.append(root_copy)

    # Copy children
    for childsec in rootsec.children():
        child_copy = duplicate_subtree(childsec, mechs_pars, tree_copy)
        child_copy.connect(root_copy, childsec.parentseg().x, 0)

    return root_copy


if __name__ == '__main__':
    print("Main of redutils.py: TODO: execute tests")