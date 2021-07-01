"""
Merging rules for Marasco (2012) reduction

@author     Lucas Koelman
"""

# Python modules
import math
PI = math.pi

# Our modules
from .redutils import EqProps, getsecref, seg_index


################################################################################
# Section-based merging
################################################################################


def merge_sec_parallel(childrefs, allsecrefs, glist):
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
        L_child, diam_child, Ra_child, ri_child, cmtot_child, gtot_child = merge_sec_sequential(childref, allsecrefs)
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


def merge_sec_sequential(rootref, allsecrefs, glist):
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
    L_br, diam_br, Ra_br, rin_br, rax_br, cmtot_br, gtot_br = merge_sec_parallel(childrefs, allsecrefs)

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


################################################################################
# Segment-based collapsing
################################################################################


def merge_seg_subtree(rootseg, allsecrefs, eligfunc=None, modfunc=None, bound_segs=None):
    """ 
    Recursively merge within-cluster connected sections in subtree
    of the given segment using equations <br> and <seq> in Marasco (2012).

    @pre        in the tree, all connections between sections must be made
                according to the convention that the 0-end of the child
                connects to the 1-end of the parent

    @param eligfunc     function(seg,jseg,secref) -> bool indicating whether
                        the segment is eligible to be absorbed (seg is the segment,
                        jseg is its index, secref is a SectionRef to its section)

    ALGORITHM
    
    - for each child segment: recursively call `equivalent_properties = merge_sec_sequential(child)`
    
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

    child_segs = []
    root_is_bound = False # whether current root is a boundary segment

    # Find absorbable child segments
    if i_rootseg == rootsec.nseg-1:
        # IF END SEGMENT: get child segments that are in same cluster
        for secref in child_refs:

            # Get first segment (0-end) of child
            childseg = next(seg for seg in secref.sec)
            
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
        cp, _ = merge_seg_subtree(child, allsecrefs, eligfunc, modfunc, bound_segs)
        ch_area = PI*cp.diam_eq*cp.L_eq

        # Update equivalent properties of all child branches in parallel
        area_br += ch_area
        L_br += ch_area*cp.L_eq     # Length: Eq. L_br
        diam_br += cp.diam_eq**2    # Diameter: Eq. rho_br
        Ra_br += cp.Ra_eq           # Axial resistivity: average Ra
        Ri_br *= cp.Ri_eq           # Axial resistance: Eq. r_a,br
        Ri_br_sum += cp.Ri_eq       # Sum of branch axial resistances

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
    L_seq = L_root + L_br           # Eq. L_seq
    Ra_seq = (Ra_root + Ra_br)/2.   # for rho_seq calculation
    Ri_seq = Ri_root + Ri_br        # Eq. r_a,seq
    Ri_seq_eqgeom = Ri_root + Ri_br_eqgeom                      # for diam_seq calculation
    diam_seq = math.sqrt(Ra_seq*L_seq*4./PI/Ri_seq_eqgeom/100.) # Eq. rho_seq
    
    eq_props_seq = EqProps(L_eq=L_seq, diam_eq=diam_seq, Ra_eq=Ra_seq, Ri_eq=Ri_seq)
    return eq_props_seq, eq_props_br
