"""
Functions for editing neuronal trees in NEURON.

@author		Lucas Koelman
"""

from neuron import h

from bgcellmodels.common import logutils
from bgcellmodels.common.nrnutil import getsecref, seg_index
from bgcellmodels.common.treeutils import next_segs

from . import redutils

logger = logutils.getLogger('redops')


def find_collapsable(
        allsecrefs, 
        i_pass, 
        Y_criterion='highest_level', 
        zips_per_pass=1e9,
        **kwargs):
    """
    Find branch points with child branches that can be collapsed.

    @return     list of SectionRef to the branchpoint-Sections with
                collapsable children
    """

    if Y_criterion=='max_electrotonic':
        
        # Find Y with longest electrotonic length to first Y among its children.
        # This corresponds to  Y with child section that has longest L/lambda.
        # (Collapsing this Y will eliminate most compartments.)

        # Get minimum distance to next branch point (determines collapsable length)
        for secref in allsecrefs:
            
            child_secs = secref.sec.children()
            secref.collapsable_L_elec = 0.0
            
            if any(child_secs):
                min_child_L = min(sec.L for sec in child_secs)
            
            for chi_sec in child_secs:
                # Get segment that is that distance away from child
                furthest_seg = chi_sec(min_child_L/chi_sec.L)
                furthest_L_elec = redutils.seg_path_L_elec(
                                    furthest_seg, kwargs['f_lambda'], kwargs['gleak_name'])
                secref.collapsable_L_elec += (secref.pathLelec1 - furthest_L_elec)

        # Get maximum collapsable length
        max_L_collapsable = max(ref.collapsable_L_elec for ref in allsecrefs)

        # Find all branch points that have +/- same collapsable length
        low, high = 0.95*max_L_collapsable, 1.05*max_L_collapsable
        candidate_Y_secs = [ref for ref in allsecrefs if (any(ref.sec.children())) and (
                                                            low<ref.collapsable_L_elec<high)]
        
        # Off these, pick the one with highest level, and select all branch points at this level
        target_level = max(ref.level for ref in candidate_Y_secs)
        target_Y_secs = [ref for ref in candidate_Y_secs if ref.level==target_level]

        # Report findings
        logger.debug("The maximal collapsable electrotonic length is %f", max_L_collapsable)
        logger.debug("Found %i sections with collapsable length within 5%% of this value", 
                        len(candidate_Y_secs))
        logger.debug(("The highest level of a parent node with this value is %i, "
                        "and the number of nodes at this level with the same value is %i"), 
                        target_level, len(target_Y_secs))


    elif Y_criterion in ['max_level', 'highest_level']:
        # Simplify find all branch points at highest level
        max_level = max(ref.level for ref in allsecrefs)

        target_Y_secs = [ref for ref in allsecrefs if 
                            (ref.level == max_level-1)
                            and (len(ref.sec.children()) >= 2)
                            and (ref.max_passes >= i_pass+1)]

    elif Y_criterion == 'exact_level':

        target_level = kwargs['level']
        target_Y_secs = redutils.find_roots_at_level(target_level, allsecrefs)

    elif Y_criterion == 'outer_branchpoints':
        # All branch points that have no child branch points
        root_sec = allsecrefs[0].root; h.pop_section() # pushes CAS
        root_ref = getsecref(root_sec, allsecrefs)

        target_Y_secs = redutils.find_outer_branchpoints(root_ref, allsecrefs)
        

    else:
        raise Exception("Unknow Y-section selection criterion '{}'".format(Y_criterion))


    # Prune branchpoints: can't collapse more branchpoints than given maximum
    n_Y = len(target_Y_secs)
    target_Y_secs = [ref for i,ref in enumerate(target_Y_secs) if i<zips_per_pass]

    logger.debug("Found {0} Y-sections that meet selection criterion. Keeping {1}/{2}".format(n_Y, min(n_Y, zips_per_pass), n_Y))
    

    # Return branch points identified for collapsing
    return target_Y_secs


def sub_equivalent_Y_sec(
        eqsec,
        parent_seg,
        dist_bound_segs,
        allsecrefs,
        mechs_pars, 
        delete_substituted
    ):
    """
    Substitute equivalent section of a collapsed Y into the tree.

    @param eqsec    equivalent Section for cluster

    @param parent_seg   parent segment (where equivalent section will be attached)

    @param dist_bound_segs   list of distal boundary segments,
                        i.e. last absorbed/collapsed segment

    """

    # Mark all sections as uncut
    for secref in allsecrefs:
        if not hasattr(secref, 'is_cut'):
            setattr(secref, 'is_cut', False)

    # Sever tree at distal cuts
    behind_cut_segs = []
    
    for dist_seg in dist_bound_segs: # last absorbed/collapsed segment

        # first segment after cut (not collapsed)
        cut_segs = next_segs(dist_seg) 
        
        for i in xrange(len(cut_segs)):
            cut_seg = cut_segs[i] # cut_seg may be destroyed by resizing so don't make loop variable
            cut_sec = cut_seg.sec
            cut_seg_index = seg_index(cut_seg)

            # If Section cut in middle: need to resize
            if cut_seg_index > 0:
                # Check if not already cut
                cut_ref = getsecref(cut_sec, allsecrefs)
                if cut_ref.is_cut:
                    raise Exception('Section has already been cut before')

                # Resize section to portion behind cut
                new_nseg = cut_sec.nseg-(cut_seg_index)
                post_cut_L = cut_sec.L/cut_sec.nseg * new_nseg
                cut_props = redutils.get_sec_range_props(cut_sec, mechs_pars)
                
                # Reassign properties behind cut
                cut_sec.nseg = new_nseg
                cut_sec.L = post_cut_L
                for jseg, seg in enumerate(cut_sec):
                    pseg = cut_props[cut_seg_index + jseg]
                    for pname, pval in pseg.iteritems():
                        seg.__setattr__(pname, pval)
                
                logger.debug("Cut section {0} and re-assigned segment properties.".format(cut_sec.name()))

            # Set child segment for equivalent section
            behind_cut_segs.append(cut_sec(0.0))

    # Connect the equivalent section
    eqsec.connect(parent_seg, 0.0) # see help(sec.connect)
    for post_seg in behind_cut_segs:
        post_seg.sec.connect(eqsec(1.0), 0.0) # this disconnects them

    # Disconnect sections absorbed into equivalent section that were not cut)
    for subbed_sec in parent_seg.sec.children():
        # NOTE: remaining children directly attached to parent are the equivalent
        #       section and any sections that were not cut
        if subbed_sec.same(eqsec):
            continue

        # Only need to disconnect proximal ends of the subtree: the distal ends have already been disconnected by connecting them to the equivalent section for the cluster
        logger.debug("Disconnecting substituted section '{}'...".format(subbed_sec.name()))
        h.disconnect(sec=subbed_sec)

        # Delete disconnected subtree if requested
        child_tree = h.SectionList()
        child_tree.subtree(sec=subbed_sec)
        for sec in child_tree: # iterates CAS
            subbed_ref = getsecref(sec, allsecrefs)
            subbed_ref.is_substituted = True
            if delete_substituted:
                logger.debug("Deleting substituted section '{}'...".format(sec.name()))
                subbed_ref.is_deleted = True # can also be tested with NEURON ref.exists()
                h.delete_section()


def disconnect_subtree(parent_ref, candidate_refs, should_disconnect, delete=False):
    """
    Disconnect substituted sections and delete them if requested.

    @param  is_substituted : callable(SectionRef) -> bool
            Function to check if Section is substituted
    """
    
    # Delete disconnected subtree if requested
    subtree_seclist = h.SectionList()
    subtree_seclist.subtree(sec=parent_ref.sec)
    
    for child_sec in subtree_seclist: # iterates CAS
        
        if child_sec.same(parent_ref.sec):
            continue
        child_ref = getsecref(child_sec, candidate_refs)
        if child_ref is None: # Section not in disconnectable candidates
            continue
        
        if should_disconnect(child_ref):
            h.disconnect(sec=child_sec)
            if delete:
                logger.debug("Deleting substituted section '{}'.".format(child_sec.name()))
                child_ref.is_deleted = True # can also be tested with NEURON ref.exists()
                h.delete_section()

    # Alternative approach: don't use SectionList.subtree() but do recursive
    # call where children are processed first and then the current




def merge_linear_sections(nodesec, max_joins=1e9):
    """
    Join unbranched sections into a single section.

    This operation preserves the number of segments and distribution of 
    diameter and electrical properties.

    @pre    the method assumes Ra is uniform in all the sections
    """
    # TODO: implement merge_linear_sections function, with considerations:
    #   - density mechanisms * area they occupy preserved
    #   - Ra must be same
    pass