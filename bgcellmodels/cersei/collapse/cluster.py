"""
Functions for clustering sections in dendritic trees.
"""

import re

from . import redutils as redtools
from .redutils import getsecref

from bgcellmodels.common import logutils
logger = logutils.getLogger(__name__)


class Cluster(object):
    """
    A cluster representing merged sections
    """

    def __init__(self, label):
        self.label = label

    def __repr__(self):
        return '{0} ({1})'.format(self.label, super(Cluster, self).__repr__())

    # def __repr__(self):
    #   desc = super(Cluster, self).__repr__()
    #   printable = ['label', 'eqL', 'eqdiam', 'eqSurf', 'totsurf', 'surfSum']
    #   for ppty in printable:
    #       if hasattr(self, ppty):
    #           desc += '\n\t|- {0}: {1}'.format(ppty, getattr(self, ppty))
    #   return desc


def clusterroot(secref, allsecrefs):
    """
    Find the highest parent/ancestor of given section that is still
    in the same cluster
    """
    if not secref.has_parent():
        return secref
    else:
        parref = getsecref(secref.parent, allsecrefs)
        if parref is None: # case where parent is not in allsecrefs
            return secref
        elif parref.cluster_label != secref.cluster_label:
            return secref
        else:
            return clusterroot(parref)


def assign_topology_attrs(
        cur_ref, 
        all_refs, 
        par_ref=None, 
        root_order=None
    ):
    """ 
    Assign numbers that reflect topological location of sections

    @pre    a precondition/assumption for the procedure is that NEURON's
            standard assumption that a Section = a continuous length of
            UNBRANCHED cable is TRUE (i.e. branch points can only occur 
            at 0/1 ends and never in middle of Section).


    @post   Each SectionRef has the following attributes:
            
            - ref.order = number of sections between root node and section
            
            - ref.level = number of branch points between root node and section

                * level = 1 for initial section passed to this function

                * level = n + 1 after each branch point
            
            - ref.strahlernumber = 
                
                for leaf sections: 
                    * n = 1, 
                
                for parent sections:
                    * n = max(n_i) for children i if n_i are different, 
                    * n = (n_i + 1) if child's n_i equal

                (i.e. assign level starting from leaf nodes, and keep minimum)

            - ref.end_branchpoint = True if the Section has a branchpoint at its end
              (1-end segment is branchpoint)


    @param  cur_ref     SectionRef to current node


    @param  all_refs    list(SectionRef) to all sections in the cell


    @param  par_order   order of parent sections (distance in #sections from soma)
    
    """
    if cur_ref is None:
        logger.debug("Node not found in SectionRef list")
        return

    # assign order and level
    if par_ref is None: # start from beginning
        cur_ref.order = 1
        cur_ref.level = 1
    
    else: # determine from parent
        
        cur_ref.order = par_ref.order + 1
        
        if par_ref.end_branchpoint:
            # increase level if parent section is branch point
            cur_ref.level = par_ref.level + 1
        else:
            cur_ref.level = par_ref.level

    # User override of inital order
    if root_order is not None:
        cur_ref.order = root_order
        cur_ref.level = root_order

    childsecs = cur_ref.sec.children()
    childrefs = [getsecref(sec, all_refs) for sec in childsecs]

    # Check if current cylinder ends in branch point (fork)
    if not any(childsecs):

        # LEAF NODE: strahler number = 1, stop descending
        cur_ref.strahlernumber = 1
        cur_ref.end_branchpoint = False
        return
    
    else:
        cur_ref.end_branchpoint = (len(childsecs) > 1)

    # NON LEAF NODE: descend tree and calculate strahler
    for childref in childrefs:
        assign_topology_attrs(childref, all_refs, par_ref=cur_ref)

    # After return from tree ascent: query Strahler numbers
    strahler_numbers = [ref.strahlernumber for ref in childrefs]
    max_strahler = max(strahler_numbers)
    
    # Only increase if all the same, else take maximum
    if all([ref.strahlernumber==max_strahler for ref in childrefs]):
        cur_ref.strahlernumber = max_strahler + 1
    
    else:
        cur_ref.strahlernumber = max_strahler


def clusterize_custom(noderef, allsecrefs, clusterlist, labelsuffix, clusterfun, cluster_args):
    """
    Cluster dendritic tree based on custom criterion

    @type   noderef     SectionRef
    @param  noderef     Section reference to current node

    @type   allsecrefs  list(SectionRef)
    @param  allsecrefs  references to all sections in the cell

    @type   clusterlist list(Cluster)
    @param  clusterlist list of existing Clusters

    @type   clusterfun  function(SectionRef) -> string
    @param  clusterfun  function mapping a SectionRef to cluster label
    """
    if noderef is None: return

    # Cluster section based on custom criterion
    cluster_args['labelsuffix'] = labelsuffix
    labels = clusterfun(noderef, **cluster_args)

    # Add new Cluster object if necessary
    for label in labels:
        if not any(clu.label==label for clu in clusterlist):
            newclu = Cluster(label)
            clusterlist.append(newclu)

    # Cluster children (iteratively)
    for childsec in noderef.sec.children():
        childref = getsecref(childsec, allsecrefs)
        clusterize_custom(childref, allsecrefs, clusterlist, labelsuffix, clusterfun, cluster_args)

def label_from_custom_regions(noderef, label_seg=True, labelsuffix='', marker_mech=None):
    """
    Assign cluster label based on functional regions identified from
    channel distribution (max conductance for each channel)

    @param marker_mech  a pair (mechanism_name, mechanism_ppty) indicating
                        a mechanism that should be inserted to flag segments
                        with their cluster label. The property will be used
                        to set an integer determined by the label.
    """
    if not hasattr(noderef, 'cluster_labels'):
        noderef.cluster_labels = [None] * noderef.sec.nseg

    if marker_mech is not None:
        noderef.sec.insert(marker_mech[0])

    for iseg, seg in enumerate(noderef.sec):
        # Compute path length
        path_L = redtools.seg_path_L(seg, endpoint='segment_end')

        # Determine label from threshold (see gbar plot)
        if path_L > 220.:
            label = 'spiny' + labelsuffix
            flag = 3
        elif path_L > 90.:
            label = 'smooth' + labelsuffix
            flag = 2
        else:
            label = 'trunk' + labelsuffix
            flag = 1
        noderef.cluster_labels[iseg] = label

        # Flag the segment based on its label
        if marker_mech is not None:
            propname = marker_mech[1] + '_' + marker_mech[0]
            seg.__setattr__(propname, flag)

    return noderef.cluster_labels


def label_from_strahler(noderef, thresholds=None, labelsuffix=''):
    """
    Cluster a tree based on Strahler numbers alone

    @effect         assign label 'trunk'/'smooth'/spiny' to sections
                    based on their Strahler number
    """
    if noderef is None: return

    # Clustering thresholds
    if thresholds is None:
        thresholds = (3,5)

    # Cluster label based on strahler's number
    if noderef.strahlernumber <= thresholds[0]: # default: <= 3
        label = 'spiny' + labelsuffix
    elif noderef.strahlernumber <= thresholds[1]: # default: <= 5
        label = 'smooth' + labelsuffix
    else:
        label = 'trunk' + labelsuffix
    noderef.cluster_label = label
    return [label]


def label_sec_electrotonic(noderef, thresholds, labelsuffix=''):
    """ 
    Cluster dendritic tree based on length divided by length constant
    (i.e. length in terms of its electrotonic length) measured from soma section.

    @type   thresholds  tuple(float)
    @param  thresholds  one or two thresholds on electrotonic path length 
                        to distinguish between trunk/smooth/spiny sections

    @pre                all section references must have their length constant
                        set using assign_electrotonic_length
    """
    if noderef is None:
        return

    # Cluster based on electronic path length at midpoint
    L_mid = noderef.pathLelec0 + (noderef.pathLelec1 - noderef.pathLelec0)/2.
    if L_mid <= thresholds[0]:
        label = 'trunk' + labelsuffix
    elif len(thresholds) > 1 and L_mid <= thresholds[1]:
        label = 'smooth' + labelsuffix
    else:
        label = 'spiny' + labelsuffix
    noderef.cluster_label = label

    logger.debug("Section with index {} assigned to cluster {}".format(
                    noderef.table_index, noderef.cluster_label))
    return [label]

    
def label_seg_electrotonic(noderef, thresholds, labelsuffix=''):
    """ 
    Cluster segments in dendritic tree based on length divided by 
    length constant (i.e. length in terms of its electrotonic length) 
    measured from soma section.

    @type   thresholds  tuple(float)
    @param  thresholds  one or two thresholds on electrotonic path length 
                        to distinguish between trunk/smooth/spiny sections

    @pre                all section references must have their length constant
                        set using assign_electrotonic_length
    """
    if noderef is None:
        return

    # Cluster each segment
    noderef.cluster_labels = ['unassigned'] * noderef.sec.nseg
    for i in range(noderef.sec.nseg):
        pathL = noderef.pathL_elec[i] # electrotonic path length of segment i
        if pathL <= thresholds[0]:
            noderef.cluster_labels[i] = 'trunk' + labelsuffix
        elif len(thresholds) > 1 and pathL <= thresholds[1]:
            noderef.cluster_labels[i] = 'smooth' + labelsuffix
        else:
            noderef.cluster_labels[i] = 'spiny' + labelsuffix
    return noderef.cluster_labels


def cluster_topology(rootref, allsecrefs, relations):
    """
    Determine cluster topology from list of clustered section references

    @param relations    a container for relations of the form (parentlabel, childlabel)
                        if this is a set, the entries are guaranteed to be unique
    """
    # Depth-first recursion of tree
    for childsec in rootref.sec.children():
        childref = getsecref(childsec, allsecrefs)
        # Add parent-child relationship
        if childref.cluster_label != rootref.cluster_label:
            relations.add((rootref.cluster_label, childref.cluster_label))
        # determine topology of subtree
        cluster_topology(childref, allsecrefs, relations)


def assign_topology(
        clusters, 
        radial_prefixes, 
        prefix_pattern=r'^[a-zA-Z0-9]+', 
        suffix_pattern=r'_[0-9]+$'
    ):
    """ 
    Assign attributes 'parent_label' and 'parent_pos' for each Cluster in clusters
    based on the given prefixes, which should be ordered from proximal to distal.
    """
    # Assign defaults
    for clu in clusters:
        if not hasattr(clu, 'order'):
            clu.order = next(i for i,prefix in enumerate(radial_prefixes) if clu.label.startswith(prefix))

    # Assign based on topology
    for clu in clusters:
        # Get current cluster prefix
        prefix_find = re.compile(prefix_pattern) # same as ^[^\W_]+ : matches anything before first underscore
        match = re.search(prefix_find, clu.label) 
        clu_prefix = match.group()

        # Defaults
        clu.parent_label = clu.label # if no parent found it is root cluster
        clu.parent_pos = 1.0

        # look for next prefix in anti-radial direction that is available
        parent_index = radial_prefixes.index(clu_prefix)

        # Root cluster
        if parent_index == 0:
            clu.parent_label = clu.label
            clu.order = 0

        # utility function to find parent cluster
        def parent_generator(child_clu, parent_index):
            parent_prefix = radial_prefixes[parent_index]
            for clu in clusters:
                if clu.label.startswith(parent_prefix):
                    match_suffix = re.search(suffix_pattern, clu.label)
                    suffix_ok = (match_suffix is None) or child_clu.label.endswith(match_suffix.group())
                    if suffix_ok:
                        yield clu

        # Non-root clusters
        while parent_index > 0:
            parent_index -= 1
            parent_gen = parent_generator(clu, parent_index)
            parent_clu = next(parent_gen, None)
            if parent_clu is not None:
                logger.debug("Found parent <%s> for cluster <%s>" % (parent_clu.label, clu.label))
                clu.parent_label = parent_clu.label
                clu.order = parent_clu.order + 1
                break

        # Sanity check
        assert(parent_index==0 or clu.parent_label!=clu.label)
