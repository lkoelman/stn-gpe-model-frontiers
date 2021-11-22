"""
Morphology reduction using method described in Marasco & Migliore (2012)

@author Lucas Koelman
@date   5-12-2016
"""

# Python modules
import re
import math
import logging

# NEURON modules
import neuron
h = neuron.h
h.load_file("stdlib.hoc") # Load the standard library

# Own modules
import bgcellmodels.common.electrotonic as electro
from bgcellmodels.common import logutils
from bgcellmodels.common.nrnutil import ExtSecRef, seg_index
from bgcellmodels.common.treeutils import subtree_topology

from . import (
    redutils,
    tree_edit as treeops,
    cluster as clutools,
    interpolation as interp
)
from .fold_algorithm import FoldingAlgorithm, ReductionMethod
from .marasco_merging import merge_seg_subtree

PI = math.pi
Cluster = clutools.Cluster

logger = logutils.getBasicLogger(
                    name='marasco', level=logging.DEBUG,
                    format="%(levelname)s@%(filename)s:%(lineno)s  %(message)s")

# Log to file
# fmtr = logging.Formatter('%(levelname)s:%(message)s @%(filename)s:%(lineno)s')
# fh = logging.FileHandler('reduce_marasco.log')
# fh.setFormatter(fmtr)
# logger.addHandler(fh)
# Log to stream
# ch = logging.StreamHandler(sys.stdout)
# ch.setFormatter(fmtr)
# logger.addHandler(ch)


################################################################################
# Interface implementations
################################################################################


class MarascoFolder(FoldingAlgorithm):
    """
    Marasco folding/collapsing algorithm.

    Original publication: Marasco, A., Limongiello, A. & Migliore, M. -
    Fast and accurate low-dimensional reduction of biophysically detailed 
    neuron models. Scientific Reports 2, (2012).
    """

    impl_algorithm = ReductionMethod.Marasco

    def __init__(self, algorithm):
        """
        @note   FoldReduction class is responsible for maintaining
                bi-directional association
        """
        self.reduction = None
        if algorithm != self.impl_algorithm:
            raise ValueError("MarascoFolder has no implementation for folding algorithm {}".format(algorithm))


    def preprocess_reduction(self):
        """
        Preprocess cell for Marasco reduction. Execute once before all
        folding passes.
        """
        pass


    def _preprocess_folding_pass(self):
        """
        Prepare next folding pass: assign topology information
        to each Section.

        (Implementation of interface declared in FoldingAlgorithm)
        """
        reduction = self.reduction

        root_ref = reduction._root_ref
        allsecrefs = reduction.all_sec_refs

        # Set properties used in calculation
        reduction.assign_new_sec_gids(root_ref, allsecrefs)
        redutils.subtree_assign_attributes(root_ref, allsecrefs, {'max_passes': 100})

        # Assign topology info (order, level, strahler number)
        clutools.assign_topology_attrs(root_ref, allsecrefs)
        reduction.fix_topology_below_roots()

        # Assign path properties
        f_lambda = self.reduction.get_reduction_param('f_lambda')
        for secref in allsecrefs:
            # Calculate path length, path resistance, electrotonic path length to each segment
            redutils.sec_path_props(secref, f_lambda, reduction.gleak_name)


    def fold_one_pass(self, i_pass, Y_criterion='outer_branchpoints'):
        """
        Do a folding pass.
        """
        logger.info("\n###############################################################"
                    "\nPreprocessing folding pass ..."
                    "\n###############################################################")

        self._preprocess_folding_pass()

        # Find collapsable branch points
        target_Y_secs = treeops.find_collapsable(
                                self.reduction.all_sec_refs, 
                                i_pass, Y_criterion)

        fold_pass = FoldingPass(
                        self.reduction, i_pass, 
                        fold_roots=target_Y_secs,
                        interp_prop='path_L',
                        interp_method='linear_neighbors',
                        gbar_scaling='area')
        
        new_refs = fold_pass.do_folds() # do collapse operation at each branch points

        # TODO: merge linear sections so next folding pass can collapse further
        #       OR BETTER: addapt merge_seg_subtree() so it jump over linear connection points
        # TODO: check what happens to non-uniform diameter

        return new_refs



    def postprocess_reduction(self):
        """
        Post-process cell after Marasco reduction. Execute once after all
        folding passes.
        """
        reduction = self.reduction
        all_sec_refs = reduction.all_sec_refs

        # Tweaking
        tweak_funcs = reduction.get_reduction_param('post_tweak_funcs')
        for func in tweak_funcs:
            func(reduction)

        # Assign identifiers (for synapse placement etc.)
        reduction.assign_new_sec_gids(reduction._root_ref, all_sec_refs)

        # Assign topology info (order, level, strahler number)
        for fold_root in reduction._fold_root_refs:
            clutools.assign_topology_attrs(fold_root, all_sec_refs)


################################################################################
# Private functions
################################################################################


class FoldingPass(object):
    """
    Encapsulate data for one folding pass (all folds in one reduction step)

    Member data
    -----------

    i_pass : int
        number of the folding pass

    clusters : list(cluster.Cluster)
        data describing a collapsed fork

    """

    def __init__(
            self,
            reduction,
            i_pass,
            fold_roots,
            interp_prop,
            interp_method,
            gbar_scaling
        ):
        self.i_pass = i_pass
        self.reduction = reduction
        self.fold_roots = fold_roots
        self.interp_prop = interp_prop
        self.interp_method = interp_method
        self.gbar_scaling = gbar_scaling

        self.f_lambda = self.reduction.get_reduction_param('f_lambda')
        self.gleak_name = self.reduction.gleak_name


    def do_folds(self):
        """
        Calculate folds, make equivalent sections, and insert them into tree.

        @return     list(SectionRef)
                    list of equivalent sections for folds
        """

        # Flag all sections as unvisited
        for secref in self.reduction.all_sec_refs:

            secref.is_substituted = False
            secref.is_deleted = False

            if not hasattr(secref, 'absorbed'):
                secref.absorbed = [False] * secref.sec.nseg
                secref.visited = [False] * secref.sec.nseg
            
            secref.zip_labels = [None] * secref.sec.nseg

        logger.info("\n###############################################################"
                    "\nCalculating folds ..."
                    "\n###############################################################")

        # Create equivalents: collapse (zip) each Y section up to length of first branch point
        clusters = []
        for j_zip, par_ref in enumerate(self.fold_roots):
            cluster = self.fold_subtree(par_ref, j_zip)
            self.calc_fold_statistics(cluster)
            clusters.append(cluster)

        # substitute equivalents
        clusters_eqrefs = {cluster: None for cluster in clusters}
        self.make_equivalent_secs(clusters_eqrefs)

        logger.info("\n###############################################################"
                    "\nDetermining active electrical properties ..."
                    "\n###############################################################")

        # Set their properties
        self.set_passive_params(clusters_eqrefs)
        self.set_conductances_interp(clusters_eqrefs)
        self.scale_conductances(clusters_eqrefs)

        logger.info("\n###############################################################"
                    "\nSubstituting new sections ..."
                    "\n###############################################################")

        # TODO FIXME: this must happen before interpolation?!
        self.sub_equivalent_secs(clusters_eqrefs)

        return clusters_eqrefs.values()


    def fold_subtree(self, par_ref, j_zip):
        """
        Calculate one fold and return equivalent properties

        @param  par_ref : SectionRef
                folding root (fork) whose children should be folded

        @param  j_zip : int
                sequence number of the fold within current folding pass

        @return cluster: Cluster
                cluster object describing properties of equivalent section
        """
        par_sec = par_ref.sec
        child_secs = par_sec.children()
        allsecrefs = self.reduction.all_sec_refs

        # Function to determine which segment can be 'zipped'
        min_child_L = min(sec.L for sec in child_secs) # Section is unbranched cable
        def can_absorb(seg, jseg, ref):
            """
            Can segment be absorbed? Yes if following criteria satisfied:
                - section is connected to fold root
                - segment is not farther from fold root than closest terminal
                  segment of any child section of the fold root
            
            @param  seg     <nrn.Segment> with x-value equal to segment midpoint
            """
            direct_child = (ref.parent.same(par_sec))
            not_too_far = (seg.x*seg.sec.L <= min_child_L) or jseg==0
            return direct_child and not_too_far
        
        # Name for equivalent zipped section
        # name_sanitized = par_sec.name().replace('[','').replace(']','').replace('.','_')
        name_sanitized = re.sub(r"[\[\]\.]", "", par_sec.name())
        alphabet_uppercase = [chr(i) for i in range(65,90+1)] # A-Z are ASCII 65-90
        zip_label = "zip{0}_{1}".format(alphabet_uppercase[self.i_pass], name_sanitized)
        zip_id = 1000*self.i_pass + j_zip

        # Function for processing zipped SectionRefs
        merged_sec_gids = set()
        merged_region_labels = set()
        def process_zipped_seg(seg, jseg, ref):
            # Tag segment with label of current zip operation
            ref.zip_labels[jseg] = zip_label
            # Save GIDs of original sections that are zipped/absorbed into equivalent section
            if ref.is_original:
                merged_sec_gids.add(ref.gid)
                merged_region_labels.add(ref.region_label)
            else:
                merged_sec_gids.update(ref.merged_sec_gids)
                merged_region_labels.update(ref.merged_region_labels)
        
        # Perform 'zip' operation
        far_bound_segs = [] # last (furthest) segments that are zipped
        eq_seq, eq_br = merge_seg_subtree(
                            par_sec(1.0), allsecrefs, can_absorb, 
                            process_zipped_seg, far_bound_segs)

        bounds_info = ["\n\t->| seg[{1}/{2}] @{0}".format(seg, seg_index(seg)+1, seg.sec.nseg) for seg in far_bound_segs]
        logger.debug("Folding @{2} up to next {0} boundary segments:{1}".format(len(far_bound_segs), "".join(bounds_info), par_sec.name()))
        logger.debug("Topology at folding root:\n"+subtree_topology(par_sec, max_depth=2))

        # Make Cluster object that represents collapsed segments
        cluster = Cluster(zip_label)
        cluster.eqL = eq_br.L_eq
        cluster.eqdiam = eq_br.diam_eq
        cluster.eq_area_sum = eq_br.L_eq * PI * eq_br.diam_eq
        cluster.eqri = eq_br.Ri_eq
        cluster.merged_sec_gids = merged_sec_gids
        cluster.merged_region_labels = merged_region_labels
        cluster.zip_id = zip_id
        cluster.parent_seg = par_sec(1.0)
        cluster.bound_segs = far_bound_segs # Save boundaries (for substitution)

        return cluster


    def calc_fold_statistics(self, cluster):
        """
        Calculate statistics on sections that are absorbed into cluster.

        @post   cluster has following attributes:

                or_area: float
                total area of original sections absorbed by cluster.

                or_cm: float
                mean specific capacitance in original sections, weighted by area

                or_cmtot: float
                total capacitance (non-specific) of original sections


        """
        
        glist = self.reduction.gbar_names
        allsecrefs = self.reduction.all_sec_refs
        
        # Gather all cluster sections & segments
        clu_secs = [secref for secref in allsecrefs if (cluster.label in secref.zip_labels)]
        clu_segs = [seg for ref in clu_secs for jseg,seg in enumerate(ref.sec) if (
                        ref.zip_labels[jseg]==cluster.label)]
        cluster.or_nseg = len(clu_segs)

        # Calculate max/min path length
        clu_path_L = [ref.pathL_seg[j] for ref in clu_secs for j,seg in enumerate(ref.sec) if (
                        ref.zip_labels[j]==cluster.label)]
        cluster.orMaxPathL = max(clu_path_L)
        cluster.orMinPathL = min(clu_path_L)

        # Calculate min/max axial path resistance
        clu_path_ri = [ref.pathri_seg[j] for ref in clu_secs for j,seg in enumerate(ref.sec) if (
                        ref.zip_labels[j]==cluster.label)]
        cluster.orMaxpathri = max(clu_path_ri)
        cluster.orMinpathri = min(clu_path_ri)

        # Calculate area, capacitance, conductances
        cluster.or_area = sum(seg.area() for seg in clu_segs)
        cluster.or_cmtot = sum(seg.cm*seg.area() for seg in clu_segs)
        cluster.or_cm = cluster.or_cmtot / cluster.or_area
        cluster.or_gtot = dict((gname, 0.0) for gname in glist)
        for gname in glist:
            try:
                gval = getattr(seg, gname, 0.0)
            except NameError: # NEURON error if mechanism not inserted in section
                gval = 0.0
            cluster.or_gtot[gname] += sum(gval*seg.area() for seg in clu_segs)

        # Equivalent axial resistance
        clu_segs_Ra = [seg.sec.Ra for seg in clu_segs]
        if min(clu_segs_Ra) == max(clu_segs_Ra):
            cluster.eqRa = clu_segs_Ra[0]
        else:
            logger.warning("Sections have non-uniform Ra, calculating average "
                            "axial resistance per unit length, weighted by area")
            cluster.eqRa = PI*(cluster.eqdiam/2.)**2*cluster.eqri*100./cluster.eqL # eq. Ra^eq

        # Calculate electrotonic path length
        cluster.or_L_elec = sum(electro.seg_L_elec(seg, self.gleak_name, self.f_lambda) 
                                for seg in clu_segs)
        
        cluster.eq_lambda = electro.calc_lambda_AC(
                                self.f_lambda, cluster.eqdiam, cluster.eqRa, 
                                cluster.or_cmtot/cluster.eq_area_sum)

        cluster.eq_L_elec = cluster.eqL/cluster.eq_lambda
        

        # Debug
        print_attrs = ['eqL', 'eqdiam', 'or_area', 'eq_area_sum', 'or_L_elec', 'eq_L_elec']
        clu_info = ("- {0}: {1}".format(prop, getattr(cluster, prop)) for prop in print_attrs)
        logger.anal("Properties of equivalent section for fold:\n\t{0}".format("\n\t".join(clu_info)))
    


    def make_equivalent_secs(self, clusters_refs):
        """
        Create equivalent section for each cluster (only geometry and connections).

        @param  clusters_refs: dict(Cluster: object)
                (passed by reference) dict to fill with equivalent SectionRef
        """

        for cluster in clusters_refs.iterkeys():

            # Create equivalent section
            if cluster.label in [sec.name() for sec in h.allsec()]:
                raise Exception('Section named {} already exists'.format(cluster.label))
            
            created = h("create %s" % cluster.label)
            if created != 1:
                raise Exception("Could not create section with name '{}'".format(cluster.label))
            
            eqsec = getattr(h, cluster.label)
            eqref = ExtSecRef(sec=eqsec)

            # Set passive properties
            eqsec.L = cluster.eqL
            eqsec.diam = cluster.eqdiam

            # Connect to tree (need to trace path from soma to section)
            eqsec.connect(cluster.parent_seg, 0.0) # see help(sec.connect)

            # Save metadata
            for prop in ['zip_id', 'or_area', 'merged_sec_gids', 'merged_region_labels']:
                setattr(eqref, prop, getattr(cluster, prop)) # copy some useful attributes
            eqref.is_original = False

            clusters_refs[cluster] = eqref


    def set_passive_params(self, clusters_refs):
        """
        Set passive parameters and insert mechanisms into equivalent sections.

        @param  clusters_refs: dict(Cluster: SectionRef)
                clusters representing collapsed branches and their equivalent section
        """
        
        # Set passive electrical properties
        for cluster, eqref in clusters_refs.items():
            eqsec = eqref.sec

            # Scale passive electrical properties
            cluster.eq_area = sum(seg.area() for seg in eqsec) # should be same as cluster eqSurf
            area_ratio = cluster.or_area / cluster.eq_area
            logger.debug("Surface area ratio is %f" % area_ratio)

            # Area scaling of passive electrical properties
            eq_cm = cluster.or_cmtot / cluster.eq_area # more accurate than cm * or_area/eq_area
            or_gleak = cluster.or_gtot[self.gleak_name] / cluster.or_area
            eq_gleak = or_gleak * area_ratio # same as reducing Rm by area_new/area_old

            # Set number of segments based on rule of thumb electrotonic length
            eqsec.Ra = cluster.eqRa
            eqsec.cm = eq_cm
            eqsec.nseg = electro.calc_min_nseg_hines(
                                    self.f_lambda, eqsec.L, eqsec.diam, 
                                    eqsec.Ra, eq_cm)

            logger.debug("Fold reduces number of segments saved by {0} (Hines rule) and reduces L/lambda by {1:.2f} \n".format(cluster.or_nseg-eqsec.nseg, cluster.eq_L_elec/cluster.or_L_elec))

            # Copy section mechanisms and properties
            absorbed_secs = redutils.find_secprops(
                                    self.reduction.orig_tree_props,
                                    lambda sec: sec.gid in eqref.merged_sec_gids)
            
            redutils.merge_sec_properties(
                            absorbed_secs, eqsec, 
                            self.reduction.mechs_params_nogbar, 
                            check_uniform=True)

            # Save Cm and conductances for each section for reconstruction
            cluster.nseg = eqsec.nseg # save for reconstruction
            cluster.eq_gbar = dict((
                (gname, [float('NaN')]*cluster.nseg) for gname in self.reduction.gbar_names))
            cluster.eq_cm = [eq_cm]*cluster.nseg

            # Set Cm and gleak (Rm) for each segment
            if self.gbar_scaling is not None:
                for j, seg in enumerate(eqsec):
                    setattr(seg, self.gleak_name, eq_gleak)
                    cluster.eq_gbar[self.gleak_name][j] = eq_gleak


    def set_conductances_interp(self, clusters_refs):
        """
        Set channel conductances by interpolating conductance values along
        a path of sections in the original cell model.

        @param  clusters_refs: dict(Cluster: SectionRef)
                clusters representing collapsed branches and their equivalent section

        @post   for each conductance in reduction.active_gbar_names, all segments
                in the cluster are assigned a value true interpolation.
        """
        if self.interp_prop == 'path_L':
            seg_prop = 'pathL_seg'
        elif self.interp_prop == 'path_ri':
            seg_prop = 'pathri_seg'
        elif self.interp_prop == 'path_L_elec':
            seg_prop = 'pathL_elec'
        else:
            raise ValueError("Unknown path property '{}'".format(self.interp_prop))

        for cluster, eqref in clusters_refs.items():
            logger.debug("Scaling properties of cluster %s ..." % cluster.label)
            eqsec = eqref.sec

            # Calculate path lengths in equivalent section
            redutils.sec_path_props(eqref, self.f_lambda, self.gleak_name)

            # Get path of original sections running through one of folded branches
            path_secs = self.reduction.get_interpolation_path_secs(eqref)

            # Find conductances at same path length (to each segment midpoint) in original cell
            for j_seg, seg in enumerate(eqsec):
                
                # Get adjacent segments along interpolation path
                path_L = getattr(eqref, seg_prop)[j_seg]
                bound_segs, bound_L = interp.find_adj_path_segs(self.interp_prop, path_L, path_secs)
                
                # DEBUG STATEMENTS:
                # bounds_info = "\n".join(("\t- bounds {0} - {1}".format(a, b) for a,b in bound_segs))
                # logger.debug("Found boundary segments at path length "
                #              "x={0:.3f}:\n{1}".format(path_L, bounds_info))

                # INTERPOLATE: Set conductances by interpolating neighbors
                for gname in self.reduction.active_gbar_names:

                    if not hasattr(eqsec, gname):
                        cluster.eq_gbar[gname][j_seg] = 0.0
                        continue
                    
                    if self.interp_method == 'linear_neighbors':
                        gval = interp.interp_gbar_linear_neighbors(path_L, gname, bound_segs, bound_L)
                    
                    else:
                        match_method = re.search(r'^[a-z]+', self.interp_method)
                        method = match_method.group() # should be nearest, left, or right
                        gval = interp.interp_gbar_pick_neighbor(path_L, gname, 
                                            bound_segs[0], bound_L[0], method)
                    
                    setattr(seg, gname, gval)
                    cluster.eq_gbar[gname][j_seg] = gval


    def scale_conductances(self, clusters_refs):
        """
        After interpolation of specific conductance values, scale values
        so that total non-specific conductance integrated over cluster is the same.

        @param  clusters_refs: dict(Cluster: SectionRef)
                clusters representing collapsed branches and their equivalent section

        @post   
        """

        for cluster, eqref in clusters_refs.items():
            eqsec = eqref.sec

            # Re-scale gbar distribution to yield same total gbar (sum(gbar*area))
            for gname in self.reduction.active_gbar_names:
                if not hasattr(eqsec, gname):
                    continue
                
                # original and current total conductance
                or_gtot = cluster.or_gtot[gname]
                eq_gtot = sum(getattr(seg, gname, 0.0)*seg.area() for seg in eqsec)
                if eq_gtot <= 0.:
                    eq_gtot = 1.
                
                
                for j_seg, seg in enumerate(eqsec):
                    
                    if self.gbar_scaling == 'area':
                        # conserves ratio in each segment but not total original conductance
                        scale = cluster.or_area/cluster.eq_area
                    
                    elif self.gbar_scaling == 'gbar_integral':
                        # does not conserve ratio but conserves gtot_or since: sum(g_i*area_i * or_area/eq_area) = or_area/eq_area * sum(gi*area_i) ~= or_area/eq_area * g_avg*eq_area = or_area*g_avg
                        scale = or_gtot/eq_gtot
                    
                    else:
                        raise Exception("Unknown gbar scaling method'{}'.".format(self.gbar_scaling))
                    
                    # Set gbar
                    gval = getattr(seg, gname) * scale
                    setattr(seg, gname, gval)
                    
                    cluster.eq_gbar[gname][j_seg] = gval # save for reconstruction

            # Debugging info:
            logger.anal("Created equivalent Section '%s' with \n\tL\tdiam\tcm\tRa\tnseg"
                         "\n\t%.3f\t%.3f\t%.3f\t%.3f\t%d\n", cluster.label, 
                         eqsec.L, eqsec.diam, eqsec.cm, eqsec.Ra, eqsec.nseg)


    def sub_equivalent_secs(self, clusters_refs):
        """
        Substitute / insert fold-equivalent sections into tree.
        """

        orsecrefs = self.reduction.all_sec_refs

        for cluster, eqref in clusters_refs.items():
            eqsec = eqref.sec

            # Disconnect substituted segments and attach segment after Y boundary
            # Can only do this now since all paths need to be walkable before this
            logger.debug("Substituting zipped section {0} "
                         "at proximal seg {1} and distal segments {2}".format(
                         eqsec, cluster.parent_seg, cluster.bound_segs))
            
            treeops.sub_equivalent_Y_sec(
                        eqsec, cluster.parent_seg, cluster.bound_segs, 
                        orsecrefs, self.reduction.mechs_gbars_dict,
                        delete_substituted=True)

        # build new list of valid SectionRef
        newsecrefs = [ref for ref in orsecrefs if not (ref.is_substituted or ref.is_deleted)]
        newsecrefs.extend(clusters_refs.values())