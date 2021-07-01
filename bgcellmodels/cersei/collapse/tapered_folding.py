"""
Morphology reduction by constructing equivalent cylinders with variable
diameter.

@author Lucas Koelman
@date   19-01-2017
"""

# Python modules
import re
import math

# NEURON modules
import neuron
h = neuron.h
h.load_file("stdlib.hoc") # Load the standard library

# Own modules
from bgcellmodels.common import logutils
from bgcellmodels.common.treeutils import subtree_topology
from bgcellmodels.common.nrnutil import ExtSecRef, getsecref
from bgcellmodels.common.stdutil import isclose
from bgcellmodels.common.electrotonic import calc_min_nseg_hines

from . import (
    cluster,
    redutils,
    tree_edit,
    fold_algorithm,
    interpolation as interp,
    tapered_merging as taper
)

PI = math.pi
logger = logutils.getLogger(__name__)


class TaperedFolder(fold_algorithm.FoldingAlgorithm):
    """
    Tapered folding procedure.

    Constructs equivalent cylinders with variable diameter by marching along
    parallel branches concurrently and, at each step, calculating an equivalent
    cylinder that is connected sequentially to the last.
    """

    def __init__(self, algorithm):
        """
        @note   FoldReduction class is responsible for maintaining
                bi-directional association
        """
        self.reduction = None
        self.impl_algorithm = algorithm


    def preprocess_reduction(self):
        """
        Preprocess cell for reduction.
        """
        cluster.assign_topology_attrs(self.reduction._root_ref, self.reduction.all_sec_refs)


    def fold_one_pass(self, i_pass, Y_criterion='use_fold_roots'):
        """
        Collapse branches at branch points identified by given criterion.
        """
        # Use user-provided folding roots, don't look for forks ourselves
        if Y_criterion != 'use_fold_roots':
            raise ValueError(Y_criterion)
        target_Y_secs = self.reduction._fold_root_refs

        settings = self.reduction.reduction_params

        collapser = AdvancingFrontCollapse(
                        self.reduction, 
                        fold_roots=target_Y_secs,
                        **settings)
        
        new_refs = collapser.collapse() # do collapse operation at each branch points

        return new_refs


    def postprocess_reduction(self):
        """
        Postprocess cell after reduction procedure. Executed once.
        """
        pass


################################################################################
# Private functions
################################################################################


class AdvancingFrontCollapse(object):
    """
    Encapsulate data for collapse operation.

    @param      gbar_init_method : str

                Method used for setting initial values of specific membrane
                conductances and passive electrical properties (cm, gleak).
                Valid values are:

                    - 'area_weighted_average': average value of absorbed
                      cylinders, weighted by their area

                    - 'interp_<METHOD>': interpolate values found at equivalent
                      path length in the full model

    @param      gbar_scale_method : str

                Method to scale specific membrane conductances. Valid values are:
                    - None/'None'
                    - 'surface_area_ratio'
                    - 'gbar_integrated_cylinder_ratio'
                    - 'gbar_integrated_subtree_ratio'

    @param      passive_scale_method : str

                Method to scale specific passive electrical parameters: membrane
                capacitance (sec.cm) and leak conductance (sec.gleak).
                Valid values: 'match_input_impedance_soma',
                'match_input_impedance_subtrees' or one of accepted values of
                gbar_scale_method


    WAYS TO SET GBAR
    ----------------

    - set them to the area-weighted average gbar of all cylinders
      that have been merged into the equivalent sequential section, with
      no further scaling

        + gbar_scale_method = None / 'None'

        + i.e. g_i = g_wavg = sum(g_j*S_j) / sum(S_j) where j are absorbed cylinders

        +  => Gtot = sum(g_wavg * Snew_cyl)
        
        + this is likely to preserve the gradient/spatial profile since
          gbar is the weighted average at each distance
        
        + this does NOT compensate for area reduction, so the total
          integrated Gbar in an equivalent section will not be the same
          as that in the original cylinders
    

    - initialize to area-weighted gbar, then multiply the gbar in each
      equivalent section with the ration of original area over new area

        + gbar_scale_method = 'surface_area_ratio'

        + i.e. g_i = g_wavg * Sorig_cyl / Snew_cyl

        + => Gtot = sum(g_wavg * Snew_cyl * Sorig_cyl / Snew_cyl)
                  = sum(g_wavg * Sorig_cyl)
                  ~= sum(Gorig_cyl)
                  ~= Gorig_subtree

        + since each cylinder has its own scaling factor, this will likely
          NOT preserve the gradient/spatial profile

        + this _approximately_ preserves Gbar in each cylinder and the
          entire subtree


    - initialize to area-weighted gbar, then multiply the gbar in each
      equivalent section with the ratio of integrated Gbar over that
      section vs over absorbed cylinders

        + gbar_scale_method = 'gbar_integrated_cylinder_ratio'

        + i.e. g_i = g_wavg * Gorig_cyl / Gnew_cyl

        + => Gtot = sum(g_wavg * Snew_cyl * Gorig_cyl / Gnew_cyl)
                  = sum(Gorig_cyl)
                  = Gorig_subtree

        + since each cylinder has its own scaling factor, this will likely
          NOT preserve the gradient/spatial profile

        + this _exactly_ preserves total Gbar in each cylinder AND over,
          entire subtree so this has good area compensation


    - initialize to area-weighted gbar, then multiply the gbar in each
      equivalent section with the ratio of integrated Gbar over entire
      subtree

        + gbar_scale_method = 'gbar_integrated_subtree_ratio'

        + i.e. g_i = g_wavg * Gorig_subtree / Gnew_subtree

        + => Gtot = sum(g_wavg * Snew_cyl * Gorig_subtree / Gnew_subtree)
                  = Gorig_subtree / Gnew_subtree * sum(g_wavg * Snew_cyl)
                  = Gorig_subtree

        + this is likely to preserve the gradient/spatial profile since all
          the g_wavg are multiplied by constant factor

        + this preserves total Gbar over entire subtree, but NOT in each
          cylinder individually, which is a form of area compensation

    """

    def __init__(
            self,
            reduction,
            **kwargs
        ):
        self.reduction = reduction
        self.fold_roots = kwargs.pop('fold_roots')

        # Settings for scaling of electrical parameters
        self.gbar_init_method = kwargs['gbar_init_method']
        self.gbar_scale_method = kwargs['gbar_scale_method']
        self.passive_scale_method = kwargs['passive_scale_method']
        
        if self.gbar_init_method.startswith('interp'):
            self.interp_prop = kwargs['interp_prop']
        
        if self.passive_scale_method.startswith('match_input_impedance'):
            self.impedance_linearize_gating = kwargs['Z_linearize_gating']
            self.match_Zin_frequency = kwargs.get('match_input_impedance_frequency', 0.0)

        self.f_lambda = self.reduction.get_reduction_param('f_lambda')
        self.gleak_name = self.reduction.gleak_name


    def collapse(self):
        """
        Calculate folds, make equivalent sections, and insert them into tree.

        @return     list(SectionRef)
                    list of equivalent sections for folds
        """

        # Flag all sections as unvisited
        for secref in self.reduction.all_sec_refs:
            secref.is_original = True
            secref.is_substituted = False
            secref.is_deleted = False
            if not hasattr(secref, 'absorbed'):
                secref.absorbed = [False] * secref.sec.nseg
                secref.visited = [False] * secref.sec.nseg

        # Measure input impedance first
        if self.passive_scale_method.startswith('match_input_impedance'):
            self.measure_input_impedance()

        logger.info("\n###############################################################"
                    "\nCollapsing subtrees ..."
                    "\n###############################################################")

        # Calculate equivalent cylinders for each folding root
        subtree_eq_cyls = []
        each_subtree_refs = []
        for collapse_index, root_ref in enumerate(self.fold_roots):

            seq_cyls = self.collapse_subtree(root_ref, collapse_index)
            subtree_eq_cyls.append(seq_cyls)

            # TODO: refactor from here until end of for-loop
            # Make NEURON Section objects and connect them
            sequential_refs = self.make_equivalent_secs(seq_cyls)
            each_subtree_refs.append(sequential_refs)
            self.insert_mechanisms(sequential_refs)

            logger.info("\n###############################################################"
                        "\nSet electrical properties ..."
                        "\n###############################################################")

            # Set initial values of electrical properties
            self.init_passive_wavg(sequential_refs)
            
            if self.gbar_init_method == 'area_weighted_average':
                self.init_conductances_wavg(sequential_refs)
            
            elif self.gbar_init_method.startswith('interp'):
                self.init_conductances_interp(sequential_refs)

            # Compensate for area reduction
            self.scale_specific_elec_params(sequential_refs)

            logger.info("\n###############################################################"
                        "\nDelete substituted sections ..."
                        "\n###############################################################")

            self.disconnect_substituted_secs(root_ref, delete=True)

            # Debug info
            logger.debug("Topology at root {0} after collapse:\n".format(
                        root_ref.sec.name()) + subtree_topology(root_ref.sec, max_depth=2))

            # Fit input impedance after each subtree substitution
            if self.passive_scale_method == 'match_input_impedance_subtrees':
                self.fit_passive_params([root_ref], [sequential_refs])

        # Fit input impedance if requested
        if self.passive_scale_method == 'match_input_impedance_soma':
            self.fit_passive_params(self.fold_roots, each_subtree_refs)

        new_sec_refs = [ref for subtree in each_subtree_refs for ref in subtree]
        return new_sec_refs


    def int_to_alphabet(self, i):
        """
        Convert integer to alphabetic character.
        """
        # NOTE: A-Z are ASCII characters 65-90
        num_char = 90-65+1
        q, r = divmod(i, num_char)
        letter = chr(65+r)
        return (q+1)*letter


    def collapse_subtree(self, par_ref, i_collapse):
        """
        Calculate one subtree and return equivalent cylinders

        @param  par_ref : SectionRef
                folding root (fork) whose children should be folded

        @param  i_collapse : int
                index (sequence number_ of collapse operation

        @return cluster: Cluster
                cluster object describing properties of equivalent section
        """
        par_sec = par_ref.sec
        allsecrefs = self.reduction.all_sec_refs

        logger.debug("Topology at folding root {0}:\n".format(par_sec.name()) + 
                        subtree_topology(par_sec, max_depth=2))

        # TODO: check in this file and make everything compatible with nseg > 1
        # TODO: - determine nseg based on split_criterion and split_dX 
        #       - or based on actual lambda, but prevent increase in nseg through rounding
        
        # Collapse operation: calculate cylinders
        merge_walk = taper.MergingWalk(
                            par_sec(1.0), allsecrefs,
                            **self.reduction.reduction_params)
                            
        seq_cyls = merge_walk.merge_cylinders_subtree()

        # Label this collapse operation
        zip_id = 1000 + i_collapse # usable as gid

        # Save metadata for each equivant cylinder
        for i_cyl, cyl in enumerate(seq_cyls):
            cyl.zip_id = zip_id
            cyl.label = "zip{0}{1}".format(self.int_to_alphabet(i_collapse), i_cyl)
            if i_cyl == 0:
                cyl.parent_seg = par_sec(1.0)
            else:
                cyl.parent_seg = None

        # Post-process subtree
        subtree_seclist = h.SectionList()
        subtree_seclist.subtree(sec=par_sec)
        for child_sec in subtree_seclist: # iterates CAS
            # Make sure every Section has been visited and fully absorbed (until 1-end)
            child_ref = getsecref(child_sec, allsecrefs)
            assert child_ref.visited and child_ref.absorbed

        return seq_cyls


    def make_equivalent_secs(self, subtree_eq_cyls):
        """
        Create equivalent section for each cylinder (only geometry and connections).

        @param  subtree_eq_cyls : list(Cylinder)
                For one collapsed subtree: list of equivalent sequential Cylinders
        """
        # Process all sequential cylinders of one collapsed subtree
        sequential_refs = []
        for cyl in subtree_eq_cyls:

            if h.section_exists(cyl.label):
                raise Exception('Section named {} already exists'.format(cyl.label))
            
            created = h("create %s" % cyl.label)
            if created != 1:
                raise Exception(
                        "Could not create section with name '{}'".format(cyl.label))
            
            eqsec = getattr(h, cyl.label)
            eqref = ExtSecRef(sec=eqsec)
            eqref.label = cyl.label
            eqref.gid = cyl.zip_id

            # Set passive properties (nseg determined after cm)
            eqsec.L = cyl.L
            eqsec.diam = cyl.diam
            eqsec.Ra = cyl.Ra
            assert isclose(cyl.area, sum(seg.area() for seg in eqsec), rel_tol=1e-9)

            # Connect to tree (need to trace path from soma to section)
            if cyl.parent_seg is not None:
                eqsec.connect(cyl.parent_seg, 0.0)
            else:
                eqsec.connect(sequential_refs[-1].sec(1.0), 0.0)

            # Save info about merged properties
            eqref.orig_props = cyl.orig_props
            eqref.is_original = False
            sequential_refs.append(eqref)

            # Debugging info
            logger.anal("Created equivalent Section '%s' with \n\tL\tdiam\tRa"
                         "\n\t%.3f\t%.3f\t%.3f\t%d\n", eqsec.name(), 
                         eqsec.L, eqsec.diam, eqsec.Ra)

        return sequential_refs


    def insert_mechanisms(self, subtree_sec_refs):
        """
        Set passive parameters and insert mechanisms into equivalent sections.

        @param  subtree_sec_refs : list(SectionRef)
                For one collapsed subtree: list of sequential equivalent Sections
        """
        for eqref in subtree_sec_refs:
            
            eqsec = eqref.sec
            orig = eqref.orig_props

            # Copy section mechanisms and properties
            absorbed_secs = redutils.find_secprops(
                                    self.reduction.orig_tree_props,
                                    lambda sec: sec.gid in orig.merged_sec_gids)
            
            redutils.merge_sec_properties(
                            absorbed_secs, eqsec, 
                            self.reduction.mechs_params_nogbar, 
                            check_uniform=True)


    def init_passive_wavg(self, subtree_sec_refs, adjust_nseg=True):
        """
        Set initial values of specific passive electrical properties
        to the area-weighted average of the values in absorbed cylinders.

        @param  subtree_sec_refs : list(SectionRef)
                For one collapsed subtree: list of sequential equivalent Sections
        """
        for eqref in subtree_sec_refs:

            eqsec = eqref.sec
            orig = eqref.orig_props
            assert all((seg.diam == eqsec(0.5).diam for seg in eqsec)) # assume uniform diameter

            # Set specific capacitance and leak conductance in each section
            eq_cm = orig.cmtot / orig.area
            eq_gleak = orig.gtot[self.gleak_name] / orig.area

            if adjust_nseg:
                eqsec.nseg = calc_min_nseg_hines(
                                    self.f_lambda, eqsec.L, eqsec.diam, 
                                    eqsec.Ra, eq_cm, round_up=False)

            # Set RANGE properties after nseg
            setattr(eqsec, 'cm', eq_cm)
            setattr(eqsec, self.gleak_name, eq_gleak)
            # TODO PROBLEM: cm and gbar are uniform across entire cylinder. One method to solve this would be to first make Cylinder objects with smaller step size, so you have orig_props with finer resolution, then merge these sequentially but keep orig_props for each segment cylinder (!for area ratios, would need to use area of new segment cylinders, not original segment cylinders!)


    def init_conductances_wavg(self, subtree_sec_refs):
        """
        Set initial values of active membrane conductances
        to the area-weighted average of the values in absorbed cylinders.
        
        @param  subtree_sec_refs : list(SectionRef)
                For one collapsed subtree: list of sequential equivalent Sections
        """
        for eqref in subtree_sec_refs:

            eqsec = eqref.sec
            orig = eqref.orig_props

            # All active membrane conductances: set uniform in cylinder
            for gname in self.reduction.active_gbar_names:

                eq_gbar = orig.gtot[gname] / orig.area

                if hasattr(eqsec, gname):
                    # assign same value to all segments
                    setattr(eqsec, gname, eq_gbar)
                else:
                    # if mechanism not inserted: assume that this is because none of the absorbed cylinders had the mechanism present
                    assert isclose(eq_gbar, 0.0, rel_tol=1e-9)


    def init_conductances_interp(self, subtree_sec_refs):
        """
        Set initial values of active membrane conductances by interpolating 
        conductance values along a path of sections in the original cell model.

        @param  subtree_sec_refs : list(SectionRef)
                For one collapsed subtree: list of sequential equivalent Sections
        """
        # Names of path-integrated properties stored on SectionRef
        path_prop_names = {
            'path_L': 'pathL_seg',
            'path_ri': 'pathri_seg',
            'path_L_elec': 'pathL_elec',
        }
        seg_prop = path_prop_names[self.interp_prop]

        for eqref in subtree_sec_refs:
            
            logger.debug("Scaling properties of cylinder %s ..." % eqref.label)
            eqsec = eqref.sec

            # Calculate path-integrated properties
            redutils.sec_path_props(eqref, self.f_lambda, self.gleak_name)

            # Find path of original sections running through one of folded branches
            path_secs = self.reduction.get_interpolation_path_secs(eqref)

            # Find conductances at same path length (to each segment midpoint) in original cell
            for j_seg, seg in enumerate(eqsec):
                
                # Get adjacent segments along interpolation path
                path_L = getattr(eqref, seg_prop)[j_seg] # value of interpolated property
                bound_segs, bound_L = interp.find_adj_path_segs(self.interp_prop, path_L, path_secs)

                # INTERPOLATE: Set conductances by interpolating neighbors
                for gname in self.reduction.active_gbar_names:
                    if not hasattr(seg, gname):
                        continue
                    
                    # get interpolation function and use it to compute gbar
                    if self.gbar_init_method == 'linear_neighbors':
                        gval = interp.interp_gbar_linear_neighbors(path_L, gname, bound_segs, bound_L)
                    else:
                        match_method = re.search(r'^[a-z]+', self.gbar_init_method)
                        method = match_method.group() # should be 'nearest', 'left', or 'right'
                        gval = interp.interp_gbar_pick_neighbor(
                                            path_L, gname, 
                                            bound_segs[0], bound_L[0],
                                            method)
                    
                    setattr(seg, gname, gval)


    def scale_specific_elec_params(self, subtree_sec_refs):
        """
        After interpolation of specific conductance values, scale values
        so that total non-specific conductance integrated over cluster is the same.

        @param  subtree_sec_refs : list(SectionRef)
                For one collapsed subtree: list of sequential equivalent Sections

        """
        g_range_names = self.reduction.gbar_names + ['cm']

        # For scaling method C: look-ahead and calculate integrated conductance over entire subtree
        subtree_gtot = {gname: 0.0 for gname in g_range_names}
        if self.gbar_scale_method == 'gbar_integrated_subtree_ratio':
            for eqref in subtree_sec_refs:
                for gname in g_range_names:
                    subtree_gtot[gname] += sum(getattr(seg, gname, 0.0)*seg.area() for seg in eqref.sec)

        for eqref in subtree_sec_refs:

            eqsec = eqref.sec
            assert eqsec.nseg == 1
            orig = eqref.orig_props
            orig.gtot['cm'] = orig.cmtot

            eq_area = sum(seg.area() for seg in eqsec)
            # logger.debug("Surface area ratio is %f" % (orig.area / eq_area))

            # Re-scale gbar distribution to yield same total gbar (sum(gbar*area))
            for gname in g_range_names:
                if not hasattr(eqsec, gname):
                    continue

                if gname in ['cm', self.gleak_name]:
                    requested_scale_method = self.passive_scale_method
                else:
                    requested_scale_method = self.gbar_scale_method
                
                # Determine scale factor (uniform in Section, like diameter)
                if ((requested_scale_method is None) or 
                    (requested_scale_method in ['None', 'none'])):
                    continue
                elif requested_scale_method == 'surface_area_ratio':
                    scale_factor = orig.area / eq_area
                elif requested_scale_method == 'gbar_integrated_cylinder_ratio':
                    g_cyl_integr = sum(getattr(seg, gname, 0.0)*seg.area() 
                                            for seg in eqsec)
                    g_orig_integr = orig.gtot[gname]
                    scale_factor = g_orig_integr / g_cyl_integr
                elif (requested_scale_method.startswith('match_input_impedance')
                        or requested_scale_method=='gbar_integrated_subtree_ratio'):
                    continue # do not scale now
                else:
                    raise ValueError("Unknown scaling method '{}'".format(
                            self.gbar_scale_method))
                
                # Set segment property
                for seg in eqsec:
                    g_init = getattr(seg, gname)
                    g_scaled = g_init * scale_factor
                    setattr(seg, gname, g_scaled)


    def measure_input_impedance(self, disconnect=False):
        """
        Measure input impedance of each subtree individually.

        @post   each SectionRef in self.fold_roots will have the attribute
                'subtree_Zin' set on its 'orig_props' object.
        """
        self.reduction.init_cell_steadystate()

        measure_points = [self.reduction._root_ref] + self.fold_roots
        for root_ref in measure_points:

            # First disconnect it
            root_sec = root_ref.sec
            parent_seg = root_sec.parentseg()
            if disconnect:
                h.disconnect(sec=root_sec)

            # Measure Zin at DC
            probe = h.Impedance()
            probe.loc(0.0, sec=root_sec) # injection site
            probe.compute(self.match_Zin_frequency, int(self.impedance_linearize_gating))
            subtree_Zin_DC = probe.input(0.0, sec=root_sec)
            root_ref.subtree_Zin_DC = subtree_Zin_DC

            # Reconnect
            if disconnect:
                root_sec.connect(parent_seg)

                    
    def fit_passive_params(self, root_refs, each_subtree_refs):
        """
        Fit passive params using measurement of input resistance in each
        subtree individually.

        Leak conductance is set so input resistance at DC matches original 
        subtree, and membrane capacitance is set to preserve time constant in 
        each cylinder (algorithm by Bush & Sejnowski, 1993), i.e.

            + Zin_orig_subtree(f=DC) == Zin_equiv_subtree(f=DC)
            
            + for cylinder in subtree:
                 Cm_old * Rm_old == Cm_new * Rm_new
        
        @post   for each collapsed subtree, sec.cm and sec.gleak are adjusted
                so that Zin measured at 0-end of subtree root is equal to
                value before reduction, within tolerance of 1e-6.
        """
        self.reduction.init_cell_steadystate()


        # Get folding root and its first child section
        root_sec = self.reduction._root_ref.sec
        subtree_sec_refs = [ref for subtree in each_subtree_refs for ref in subtree]

        # Initialize optimization
        g_range_names = self.reduction.gbar_names + ['cm']
        for ref in subtree_sec_refs: # save initial values of conductances and capacitance
            for gname in g_range_names:
                if hasattr(ref.sec, gname):
                    setattr(ref, gname+'_init', [getattr(seg, gname) for seg in ref.sec])
        
        probe = h.Impedance()
        probe.loc(0.0, sec=root_sec) # injection site
        probe.compute(self.match_Zin_frequency, int(self.impedance_linearize_gating))
        initial_Zin = probe.input(0.0, sec=root_sec)

        # Cost function for optimization
        target_Zin = self.reduction._root_ref.subtree_Zin_DC
        def cost_fun(param_hvec):
            """ Cost function/error to minimize """
            # Adjust leak conducance
            scale_factor = param_hvec.x[0]
            for ref in subtree_sec_refs:

                # Only scale gleak and cm
                # for i, seg in enumerate(ref.sec):
                #     new_gleak = ref.gleak_init[i] * scale_factor
                #     new_cm = ref.cm_init[i] * scale_factor # preserve time constant R*C
                #     setattr(seg, self.gleak_name, new_gleak)
                #     setattr(seg, 'cm', new_cm)

                # Scale all conductances
                for gname in g_range_names:
                    if not hasattr(ref.sec, gname):
                        continue
                    for i, seg in enumerate(ref.sec):
                        g_init = getattr(ref, gname+'_init')[i]
                        g_scaled = g_init * scale_factor
                        setattr(seg, gname, g_scaled)

            # Measure input inpedance
            probe.compute(self.match_Zin_frequency, int(self.impedance_linearize_gating))
            this_Zin = probe.input(0.0, sec=root_sec)
            return (this_Zin - target_Zin)**2 + float(scale_factor < 0) * 1e9

        # Parameters for optimization
        param_tolerance = 1e-3 # praxis attempt to return cost(x) such that if x0 is the true local minimum then norm(x-x0) < tolerance
        max_step = 10.0 # maximum distance from initial scale factor to optimal one
        praxis_verbosity = 1
        h.attr_praxis(param_tolerance, max_step, praxis_verbosity)
        
        # Optimize
        scale_init = 1.0
        param_vec = h.Vector([scale_init])
        Zin_err_sqrd = h.fit_praxis(cost_fun, param_vec)

        # Report results
        Zin_final = probe.input(0.0, sec=root_sec)
        Zin_error = math.sqrt(Zin_err_sqrd)
        scale_fitted = param_vec.x[0]
        if scale_fitted < 0.0:
            raise Exception("Got negative scale factor!")
        if not isclose(Zin_final, target_Zin, rel_tol=.05):
            logger.warning("Could not fit input impedance to within tolerance")

        print("""\
Optimization finished:
scale_fitted = {}
Zin_orig     = {}
Zin_initial  = {}
Zin_final    = {}
Zin error    = {}""".format(scale_fitted, target_Zin, initial_Zin, Zin_final, Zin_error))


    def disconnect_substituted_secs(self, fold_root, delete=True):
        """
        Disconnect substituted sections and delete them if requested.
        """
        def section_is_substituted(sec_ref):
            """ Check if Section was substituted and should be deleted. """
            return sec_ref.absorbed

        tree_edit.disconnect_subtree(
                    fold_root,
                    self.reduction.all_sec_refs,
                    should_disconnect=section_is_substituted,
                    delete=True)