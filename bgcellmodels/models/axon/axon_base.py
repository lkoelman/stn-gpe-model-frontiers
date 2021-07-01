"""
Base class for axon builders.
"""

from __future__ import division # float division for int literals, like Hoc
import math
import re
import logging

import numpy as np
import neuron

from bgcellmodels.morphology import morph_3d
from bgcellmodels.common import logutils, nrnutil
from bgcellmodels.common.stdutil import dotdict

from transforms3d import axangles

h = neuron.h
PI = math.pi

logger = logging.getLogger('AxonBuilder')
logger.setLevel(logging.DEBUG)
logging.basicConfig(format=logutils.DEFAULT_FORMAT)

def normvec(a):
    """
    Normalized vector pointing from a to b.
    """
    return a / np.sqrt(np.dot(a, a))


def veclen(a):
    """
    Vector length of a (Euclidean norm).
    """
    return np.sqrt(np.dot(a, a))


class AxonBuilder(object):
    """
    Base class for axon model

    The attributes that must be set are described below. You can either set
    them in a method or define them as class attributes.

    Attributes
    ----------

    @attr   'initial_comp_sequence' : list(str)
            Sequence of compartment types that define initial structure of axon

    @attr   'repeating_comp_sequence' : list(str)
            Sequence of compartment types that define repeating structure of axon

    @attr   'nodal_compartment_type' : str
            Name of compartment type representing node of Ranvier

    @attr   'compartment_defs' : dict[str, dict[str, object]]
            Dictionary mapping names of compartment types to the keys
            - 'mechanisms' : dict[<mechanism name> , <dict of parameters>]
            - 'passive' : dict[<parameter name>, <value>]
            - 'morphology': dict[<parameter name>, <value>]
    """

    def __init__(
            self,
            streamline_coords,
            interp_method='cartesian',
            tolerance_mm=1e-6,
            parent_cell=None,
            parent_sec=None,
            without_extracellular=False,
            termination_method='terminal_sequence',
            unmyelinated_terminal_length=0.0,
            connection_method='translate_axon',
            connect_gap_junction=False,
            gap_conductances=None,
            raise_if_existing=True,
            use_initial_segment=True,
            rng=None):
        """
        Define compartments types that constitute the axon model.

        @post   attributes 'compartment_defs' and 'repeating_comp_sequence'
                have been set.

        Arguments
        ---------

        @param  termination_method : str

                - 'any' to terminate as soon as last compartment extends beyond
                  endpoint of streamline

                - 'unmyelinated' to end with an unmyelinated segment of length
                  <unmyelinated_terminal_length>

                - 'terminal_sequence' to use the axon class' <terminal_comp_sequence>
                  variable for building the final segment

                - 'nodal_extend' to terminate axon with nodal compartment
                   extended beyond streamline endpoint if necessary

                - 'nodal_cutoff' to terminate axon with nodal compartment
                   before streamline endpoint is reached


        @param  tolerance_mm : float

                Tolerance on length of compartments caused by discrepancy
                between streamline arclength and length of compartments spanning
                a streamline node.

        @param  connection_method : str

                Method for connecting the reconstructed acon to the parent
                sections. One of the following:

                'orient_coincident': find streamline end that connects to the
                parent section and start building from this point. Raise
                exception if not coincident.

                'translate_axon_<start/end>': translate starting or endpoint
                of axon to endpoint of parent section.

                'translate_cell_<start/end>': Translate cell so that connection
                point is coincident with start or end of streamline

        @param  connect_gap_junction : bool

                If true, connect axon to cell using gap junction instead of
                electrical connection between compartments. In this case the
                new axonal compartmetns are not added to SectionLists in
                the parent cell object (wrapper for compartments).
        """
        self.without_extracellular = without_extracellular

         # Parse arguments and check preconditions
        if termination_method == 'unmyelinated' and \
           unmyelinated_terminal_length <= 0.0:
            raise ValueError('Must use positive unmyelinated terminal length'
                             ' for termination method "unmyelinated"')
        if termination_method == 'terminal_sequence' and \
           len(self.terminal_comp_sequence) == 0:
            raise ValueError('No terminal compartment sequence defined '
                             'in axon class {}'.format(self.__class__))

        # Save properties for building algorithm
        self.tolerance_mm = tolerance_mm
        self.interp_method = interp_method
        self.termination_method = termination_method
        self.use_initial_segment = use_initial_segment
        self.unmyelinated_terminal_length = unmyelinated_terminal_length
        self.streamline_pts = np.array(streamline_coords)
        self.parent_cell = parent_cell
        self.parent_sec = parent_sec
        self.connect_gap_junction = connect_gap_junction
        self.gap_conductances = gap_conductances

        if rng is None:
            self.rng = np.random
        else:
            self.rng = rng

        # Connect to parent cell
        if parent_sec is not None:
            self._join3d_axon2cell(parent_cell, parent_sec, connection_method,
                               tolerance_mm=1e-3)

        # Save streamline info
        self.num_streamline_pts = len(self.streamline_pts)
        self._set_streamline_length()

        # Set interpolation / walking function
        if self.interp_method == 'arclength':
            self.walk_func = self._walk_arclength
        elif self.interp_method == 'cartesian':
            self.walk_func = self._walk_cartesian_length
        else:
            raise ValueError('Unknown interpolation method: ', interp_method)

        # Data structures for collateral definitions
        self.colt_branch_points = []
        self.colt_target_points = []
        self.colt_step_lengths = []
        self.colt_lvl_numbranch = []
        self.colt_lvl_angles = []
        self.colt_lvl_lengths = []
        self.colt_subtrees_seclists = []


    def add_collateral_definition(
            self,
            branch_point_um,
            target_point_um,
            levels_lengths_um,
            step_length_um,
            levels_num_branches,
            levels_angles_deg):
        """
        Add axon collateral definition at branch point coordinates,
        branching out to the given target point.

        Arguments
        ---------

        @param  branch_point_um : np.array[float] of shape (3,)

        @param  target_pt_um : np.array[float] of shape (3,)

        @param  max_length_um : float

        @param  step_length_um : float

        @param  levels_num_steps : tuple[float]

        @param  levels_num_branches : tuple[int]

        @oaran  levels_angles_deg : float
                Maximum angular deviation for step in direction of target point
        """
        self.colt_branch_points.append(branch_point_um)
        self.colt_target_points.append(target_point_um)
        self.colt_step_lengths.append(step_length_um)
        self.colt_lvl_lengths.append(levels_lengths_um)
        self.colt_lvl_numbranch.append(levels_num_branches)
        self.colt_lvl_angles.append(levels_angles_deg)


    def get_terminal_section(self, terminal_spec):
        """
        Get a terminal section based on string specification

        'main_branch:-1' -> last section of main branch
        'branch_point:2:collateral' -> section of type 'collateral' closest to second branch point
        'coordinate:[3.14,4.56,8.234]:collateral' -> section of type 'collateral' closest to coordinate
        """
        terminal_spec = terminal_spec.split(':')
        structure = terminal_spec[0]
        if structure == 'main_branch':
            # Index into the axon's main branch
            index = int(terminal_spec[1])
            return self.built_sections[structure][index]

        elif structure in ('branch_point', 'target_point'):
            # Section close to one of the branch points
            index = int(terminal_spec[1])
            sec_type = terminal_spec[2]
            measure_point = getattr(self, 'colt_%ss' % structure)[index]

        elif structure == 'coordinate':
            coord_spec = terminal_spec[1]
            sec_type = terminal_spec[2]
            float_pattern = r'([0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
            grp_occurrences = re.findall(float_pattern, coord_spec)
            measure_point = np.array([float(grp[0]) for grp in grp_occurrences])

        else:
            raise ValueError(terminal_spec)

        # find section closest to the branch point
        candidate_secs = self.built_sections[sec_type]
        i_sec, i_node = morph_3d.find_closest_section(measure_point,
                            candidate_secs, measure_from='segment_centers')
        return candidate_secs[i_sec]


    def _build_all_collaterals(self):
        """
        Build all axon collaterals according to their definitions
        added with add_collateral_definition().
        """
        # Get 3D coordinates of all nodal sections
        node_type = self.collateral_comp_sequence[0]
        node_secs = self.built_sections[node_type]
        node_pt3d, node_n3d, sl_nsec = morph_3d.get_segment_centers([node_secs],
                                                    samples_as_rows=True)
        node_pt3d = np.array(node_pt3d)

        # Initialize data structures
        self.colt_subtrees_seclists = [[] for i in range(len(self.colt_branch_points))]

        # Branch out from node closest to each branch point
        node_pt_idx_upper = np.cumsum(node_n3d)
        for i_colt, branch_pt_um in enumerate(self.colt_branch_points):

            # Find nodal section that is closest to branch point
            node_pt_dists = np.linalg.norm(node_pt3d - branch_pt_um, axis=1)
            i_pt3d = np.argmin(node_pt_dists)
            i_sec = next((i for i, hi in enumerate(node_pt_idx_upper) if (i_pt3d < hi)))
            parent_sec = node_secs[i_sec]
            parent_pt0_idx = node_pt_idx_upper[i_sec-1] if i_sec > 0 else 0
            parent_node_idx = i_pt3d - parent_pt0_idx
            parent_seg = nrnutil.seg_at_node(parent_sec, parent_node_idx)
            parent_pt_um = node_pt3d[i_pt3d]

            # Initialize distance measurement
            h.distance(0, parent_seg.x, sec=parent_sec)

            # Build the collateral
            for i in range(self.colt_lvl_numbranch[i_colt][0]):
                self._build_collateral_subtree(i_colt, 0, parent_pt_um, parent_seg)


    def _build_collateral_subtree(self, colt_idx, level, parent_pt_um, parent_seg):
        """
        Build one unbranched section between two collateral branch points
        at subsequent levels of collateral subtree and recurse.
        """

        # Get data defining the collateral
        target_pt_um = self.colt_target_points[colt_idx]
        step_length_um = self.colt_step_lengths[colt_idx]
        angle_deg = self.colt_lvl_angles[colt_idx][level]

        level_length_um = self.colt_lvl_lengths[colt_idx][level]
        max_length_um = sum(self.colt_lvl_lengths[colt_idx])
        num_steps = int(np.ceil(level_length_um / step_length_um))
        max_level = len(self.colt_lvl_lengths[colt_idx])-1

        # Build the Section that will hold all the 3D points
        sec_type = self.collateral_comp_sequence[-1]
        ax_sec = self._make_sec(sec_type)
        sec_attrs = self.compartment_defs[sec_type]
        diam = sec_attrs['morphology']['diam']
        self._set_comp_attributes(ax_sec, sec_attrs)
        ax_sec.connect(parent_seg, 0.0)

        # Save the compartment
        self.built_sections[sec_type].append(ax_sec)
        self.colt_subtrees_seclists[colt_idx].append(ax_sec)

        # Add 3D point to Section by stepping towards target point
        h.pt3dclear(sec=ax_sec)
        x, y, z = parent_pt_um
        h.pt3dadd(x, y, z, diam, sec=ax_sec)
        i_step = 0
        u_step = normvec(target_pt_um - parent_pt_um)
        reached_target = False
        while True:
            i_step += 1

            # Step in the direction of target
            v_target = target_pt_um - parent_pt_um
            if np.linalg.norm(v_target) > step_length_um and not reached_target:
                u_target = normvec(v_target)
            else:
                # If we are about to bump into target point, walk in direction of previous step
                reached_target = True
                u_target = u_step
            next_pt_pre_rot = parent_pt_um + step_length_um * u_target

            # Choose random axis perpendicular to sample axis
            rvec = self.rng.rand(3) - 0.5
            rvec_parallel = np.dot(rvec, u_target) * u_target
            rvec_perpendicular = rvec - rvec_parallel

            # Make rotation matrix for all samples in subtree of current sample
            x = self.rng.rand(1) - 0.5
            angle = np.deg2rad(2 * x * angle_deg)
            A = axangles.axangle2aff(rvec_perpendicular, angle, point=parent_pt_um)
            next_pt_post_rot = np.dot(A, np.append(next_pt_pre_rot, 1.0))[:3]

            # Add 3d information
            x, y, z = next_pt_post_rot
            h.pt3dadd(x, y, z, diam, sec=ax_sec)

            # Update parent point
            u_step = normvec(next_pt_post_rot - parent_pt_um)
            parent_pt_um = next_pt_post_rot


            # Measure distance to collateral endpoint
            collateral_length_um = h.distance(1.0, sec=ax_sec)
            if level == max_level:
                if collateral_length_um >= max_length_um:
                    return
                else:
                    continue
            elif i_step == num_steps:
                break


        # Reached end of unbranched section without meeting length criterion
        for i in range(self.colt_lvl_numbranch[colt_idx][level+1]):
            self._build_collateral_subtree(colt_idx, level+1, parent_pt_um, ax_sec(1.0))


    def build_axon(self):
        """
        Build NEURON axon along a sequence of coordinates.

        Returns
        -------

        @return     sections : list[nrn.Section]
                    List of axonal sections
        """
        # State variables for building algorithm
        self.interp_pts = []        # interpolated points
        self.num_passed = 1         # index of last passed streamline point
                                    # (we already 'passed' starting point)

        self.interp_pts = [self.streamline_pts[0]]        # interpolated points
        self.last_coord = self.streamline_pts[0]
        self.last_tangent = normvec(self.streamline_pts[1] - self.streamline_pts[0])
        self.built_length = 0.0
        self.i_sequence_offset = 0


        # Keep references to newly constructed Sections
        self.built_sections = {sec_type: [] for sec_type in self.compartment_defs.keys()}
        self.built_sections['main_branch'] = sec_ordered = []

        # Estimate number of Sections needed to build axon
        MAX_NUM_COMPARTMENTS = int(1e9)
        est_num_comp = self.estimate_num_sections()
        tot_num_seg = 0
        logger.debug("Estimated number of sections to build axon "
                     " of length {} mm: {:.1f}".format(self.streamline_length, est_num_comp))
        if est_num_comp > MAX_NUM_COMPARTMENTS:
            raise ValueError('Streamline too long (estimated number of '
                'compartments needed is {}'.format(est_num_comp))

        # Keep track of compartments that exceed geometry tolerance
        num_tol_exceeded = 0
        diffs_tol_exceeded = []

        # Build axon progressively by walking along streamline path
        prev_sec = self.parent_sec
        for i_compartment in xrange(MAX_NUM_COMPARTMENTS):

            # Find properties of next Section (compartment type)
            sec_type, sec_attrs = self.get_next_compartment_def(i_compartment)
            sec_L_mm = sec_attrs['morphology']['L'] * 1e-3 # um to mm


            # Create the compartment
            ax_sec = self._make_sec(sec_type)
            self._set_comp_attributes(ax_sec, sec_attrs)

            # Connect the compartment
            if prev_sec is not None:
                if i_compartment == 0 and self.connect_gap_junction:
                    # use mechanism gap.mod, comes with PyNN
                    seg_pre = prev_sec(0.5)
                    seg_post = ax_sec(0.5)
                    gap_pre = h.Gap(seg_pre)
                    gap_post = h.Gap(seg_post)
                    gap_pre.g = self.gap_conductances[0]
                    gap_post.g = self.gap_conductances[1]
                    h.setpointer(seg_pre._ref_v, 'vgap', gap_post)
                    h.setpointer(seg_post._ref_v, 'vgap', gap_pre)
                    # Save reference to gap junctions
                    self.built_sections['gap_junctions'] = (gap_pre, gap_post)
                else:
                    ax_sec.connect(prev_sec(1.0), 0.0)

            prev_sec = ax_sec
            self.built_sections[sec_type].append(ax_sec)
            sec_ordered.append(ax_sec)
            tot_num_seg += ax_sec.nseg

            # Find section endpoint by walking along streamline for sec.L
            num_passed, stop_coord, next_tangent = self.walk_func(sec_L_mm)
            self.interp_pts.append(stop_coord)

            # Check compartment length vs tolerance
            real_length = veclen(stop_coord - self.last_coord)
            if not np.isclose(real_length, sec_L_mm, atol=self.tolerance_mm):
                num_tol_exceeded += 1
                diffs_tol_exceeded.append(abs(real_length - sec_L_mm))
                # logger.warning('exceed length tolerance ({}) '
                #                ' in compartment compartment {} : L = {}'.format(
                #                     tolerance_mm, ax_sec, real_length))

            # Add the 3D start and endpoint
            sec_endpoints = self.last_coord, stop_coord
            h.pt3dclear(sec=ax_sec)
            for coords_mm in sec_endpoints:
                coords_um = coords_mm * 1e3
                x, y, z = coords_um
                h.pt3dadd(x, y, z, sec_attrs['morphology']['diam'], sec=ax_sec)

            # Update state variables
            self.built_length += real_length
            self.num_passed += num_passed
            self.last_coord = stop_coord
            self.last_tangent = next_tangent

            # If terminating axon with nodal compartment: either cutoff or extrapolate
            # remaining_length = self._get_remaining_arclength() # FIXME: fix bug
            remaining_length = self.streamline_length - self.built_length

            if self.termination_method == 'terminal_sequence':
                if self.current_compartment_sequence == 'terminal_comp_sequence' and \
                   (self.i_sequence == len(self.terminal_comp_sequence)-1):
                   break
            elif self.termination_method == 'unmyelinated':
                if remaining_length <= self.unmyelinated_terminal_length:
                    # This check was also done before building previous compartment
                    # and should have built a compartment of exactly the remaining length
                    break
            elif self.termination_method == 'any_extend':
                if self.num_passed >= len(self.streamline_pts):
                    break
            elif self.termination_method == 'any_cutoff':
                next_type, next_attrs = self.get_next_compartment_def(
                                            i_compartment + 1, update_state=False)
                next_length = next_attrs['morphology']['L'] * 1e-3
                if remaining_length - next_length <= 0:
                    break
            elif self.termination_method == 'nodal_cutoff' and sec_type == self.nodal_compartment_type:
                i_seq = (i_compartment - self.i_sequence_offset) % len(self.repeating_comp_sequence)
                next_node_dist = self.next_node_distance(i_seq, measure_from=1)
                if next_node_dist > remaining_length:
                    break
            elif self.termination_method == 'nodal_extend' and sec_type == self.nodal_compartment_type:
                if self.num_passed >= len(self.streamline_pts):
                    break

            # Sanity check: is axon too long?
            if i_compartment >= MAX_NUM_COMPARTMENTS-1:
                raise ValueError("Axon too long.")
            elif i_compartment >= 1.1 * est_num_comp:
                logger.warning("Created {}-th section, more than estimate {}".format(
                                i_compartment, est_num_comp))

        # Connext axon collaterals to main axon
        self._build_all_collaterals()

        # Status report
        logger.debug("Created %i axonal segments (%i sections)",
                     tot_num_seg, len(sec_ordered))

        if num_tol_exceeded > 0:
            logger.debug('exceeded length tolerance in %d compartments '
                            '(worst = %f mm, tolerance = %f mm)', num_tol_exceeded,
                            max(diffs_tol_exceeded), self.tolerance_mm)

        # Add to parent cell
        if (self.parent_cell is not None) and not self.connect_gap_junction:
            # NOTE: References to Python-created sections must be maintained
            # in Python. References are not kept alive by appending them
            # to one of the Hoc template's SectionList.

            # Add to axonal SectionList in order of connection
            append_to = ['all', 'axonal']
            for seclist_name in append_to:
                seclist = getattr(self.parent_cell, seclist_name, None)
                if seclist is not None:
                    for sec in self.built_sections[self.collateral_comp_sequence[-1]]:
                        seclist.append(sec=sec)
                    for sec in sec_ordered:
                        seclist.append(sec=sec)
                    logger.debug("Updated SectionList '%s' of %s", seclist_name,
                                 self.parent_cell)

        # Make SectionList objects with names conforming to MorphologyImporter interface
        all_seclist = h.SectionList()
        for sec in self.built_sections['main_branch']:
            all_seclist.append(sec=sec)
        for sec in self.built_sections[self.collateral_comp_sequence[-1]]:
            all_seclist.append(sec=sec)
        self.built_sections['all'] = all_seclist
        self.built_sections['axonal'] = all_seclist
        self.built_sections['somatic'] = h.SectionList()

        return dotdict(self.built_sections)


    def _make_sec(self, sec_type):
        sec_name = "{:s}[{:d}]".format(sec_type, len(self.built_sections[sec_type]))
        if (self.parent_cell is not None) and not self.connect_gap_junction:
            return h.Section(name=sec_name, cell=self.parent_cell) # see module neuron.__init__
        else:
            return h.Section(name=sec_name)


    def estimate_num_sections(self):
        """
        Estimate number of of Sections needed to build axon along streamline.
        """
        tck_length_mm = np.sum(np.linalg.norm(np.diff(self.streamline_pts, axis=0), axis=1))
        rep_length_um = sum((self.compartment_defs[sec]['morphology']['L']
                                    for sec in self.repeating_comp_sequence))
        return 1e3 * tck_length_mm / rep_length_um * len(self.repeating_comp_sequence)


    def get_initial_diameter(self):
        """
        Get diameter of first axonal section
        """
        if len(self.initial_comp_sequence) > 0:
            comp_type = self.initial_comp_sequence[0]
        else:
            comp_type = self.repeating_comp_sequence[0]
        return self.compartment_defs[comp_type]['morphology']['diam']


    def _set_comp_attributes(self, sec, sec_attrs):
        """
        Create compartment from properties.
        """
        # Insert mechanisms and set their parameters
        for mech_name, mech_params in sec_attrs['mechanisms'].items():
            if mech_name == 'extracellular' and self.without_extracellular:
                continue
            sec.insert(mech_name)
            for pname, pval in mech_params.items():
                if mech_name == 'extracellular':
                    for i in range(2): # default nlayer = 2
                        getattr(sec, pname)[i] = pval
                else:
                    setattr(sec, '{}_{}'.format(pname, mech_name), pval)

        # Set passive parameters
        for pname, pval in sec_attrs['passive'].items():
                setattr(sec, pname, pval)

        # Number of segments (discretization)
        sec.nseg = sec_attrs['morphology'].get('nseg', 1)
        return sec


    def next_node_distance(self, i_sequence, measure_from=1):
        """
        Distance to next node in mm.
        """
        if measure_from==1 and i_sequence==len(self.repeating_comp_sequence)-1:
            return 0.0
        i_start = i_sequence+1 if measure_from==1 else i_sequence
        dist_microns = sum((self.compartment_defs[t]['morphology']['L']
                                for t in self.repeating_comp_sequence[i_start:]))
        return 1e-3 * dist_microns


    def compartment_sequence_length(self, compartment_sequence):
        """
        Get length of a sequence of compartments in mm.

        @param  compartment_sequence : list[str]
                Sequence of compartment types.
        """
        return 1e-3 * sum((self.compartment_defs[sec]['morphology']['L'] for sec in compartment_sequence))


    def _set_streamline_length(self):
        """
        Calculate length of streamline and set array of distances
        from start of streamline to each intermediate point.
        """
        # length = np.sum(np.linalg.norm(np.diff(self.streamline_pts, axis=0), axis=1))
        dvecs = np.diff(self.streamline_pts, axis=0)
        dvecs_norms = np.sqrt(np.sum(dvecs*dvecs, axis=1))
        self.streamline_length = np.sum(dvecs_norms)
        self.streamline_segment_lengths = np.concatenate(([0.0], dvecs_norms))


    def _get_remaining_arclength(self):
        """
        Get remaining length of axon to be built by counting length
        alone streamline from the last compartment endpoint.
        """
        if self.num_passed >= self.num_streamline_pts:
            return 0.0
        dist_stop2next = veclen(self.last_coord - self.streamline_pts[self.num_passed])
        if self.num_passed >= self.num_streamline_pts-1:
            return dist_stop2next
        dist_next2end = np.sum(self.streamline_segment_lengths[self.num_passed+1:])
        return dist_stop2next + dist_next2end


    def _get_segment_vec(self, i):
        """
        Get the vector (b - a) of two consecutive streamline points
        starting at index i.

        For the last point the vector between the previous point is returned.
        """
        if i == len(self.streamline_pts) - 1:
            i = i - 1
        return np.diff(self.streamline_pts[i:i+2], axis=0).reshape((-1,))


    def _walk_arclength(self, dist):
        """
        Walk for given distance along streamline.

        If you walk past a streamline point, the arclength between
        start and stop point will be larger than the cartestian distance.
        Hence the actual compartment length will be shorter than intended.
        This shouldn't be a problem though since action potentials will just
        be slightly less attenuated.

        Returns
        -------

        @return     num_passed : int
                    Number of streamline points passed during walk.

        @return     stop_coord : np.array[float] (3x1)
                    Coordinates of endpoint of walk.

        @return     next_tangent : np.array[float] (3x1)
                    Unit tangent vector at endpoint of walk
        """
        remaining_pts = self.streamline_pts[self.num_passed:]
        last_walk_coord = self.last_coord
        dist_walked = 0.0
        num_passed = 0

        # Keep walking until we've walked for <dist> mm
        while dist_walked < dist:
            # As long as we are not at the end, set waypoint to next streamline coordinate
            if num_passed < len(remaining_pts):
                waypoint_coord = remaining_pts[num_passed]
            else:
                # walking distance extends beyond end of streamline
                print('Extending axon beyond streamline endpoint for {} mm'.format(
                       dist-dist_walked))
                waypoint_coord = last_walk_coord + 1000.0 * self.last_tangent

            # Walk up to next stopping point
            path_vec = waypoint_coord - last_walk_coord
            dist_waypoint = veclen(path_vec) # distance to next streamline point

            # Choose to stop before or after next route point
            if dist_walked + dist_waypoint < dist:
                # end of walk is beyond next streamline point
                last_walk_coord = waypoint_coord
                dist_walked += dist_waypoint
                num_passed += 1
            else:
                # end of walk is before next streamline point
                tangent_vec = path_vec / dist_waypoint
                stop_pt = last_walk_coord + (dist - dist_walked) * tangent_vec
                break

        return num_passed, stop_pt, tangent_vec


    def _walk_cartesian_length(self, dist):
        """
        Walk along streamline until cartesian distance to last stopping point
        is equal to <dist>.

        Returns
        -------

        @return     stop_coord : np.array[float] (3x1)
                    Coordinates of endpoint of walk.

        @return     next_tangent : np.array[float] (3x1)
                    Unit tangent vector at endpoint of walk

        TODO: fix bug where num_passed not updated correctly (case 1 & 2)
        """
        if self.num_passed >= self.num_streamline_pts:
            remaining_pts = []
        else:
            remaining_pts = self.streamline_pts[self.num_passed:]

        start_walk_coord = self.last_coord
        walk_passed = 0 # streamline points passed during this walk

        # Walk along remaining points until next point is farther than walking distance
        i_pre = -1  # index (offset) of next point _closer_ than distance
        i_post = -1 # index (offset) of next point _farther_ than distance
        for i, waypoint in enumerate(remaining_pts):
            dist_waypoint = veclen(waypoint - start_walk_coord)
            if dist_waypoint <= dist:
                i_pre = i
            else:
                i_post = i
                break

        # Identify one of four cases
        colinear = False
        if (i_pre >= 0) and (i_post >= 0):
            # We are moving from one line segment to one of the following
            # line segments.
            colinear = False
            walk_passed = i_pre + 1
            p1 = remaining_pts[i_pre]
            p2 = remaining_pts[i_post]
        elif i_post == 0:
            # Stopping point is on the current line segment
            colinear = True
            walk_passed = 0
            tangent = normvec(self._get_segment_vec(self.num_passed))
        elif (i_pre == -1) and (i_post == -1):
            # We are past the end of the streamline (extending it)
            assert self.num_passed == self.num_streamline_pts
            colinear = True
            walk_passed = 0
            tangent = self.last_tangent
        elif (i_post == -1) and (i_pre == 0):
            # We are on last line segment and stopping point is beyond endpoint
            # of streamline.
            assert self.num_passed == self.num_streamline_pts - 1
            colinear = True
            walk_passed = 1
            tangent = normvec(self._get_segment_vec(self.num_streamline_pts-1))
        elif (i_pre != -1) and (i_post == -1):
            # We are moving from a line segment to beyond the endpoint of
            # the streamline
            colinear = False
            streamline_i_pre = self.num_passed + i_pre + 1
            assert streamline_i_pre == self.num_streamline_pts - 1
            walk_passed = i_pre + 1
            tangent = normvec(self._get_segment_vec(streamline_i_pre))
            p1 = remaining_pts[i_pre]
            p2 = p1 + 100.0 * tangent
        else:
            assert False, "Unknown condition, should not occur."

        # Choose the interpolation method based on the case
        if colinear:
            # Stopping point is colinear with the current line segment
            stop_pt = start_walk_coord + dist * tangent

        else:
            # Stopping point lies on [p1, p2) and is not necessarily colinear
            # with the line segment
            p0 = start_walk_coord

            # Find point on [p1, p2) that yields cartesian distance to p0
            # TODO: check calculation
            u12 = normvec(p2 - p1)
            v01 = p1 - p0
            coefficients = [                # a*x^2 + b*x + c
                np.dot(u12, u12),           # a
                2 * np.dot(u12, v01),       # b
                np.dot(v01, v01) - dist**2, # c
            ]
            roots_all = np.roots(coefficients)
            roots_real = roots_all[np.isreal(roots_all)]
            alpha = np.max(roots_real)
            assert alpha > 0

            # Calculate stop point
            stop_pt = p1 + alpha * u12
            tangent = u12

        return walk_passed, stop_pt, tangent


    def _join3d_axon2cell(self, parent_cell, parent_sec, connection_method,
                      tolerance_mm=1e-3):
        """
        Ensure axon is connected to parent cell in 3D space.

        @see    build_along_streamline.

        @post   self.streamline_translation_mm is the translation vector applied
                to the original streamline points
        """
        # Get connection point on parent cell (assume last 3D point)
        n3d = int(h.n3d(sec=parent_sec))
        parent_coords = np.array([h.x3d(n3d-1, sec=parent_sec),
                                  h.y3d(n3d-1, sec=parent_sec),
                                  h.z3d(n3d-1, sec=parent_sec)]) * 1e-3 # um to mm

        # Connect axon according to method
        if connection_method == 'orient_coincident':
            # Check which end of streamline is coincident with connection
            # point and build in appropriate direction
            if np.allclose(parent_coords, self.streamline_pts[-1], atol=tolerance_mm):
                self.streamline_pts = self.streamline_pts[::-1] # reverse
            elif not np.allclose(parent_coords, self.streamline_pts[0], atol=tolerance_mm):
                raise ValueError("Start or end of streamline must be coincident"
                        "with endpoint of parent section ({})".format(
                            parent_coords))
            self.streamline_translation_mm = np.zeros(3)

        elif connection_method.startswith('translate_axon'):
            # Translate axon start or end to connection point
            if connection_method.endswith('start'):
                streamline_origin = self.streamline_pts[0]
            elif connection_method.endswith('end'):
                streamline_origin = self.streamline_pts[-1]
            elif connection_method.endswith('closest'):
                dist_ax_start = np.linalg.norm(parent_coords - self.streamline_pts[0])
                dist_ax_end = np.linalg.norm(parent_coords - self.streamline_pts[-1])
                if dist_ax_end < dist_ax_start:
                    self.streamline_pts = self.streamline_pts[::-1]
                streamline_origin = self.streamline_pts[0]
            else:
                raise ValueError(connection_method)

            translate_vec = parent_coords - streamline_origin
            self.streamline_pts = self.streamline_pts + translate_vec # broadcasts
            self.streamline_translation_mm = translate_vec
            logger.debug('Axon coordinates were translated by vector {}'.format(translate_vec))

        elif connection_method.startswith('translate_cell'):
            # Translate cell so that connection point is coincident with
            # start or end of streamline
            if connection_method.endswith('start'):
                target_pt = self.streamline_pts[0]
            elif connection_method.endswith('end'):
                target_pt = self.streamline_pts[-1]
            else:
                raise ValueError(connection_method)

            translate_vec = target_pt - parent_coords
            translate_mat = np.eye(4)
            translate_mat[:3,3] = translate_vec
            morph_3d.transform_sections(parent_cell.all, translate_mat)
            self.streamline_translation_mm = np.zeros(3)
            logger.debug('Cell morphology was transformed using matrix {}'.format(translate_mat))


    def get_next_compartment_def(self, i_compartment, update_state=True):
        """
        Get type and attributes of next compartment to be built.

        @return     sec_type, sec_attrs : tuple[str, dict]
        """

        # Find properties of next Section (compartment type)
        self.i_initial_offset = 0
        if self.use_initial_segment and (i_compartment < len(self.initial_comp_sequence)):
            # We are in the initial section of the axon (non-repeating structure)
            current_compartment_sequence = 'initial_comp_sequence'
            i_sequence = i_compartment
            sec_type = self.initial_comp_sequence[i_sequence]
            self.i_initial_offset = len(self.initial_comp_sequence)

        else:
            # We are in the repeating part of the axon
            current_compartment_sequence = 'repeating_comp_sequence'
            i_sequence = (i_compartment - self.i_initial_offset) % len(self.repeating_comp_sequence)
            sec_type = self.repeating_comp_sequence[i_sequence]

        sec_attrs = self.compartment_defs[sec_type]

        remaining_length = self.streamline_length - self.built_length
        next_remaining_length = remaining_length - sec_attrs['morphology']['L'] * 1e-3

        # Check termination criteria that change sequence
        if (self.termination_method == 'unmyelinated') and \
           (self.unmyelinated_terminal_length > 0.0):
            if next_remaining_length <= self.unmyelinated_terminal_length:
                # Next compartment will be an unmyelinated section of length 'remaining_length'
                current_compartment_sequence = 'repeating_comp_sequence'
                sec_type = self.nodal_compartment_type
                sec_attrs = dict(self.compartment_defs[sec_type])
                seg_per_mm = sec_attrs['morphology']['nseg'] / sec_attrs['morphology']['L'] * 1e3
                sec_attrs['morphology']['nseg'] = int(seg_per_mm * remaining_length) + 1
                sec_attrs['morphology']['L'] = remaining_length * 1e-3

        elif self.termination_method == 'terminal_sequence':
            if next_remaining_length <= self.compartment_sequence_length(self.terminal_comp_sequence):
                current_compartment_sequence = 'terminal_comp_sequence'
                i_sequence = (i_compartment - self.i_terminal_offset) % len(self.terminal_comp_sequence)
                sec_type = self.terminal_comp_sequence[i_sequence]
                sec_attrs = self.compartment_defs[sec_type]
            else:
                self.i_terminal_offset = len(self.built_sections['main_branch'])

        if update_state:
            self.current_compartment_sequence = current_compartment_sequence
            self.i_sequence = i_sequence

        return sec_type, sec_attrs


# Make logging functions
# for level in logging.INFO, logging.DEBUG, logging.WARN:
#     def make_func(level=level):
#         # inner function solves late binding issue
#         def log_func(self, *args, **kwargs):
#             if self.logger is None:
#                 return
#             self.logger.log(level, *args, **kwargs)
#         return log_func
#     level_name = logging.getLevelName(level).lower()
#     setattr(AxonBuilder, level_name, make_func(level))
