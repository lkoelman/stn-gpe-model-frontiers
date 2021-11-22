"""
Recursive merging procedures for cylindrical compartments

@author     Lucas Koelman
@date       17-01-2017
"""

# Python modules
import math
PI = math.pi
from copy import deepcopy
import functools

# Our modules
from bgcellmodels.common.treeutils import next_segs
from bgcellmodels.common.nrnutil import seg_xmin, seg_xmax, getsecref
from bgcellmodels.common.electrotonic import seg_lambda, calc_lambda_AC
from . import redutils

from numpy import interp

################################################################################
# Rewrite from formulas
################################################################################

class MergedProps(object):
    """
    Information about merged properties

    NOTE:   This is the 'Bunch' recipe from the python cookbook.
            An alternative would be `myobj = type('Bunch', (object,), {})()`
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class Cylinder(object):
    """
    Cylindrical compartment, the fundamental element of the cable
    representation of a neurite.
    """

    def __init__(
            self,
            nrn_seg=None, use_seg_x=None, x_bound=None,
            diam=None, L=None, Ra=None
        ):
        """
        Create new cylindrical compartment.

        @param  nrn_seg : nrn.Segment
                A NEURON segment to copy properties from. Provide either this
                or explicit properties.

        @param  use_seg_x : (optional) str
                If not None: split given segment depending on value:
                "x_min" use nrn_seg.x as its left boundary
                "x_max" use nrn_seg.x as its right boundary

        @param  x_bound: (optional) float
                
                A second x-value that should lie within the segment.
                The value of 'use_seg_x' determines how 'x_bound' is used:
                
                - 'x_min': segment.x is left boundary and x_bound is right boundary
                - 'x_max': segment.x is right boundary and x_bound is left boundary
                - other: x_bound is not used

        @param  diam : float
                Cylinder diameter in units of micron

        @param  L : float
                Cylinder length in units of micron

        @param  Ra : float
                Cytoplasmic resistivity in units of [Ohm * cm]
        """
        if nrn_seg is not None:

            self.diam = nrn_seg.diam
            self.Ra = nrn_seg.sec.Ra
            
            if use_seg_x is None:
                self.L = nrn_seg.sec.L / nrn_seg.sec.nseg
                self.x_lo = seg_xmin(nrn_seg, side='inside')
                self.x_hi = seg_xmax(nrn_seg, side='inside')
            elif use_seg_x == "x_min":
                self.x_lo = nrn_seg.x
                self.x_hi = x_bound if x_bound is not None else seg_xmax(nrn_seg, side='inside')
                self.L = nrn_seg.sec.L * (self.x_hi - self.x_lo)
            elif use_seg_x == "x_max":
                self.x_hi = nrn_seg.x
                self.x_lo = x_bound if x_bound is not None else seg_xmin(nrn_seg, side='inside')
                self.L = nrn_seg.sec.L * (self.x_hi - self.x_lo)
            else:
                raise ValueError("Argument 'use_seg_x': invalid value {}".format(use_seg_x))
        else:
            self.diam = diam
            self.L = L
            self.Ra= Ra

            if any([nrn_seg, use_seg_x]):
                raise ValueError("Provide either a NEURON segment or cylinder properties "
                                 "but not both.")


    @property
    def area(self):
        return PI * self.diam * self.L


    @property
    def Ri(self):
        """
        Absolute axial resistance in units of Ohm
        """
        # units = ([Ohm*cm] * [um]) / ([um]^2)
        # units = ([Ohm] * [m]*1e-2 * [m]*1e-6) / ([m]*1e-6 * [m]*1e-6)
        # units = [Ohm] * [m] * 1e4
        # units = [Ohm] * [cm] * 1e2
        unit_factor = 1e2
        return (4.0 * self.Ra * self.L) / (PI * self.diam**2 * unit_factor)


    def save_merged_properties(self, ref, glist):
        """
        Use portion of segment between sec(self.x_lo) and sec(self.x_hi)
        to save area-specific properties.

        @post   self.orig_props is SegProps object with attributes 
                Ra: list, L: list, diam: list, cm: list,
                area: float, cmtot: float, gtot: dict(str,float)
        """
        # Get original segment
        sec = ref.sec
        xmid = (self.x_hi + self.x_lo) / 2.0
        seg = sec(xmid)

        pmerge = MergedProps()
        pmerge.Ra = [sec.Ra]
        pmerge.L = [self.L]
        pmerge.diam = [self.diam]
        pmerge.cm = [seg.cm]

        # Calculate area, capacitance, conductances
        pmerge.area = self.area
        pmerge.cmtot = seg.cm * self.area
        pmerge.gtot = dict((gname, 0.0) for gname in glist)
        for gname in glist:
            try:
                gval = getattr(seg, gname, 0.0)
            except NameError: # NEURON error if mechanism not inserted in section
                gval = 0.0
            pmerge.gtot[gname] = gval * pmerge.area

        # Information about absorbed segments
        pmerge.merged_sec_gids = set([ref.gid])
        pmerge.merged_sec_names = set([sec.name()])
        pmerge.merged_region_labels = set([ref.region_label])

        self.orig_props = pmerge


    @staticmethod
    def update_merged_properties(pmerge, pother):
        """
        Update merged properties

        @param  cyl : Cylinder
                Cylinder object of which orig_props should be merged
        """

        pmerge.Ra.extend(pother.Ra)
        pmerge.L.extend(pother.L)
        pmerge.diam.extend(pother.diam)
        pmerge.cm.extend(pother.cm)

        # Calculate area, capacitance, conductances
        pmerge.area += pother.area
        pmerge.cmtot += pother.cmtot
        for gname, gsum in pother.gtot.items():
            pmerge.gtot[gname] += gsum

        # Information about absorbed segments
        pmerge.merged_sec_gids.update(pother.merged_sec_gids)
        pmerge.merged_region_labels.update(pother.merged_region_labels)


    @staticmethod
    def combine_merged_properties(all_pmerge):
        """
        Combine MergedProps objects.
        """
        if isinstance(all_pmerge[0], Cylinder):
            pinit = deepcopy(all_pmerge[0].orig_props)
            prest = (cyl.orig_props for cyl in all_pmerge[1:])
        elif isinstance(all_pmerge[0], MergedProps):
            pinit = deepcopy(all_pmerge[0])
            prest = all_pmerge[1:]
        else:
            raise ValueError()

        for pmerge in prest:
            Cylinder.update_merged_properties(pinit, pmerge)

        return pinit

    

def merge_cylinders_sequential(cyls):
    """
    Merge sequential cylinders sequentially using Bush & Sejnowski method,
    i.e. preserving axial resistance (absolute & per unit length).
    
    @param      cyls: list(Cylinder)
                Cilinders to merge sequentially

    @return     Cylinder
    """
    # Trivial case: only a single cylinder
    if len(cyls) == 1:
        return cyls[0]

    # Intermediate calculations
    Ri_seq = sum((cyl.Ri for cyl in cyls))

    # Cylinder geometrical & electrical properties
    L_seq = sum((cyl.L for cyl in cyls))
    Ra_seq = sum((cyl.Ra for cyl in cyls)) / len(cyls)
    diam_seq = math.sqrt(Ra_seq*L_seq*4./PI/Ri_seq/100.) # ensures that Ri_seq will be the total absolute axial resistance of the equivalent cylinder

    merged_cyl = Cylinder(diam=diam_seq, L=L_seq, Ra=Ra_seq)

    # Save properties of original NEURON segments
    merged_cyl.orig_props = Cylinder.combine_merged_properties(cyls)
    return merged_cyl


def merge_cylinders_parallel(cyls):
    """
    Merge parallel cylinders, using Bush & Sejnowski method,
    i.e. preserving axial resistance (absolute & per unit length).

    @param      cyls: list(Cylinder)
                Cilinders to merge

    @return     Cylinder
    """
    # Trivial case: only a single cylinder
    if len(cyls) == 1:
        return cyls[0]

    # Cylinder geometrical & electrical properties
    L_br = sum((cyl.area*cyl.L for cyl in cyls)) / sum((cyl.area for cyl in cyls))
    Ra_br = sum((cyl.Ra for cyl in cyls)) / len(cyls)
    diam_br = math.sqrt(sum((cyl.diam**2 for cyl in cyls))) # ensures that absolute axial resistance Ri (also per unit length, Ri/L) is preserved

    merged_cyl = Cylinder(diam=diam_br, L=L_br, Ra=Ra_br)

    # Save properties of original NEURON segments
    merged_cyl.orig_props = Cylinder.combine_merged_properties(cyls)
    return merged_cyl


def merge_until_distance(
        start_seg,
        stop_dist,
        distance_func,
        all_refs,
        gbar_names,
        lookahead=False):
    """
    Keep merging cylindrical compartments in subtree starting at start_seg until
    distance_func(end_seg) is equal to the desired stopping distance.

    @param  start_seg : nrn.Segment
            Segment to star merging procedure with x-value equal to 
            the last splitting point.
    
    @param  distance_func : callable(nrn.Segment) -> float
            Function that measures distance of a NEURON segment with an
            associated x-value, e.g. from the root of the tree

    @param  stop_dist : float
            The distance at which merging should be stopped, measured using
            distance_func

    @param  lookahead : bool
            If true, this is a look-ahead call, and the visited Sections will
            not be marked as visited or absorbed

    @return tuple(Cylinder, list(nrn.Segment))
            Equivalent cylinder until stopping distance, and the list of
            stopping points
    """
    
    start_ref = getsecref(start_seg.sec, all_refs)
    beg_seg_x = start_seg
    end_seg_x = start_seg.sec(seg_xmax(start_seg, side='boundary'))

    # Keep absorbing segments until the 1-end boundary is larger than distance
    X_beg = distance_func(beg_seg_x)
    X_end = distance_func(end_seg_x)
    if stop_dist < X_beg:
        raise ValueError("Start segment is already farther than stopping distance")

    # Scenario 1: far boundary of starting segment is farther that stop_dist
    # CONSTRAINT: can never return a segment with distance_func(seg) > stop_dist
    if stop_dist <= X_end:
        # Interpolate distance values at start and end of segment
        x_interp = interp(
                    stop_dist, 
                    [X_beg, X_end],
                    [beg_seg_x.x, end_seg_x.x])
        
        # Cut segment: Cylinder runs from segment.x to stopping distance
        start_cyl = Cylinder(nrn_seg=beg_seg_x,
                             use_seg_x='x_min',
                             x_bound=x_interp)
        start_cyl.save_merged_properties(start_ref, gbar_names)
        if not lookahead:
            start_ref.visited = True

        term_seg = beg_seg_x.sec(x_interp)
        return start_cyl, [term_seg]
    
    else: # Scenario 2: stopping distance is beyond far boundary of segment

        # get cylinder from segment.x to right boundary of segment
        start_cyl = Cylinder(nrn_seg=beg_seg_x, use_seg_x="x_min")
        start_cyl.save_merged_properties(start_ref, gbar_names)
        if not lookahead:
            start_ref.visited = True
            start_ref.absorbed = True # fully absorbed when walked past end

        # Scenario 2.1: no child segments
        child_segments = next_segs(beg_seg_x, x_loc="min")
        if not any(child_segments):
            return start_cyl, []

        # Scenario 2.2: child segment can be merged
        # Get child segment(s), absorb their cylinders into starting cylinder
        parallel_cyls = []
        parallel_x_stop = []
        for child_seg in child_segments:
            child_eq_cyl, child_x_term = merge_until_distance(
                                            child_seg, stop_dist, distance_func,
                                            all_refs, gbar_names, lookahead)
            parallel_cyls.append(child_eq_cyl)
            parallel_x_stop.extend(child_x_term)

        # Merge parallel child cylinders (equivalents)
        next_seq_cyl = merge_cylinders_parallel(parallel_cyls)
        
        # Merge children equivalent sequentially into own cylinder
        sequential_cyls = [start_cyl, next_seq_cyl]
        subtree_eq_cyl = merge_cylinders_sequential(sequential_cyls)

        return subtree_eq_cyl, parallel_x_stop


def next_splitpoints_seg_dlambda(current_seg, lambda_start=None, fraction=0.1):
    """
    Ascend tree and return first segment along each path where lambda
    differs by more than given fraction from starting value.

    @param  start_segs: list(nrn.Segment)
            segments with x-values where last split occurred
    """
    # Merge up to first segment where lambda differs more than 10% from start
    # i.e. in each sibling: ascend until you find segment with d_lambda > 10 or branch point
    if lambda_start is None:
        lambda_start = seg_lambda(current_seg, None, 100.0)

    # Check current node
    lambda_current = seg_lambda(current_seg, None, 100.0)
    if not (1.0-fraction <= lambda_current/lambda_start <= 1.0+fraction):
        # lambda changed more than theshold: split before segment
        return [current_seg]

    # Check child nodes
    child_segments = next_segs(current_seg, x_loc="min")
    child_splitpoints = []
    for child_seg in child_segments:
        child_splitpoints.extend(next_splitpoints_seg_dlambda(
                                    child_seg,
                                    lambda_start=lambda_start,
                                    fraction=fraction))

    return child_splitpoints


def next_splitpoints_fixed_distance(cur_seg, distance_func=None, dX=None, start_X=None):
    """
    Ascend tree and return first segment where distance_func(segment) is equal
    to distance_func(start_segment) + dX.

    @note       Algorithm is basically the same as merge_until_distance() except
                this function is purely for the purpose of looking ahead.
    """
    if start_X is None:
        start_X = distance_func(cur_seg)
    stop_X = start_X + dX

    # Keep absorbing segments until the 1-end boundary is larger than distance
    far_seg = cur_seg.sec(seg_xmax(cur_seg, side='boundary'))
    X_beg = distance_func(cur_seg)
    X_end = distance_func(far_seg)

    # Scenario 0: prox boundary is already beyond stopping distance
    if stop_X < X_beg:
        raise ValueError("Start segment is already farther than stopping distance")

    # Scenario 1: far boundary of starting segment is farther that stop_dist
    if stop_X <= X_end:
        # Interpolate distance values at start and end of segment
        x_interp = interp(
                    stop_X, 
                    [X_beg, X_end],
                    [cur_seg.x, seg_xmax(cur_seg, side='boundary')])
        

        term_seg = cur_seg.sec(x_interp)
        return [term_seg]
    
    else: # stop_X >= X_end : Scenario 2: stopping distance is beyond far boundary of segment
        child_segments = next_segs(cur_seg, x_loc="min") # get left bound of child segments

        # Scenario 2.1: no child segments
        if not any(child_segments):
            if cur_seg.x == 1.0: # this was the end of this branch
                return []
            else:
                return [cur_seg.sec(1.0)] # last point returned on this branch

        # Scenario 2.2: child segment present
        child_splitpoints = []
        for child_seg in child_segments:
            child_splitpoints.extend(next_splitpoints_fixed_distance(
                                        child_seg, distance_func, dX, start_X))

        return child_splitpoints


def approximate_lambda_AC(cylinder, freq):
    """
    Approximate length constant for given cylinder,
    based on surface area scaling of the membrane capacitance.

    @param      cylinder : Cylinder
                Cylinder object

    @return     float
                length constant 'lambda'
    """
    # Calculate equivalent membrane capacitance
    cm_weighted_avg = cylinder.orig_props.cmtot / cylinder.area
    area_scale_factor = cylinder.orig_props.area / cylinder.area
    cm_scaled = cm_weighted_avg * area_scale_factor
    
    return calc_lambda_AC(freq, cylinder.diam, cylinder.Ra, cm_scaled)


class MergingWalk(object):

    def __init__(
            self,
            start_node,
            allsecrefs,
            **kwargs
        ):
        """
        @param  start_node : nrn.Segment
                Segment where children will be merged.

        @param  merge_method: str
                TODO: currently the merge method is always Bush & Sejnowski method.
                Adapt this to allow using Stratford & Larkman method. Also only
                a restricted subset of 'split_criterion' should be allowed for 
                any of those two methods.

        @param  split_criterion : str
                
                Criterion for splitting cylinders. One of the following:
                
                    - 'micron_distance': split cylinders after walking 
                      'split_dX' microns along each branch.
                    
                    - 'eq_electrotonic_distance': walk along branch until the 
                      approximate change in L/lambda in the equivalent cylinder
                      satisfies d(L/lambda) == split_dX. Requires:

                        + additional argument: 'lookahead_units', determining the
                          distance function and units used to interpret 'lookahead_dX'.
                          Either 'micron' or 'lambda'
                        
                        + additional argument: 'lookahead_dX', which is the lookahead
                          distance for approximating lambda
                    
                    - 'next_seg_dlambda': split cylinder when change in length constant is 
                      larger than fixed fraction of starting value. Requires additional
                      arguments: 'lambda_fraction'

                    - (NOT IMPLEMENTED) "next_seg_ddiam": merge until diam has changed 
                      by given fraction of starting value

                    - (NOT IMPLEMENTED) "branch_electrotonic_distance": walk along
                      each branch for fixed step size measures in local branch
                      d(L/lambda)

        @param  (optional, keyword) split_dX: float
                The step size in units determined by argument 'split_criterion'.
                This argument is only required if the 'split_criterion' argument is
                'micron_distance' or 'electrotonic_distance'

        @param  gbar_names : list(str)
                list of mechanism conductances
        
        @param  kwargs : keyword argumments (**dict)
                Additional splitting parameters passed to next_splitpoints() function
        """
        self.start_node = start_node
        self.allsecrefs = allsecrefs

        # Set distance functions
        self.micron_distance_func = functools.partial(redutils.seg_path_L,
                                                      endpoint='segment_x')

        self.elec_distance_func = functools.partial(
                                        redutils.seg_path_L_elec,
                                        f=kwargs['f_lambda'],
                                        gleak=kwargs['gleak_name'],
                                        endpoint='segment_x')

        self.walking_distance_func = self.micron_distance_func
        self.split_criterion = kwargs.pop('split_criterion')
        
        self.gbar_names = kwargs.pop('gbar_names')
        self.split_params = kwargs


    def merge_cylinders_subtree(self):
        """
        Merge all cylinders in subtree in one shot.

        NOTE: the reason that we write a wrapper around merge_until_distance() is
        that in the top-level loop, we don't want to do sequential merging since
        this does not preserve diameter tapering.
        """

        # Initial cable splitting points
        current_splitpoints = [sec(0.0) for sec in self.start_node.sec.children()]
        # Distance (from subtree root) of starting point
        target_dist_X = self.walking_distance_func(self.start_node.sec(1.0))
        sequential_cyls = []

        # Continue as long as there is any branch left where we haven't marched until terminal segment
        while len(current_splitpoints) > 0:
        
            # Look-ahead: find distance of next cable splitting points on each branch.
            # Calculate actual walking distance based on the look-ahead.
            target_dist_X = self.get_walking_distance_lookahead(
                                    current_splitpoints,
                                    target_dist_X)

            # On each sibling branch: march up tree until split distance,
            # and return the equivalent cylinder until that point, as well as 
            # the new splitting points (end of merged cylinders)
            parallel_splitpoints = []
            parallel_cyls = []

            for split_seg in current_splitpoints:
                # TODO: make sure it works with above functions. Specifically:
                #       - recompute target_dist_X based on actual splitpoints?
                #       - end of segment?
                eq_cyl, end_splitpoints = merge_until_distance(
                                split_seg, target_dist_X, 
                                self.walking_distance_func,
                                self.allsecrefs,
                                self.gbar_names)

                if len(end_splitpoints) == 0 and eq_cyl.area < 1e-6 and eq_cyl.L < 1e-3:
                    # walked to end: assert len(split_seg.sec.children()) == 0
                    # do not add tiny leaf cylinders
                    continue
                
                parallel_cyls.append(eq_cyl) # always one cylinder
                parallel_splitpoints.extend(end_splitpoints) # split points only if end not reached


            # Merge the equivalent cylinders in parallel
            if len(parallel_cyls) > 0:
                seq_cyl = merge_cylinders_parallel(parallel_cyls)
                sequential_cyls.append(seq_cyl)
            
            # Update cable splitting points (new starting points)
            current_splitpoints = parallel_splitpoints

        # TODO: postprocess by merging sequential cylinders?
        return sequential_cyls


    def get_walking_distance_lookahead(
            self,
            current_splitpoints,
            last_walk_distance,
        ):
        """
        Look ahead and determine walking distance from last stop to next
        splitting points.

        @param      last_walk_distance : float
                    Walking distance to last stop, measured in units of 
                    self.walking_distance_func

        @return     next_walk_distance : float
                    Walking distance to next stop, measured in units of 
                    self.walking_distance_func
        """
        
        if self.split_criterion == 'micron_distance':
            # No lookahead: step with fixed step size
            return last_walk_distance + self.split_params['split_dX']
        
        elif self.split_criterion == 'eq_electrotonic_distance':

            lookahead_dX = self.split_params['lookahead_dX']
            lookahead_units = self.split_params['lookahead_units']
            lookahead_cyls = []

            # calculate the equivalent cylinder up to 'lookahead_dX'
            for split_seg in current_splitpoints:

                if lookahead_units == 'micron':
                    target_X = last_walk_distance + lookahead_dX
                elif lookahead_units == 'lambda':
                    target_X = self.elec_distance_func(split_seg) + lookahead_dX
                
                branch_eqcyl, _ = merge_until_distance(
                                    split_seg, target_X, self.elec_distance_func,
                                    self.allsecrefs, self.gbar_names, lookahead=True)
                lookahead_cyls.append(branch_eqcyl)

            # Estimate lambda in next equivalent cylinder
            seq_cyl = merge_cylinders_parallel(lookahead_cyls)
            seq_lambda_est = approximate_lambda_AC(seq_cyl, 100.0)

            # split_dX is in units L/lambda -> convert to estimated dL
            return last_walk_distance + self.split_params['split_dX'] * seq_lambda_est

        elif self.split_criterion == 'next_seg_dlambda':

            # Find first point along each branch where change > tolerance occurs
            candidate_x_split = []
            for split_seg in current_splitpoints:
                branch_x_split = next_splitpoints_seg_dlambda(
                                    split_seg,
                                    fraction=self.split_params['split_dX'])
                candidate_x_split.extend(branch_x_split)

            # Find closest point where change > tolerance occurs
            candidate_dists = [self.walking_distance_func(x_seg) for x_seg in candidate_x_split]
            return min(candidate_dists)

        else:
            raise ValueError(self.split_criterion)

