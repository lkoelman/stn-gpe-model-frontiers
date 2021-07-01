"""
Tools for working with extracellular layers in NEURON.

@author     Lucas Koelman
@date       26-04-2019
"""

from __future__ import division
import math
import logging

from neuron import h
import numpy as np

from bgcellmodels.common.treeutils import prev_seg, next_segs

PI = math.pi
sqrt = math.sqrt
logger = logging.getLogger('emfield.xtra_utils')


def set_transfer_impedances(seclist, impedance_lookup_func):
    """
    Set transfer resistances of 'xtra' mechanism using impedance
    lookup function for segment coordinates stored in 'xtra'.

    @param  seclist : neuron.SectionList
            Sections where impedances will be set

    @param  impedance_lookup_function : callable(x, y, x)
            Function returning transfer impedance for with spatial coordinates.
    """
    for sec in seclist:
        for seg in sec:
            x, y, z = (seg.x_xtra, seg.y_xtra, seg.z_xtra)
            R_transfer = impedance_lookup_func(x, y, z)
            seg.rx_xtra = R_transfer


def set_transfer_impedances_nearest(seclist, Z_coords, Z_values,
                                    max_dist, warn_dist, min_electrode_dist,
                                    electrode_coords, Z_intersect=0.0):
    """
    Set transfer impedances using nearest neighbor interpolation, or using
    matching coordinates (set max_dist=eps).

    Uses KD-tree from scipy.spatial for fast lookups.

    @param  Z_coords : array_like of shape (N x 3)
            Coordinates of transfer impedances in Z_values (micron)

    @param  Z_values : arrayLike of length N
            Transfer impedances at coordinates in Z_coords
    """
    import scipy.spatial

    # Construct KD-tree using coordinates of transfer impedance values
    tree = scipy.spatial.KDTree(Z_coords)

    for sec in seclist:
        if not h.ismembrane('xtra', sec=sec):
            logger.debug("Skipping section '%s', no mechanism 'xtra' found.", sec.name())
            continue

        # Get coordinates of compartment centers
        num_samples = int(h.n3d(sec=sec))
        nseg = sec.nseg

        # Get 3D sample points for section
        xx = h.Vector([h.x3d(i, sec=sec) for i in xrange(num_samples)])
        yy = h.Vector([h.y3d(i, sec=sec) for i in xrange(num_samples)])
        zz = h.Vector([h.z3d(i, sec=sec) for i in xrange(num_samples)])

        # Length in micron from start of section to sample i
        pt_locs = h.Vector([h.arc3d(i, sec=sec) for i in xrange(num_samples)])
        L = pt_locs.x[num_samples-1]

        # Normalized location of 3D sample points (0-1)
        pt_locs.div(L)

        # Normalized locations of nodes (0-1)
        node_locs = h.Vector(nseg + 2)
        node_locs.indgen(1.0 / nseg)
        node_locs.sub(1.0 / (2 * nseg))
        node_locs.x[0] = 0.0
        node_locs.x[nseg+1] = 1.0

        # Now calculate 3D locations of nodes (segment centers + 0 + 1)
        # by interpolating 3D locations of samples
        node_xlocs = h.Vector(nseg+2)
        node_ylocs = h.Vector(nseg+2)
        node_zlocs = h.Vector(nseg+2)
        node_xlocs.interpolate(node_locs, pt_locs, xx)
        node_ylocs.interpolate(node_locs, pt_locs, yy)
        node_zlocs.interpolate(node_locs, pt_locs, zz)
        node_coords = np.array(zip(node_xlocs, node_ylocs, node_zlocs))

        # Query KD tree
        nn_dists, nn_idx = tree.query(node_coords, k=1, distance_upper_bound=max_dist)
        for i, c in enumerate(node_locs):
            node_xyz = node_coords[i]
            if np.linalg.norm(node_xyz - electrode_coords) <= min_electrode_dist:
                # Node is too close to electride
                Z_node = Z_intersect
            elif nn_dists[i] > max_dist:
                # Node is too far from nearest neighbor
                raise ValueError("No transfer impedance value found within distance "
                                 "{} of point {}".format(max_dist, node_coords[i]))
            elif nn_dists[i] > warn_dist:
                # Node is too far from nearest neighbor
                Z_node = Z_values[nn_idx[i]]
                logger.debug("Transfer impedance at point {} exceeds warning "
                             "distance of {} um".format(node_coords[i], warn_dist))
            else:
                # No issues, nearest neighbor and electrode distance OK
                Z_node = Z_values[nn_idx[i]]

            # Assign transfer impedance
            sec(c).rx_xtra = Z_node


def set_transfer_impedances_interp(seclist, Z_coords, Z_values,
                                   min_electrode_dist, electrode_coords,
                                   method='linear', Z_intersect=0.0):
    """
    Set transfer impedances using linear or cubic interpolation.

    Similar to function scipy.interpolate.griddata, except the interpolator
    object is saved between queries.

    @param  method : str
            Interpolation method: 'linear' or 'cubic'

    @param  Z_coords : array_like of shape (N x 3)
            Coordinates of transfer impedances in Z_values (micron)

    @param  Z_values : arrayLike of length N
            Transfer impedances at coordinates in Z_coords
    """
    if method == 'linear':
        from scipy.interpolate import LinearNDInterpolator
        interp = LinearNDInterpolator(Z_coords, Z_values)
    elif method == 'cubic':
        from scipy.interpolate import CloughTocher2DInterpolator
        interp = CloughTocher2DInterpolator(Z_coords, Z_values)
    else:
        raise ValueError(method)


    for sec in seclist:
        if not h.ismembrane('xtra', sec=sec):
            logger.debug("Skipping section '%s', no mechanism 'xtra' found.", sec.name())
            continue

        # Get coordinates of compartment centers
        num_samples = int(h.n3d(sec=sec))
        nseg = sec.nseg

        # Get 3D sample points for section
        xx = h.Vector([h.x3d(i, sec=sec) for i in xrange(num_samples)])
        yy = h.Vector([h.y3d(i, sec=sec) for i in xrange(num_samples)])
        zz = h.Vector([h.z3d(i, sec=sec) for i in xrange(num_samples)])

        # Length in micron from start of section to sample i
        pt_locs = h.Vector([h.arc3d(i, sec=sec) for i in xrange(num_samples)])
        L = pt_locs.x[num_samples-1]

        # Normalized location of 3D sample points (0-1)
        pt_locs.div(L)

        # Normalized locations of nodes (0-1)
        node_locs = h.Vector(nseg + 2)
        node_locs.indgen(1.0 / nseg)
        node_locs.sub(1.0 / (2 * nseg))
        node_locs.x[0] = 0.0
        node_locs.x[nseg+1] = 1.0

        # Now calculate 3D locations of nodes (segment centers + 0 + 1)
        # by interpolating 3D locations of samples
        node_xlocs = h.Vector(nseg+2)
        node_ylocs = h.Vector(nseg+2)
        node_zlocs = h.Vector(nseg+2)
        node_xlocs.interpolate(node_locs, pt_locs, xx)
        node_ylocs.interpolate(node_locs, pt_locs, yy)
        node_zlocs.interpolate(node_locs, pt_locs, zz)
        node_coords = np.array(zip(node_xlocs, node_ylocs, node_zlocs))

        # Interpolate node coordinates
        z_vals_nodes = interp(node_coords)

        # Check electrode distance
        for i, c in enumerate(node_locs):
            node_xyz = node_coords[i]
            if np.linalg.norm(node_xyz - electrode_coords) <= min_electrode_dist:
                # Node is too close to electrode
                Z_node = Z_intersect
            else:
                # No issues, nearest neighbor and electrode distance OK
                Z_node = z_vals_nodes[i]

            # Assign transfer impedance
            sec(c).rx_xtra = Z_node


def transfer_resistance_pointsource(seg, seg_coords, source_coords, rho):
    """
    Analytical transfer resistance between electrical point source
    and target coordinate in isotropic non-dispersive extraceullular medium.

    @param  rho : float
            Resistivity of extracullular medium (Ohm * cm).

    @return Z : float
            Transfer impedance (MOhm)
    """
    x1, y1, z1 = seg_coords
    x2, y2, z2 = source_coords
    dist = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    r_min = seg.diam / 2.0
    if dist < r_min:
        dist = r_min

    # 0.01 converts rho's cm to um and ohm to megohm
    return (rho / 4 / PI) * (1 / dist) * 0.01


def get_rattay_activating_function(
        icell,
        seclist_names,
        warn_high_act=False,
        raise_on_warn=False):
    """
    Get activation function values for each compartment (segment),
    grouped by default morphological section lists.
    
    @param  warn_high_act : bool or dict[str, float]
            If true, use default threshold values for determining if activating
            function value is too high, and give warning. If dict is given,
            use custom threshold values and give warning.
    """
    # Parse arguments
    if warn_high_act == True:
        act_warn_thresh = {
            'somatic': 1e6,
            'basal': .5e6,
            'axonal': 4e6,
        }
    elif isinstance(warn_high_act, dict):
        act_warn_thresh = warn_high_act
    else:
        warn_high_act = False

    # Preconditions for algorithm
    stim_amp_mA = 1.0
    h.finitialize()             # assigns seg.node_index()
    test_sec = list(icell.all)[0]
    root_sec = h.SectionRef(sec=test_sec).root
    h.distance(0, 0.5, sec=root_sec)

    # First collect required data for each compartment
    nodes_V_ext = {}            # mV
    nodes_R_mid = {}            # Ohm
    nodes_C_mid = {}            # mF
    all_sec = h.SectionList()
    all_sec.wholetree(sec=root_sec)
    for sec in all_sec:
        if not h.ismembrane('xtra', sec=sec):
            continue
        for seg in sec:
            seg_id = seg.node_index()
            # Extracellular voltage due to electrode
            nodes_V_ext[seg_id] = stim_amp_mA * seg.rx_xtra * seg.scale_stim_xtra # to mV
            # Resistance to parent segment
            nodes_R_mid[seg_id] = seg.ri() * 1e6 # MOhm to Ohm
            # Membrane capacitance
            nodes_C_mid[seg_id] = seg.cm * seg.area() * 1e-11 # uF/cm2 * um2 * cm2/um2 * mF/uF


    def activating_function(seg_mid, parent_seg, child_segs):
        """
        Value of the activating function at segment (compartent), based on
        stimulus-induced current flow from neighboring compartments

        Units for activation function are:
        mV / Ohm / (mA * s / V) = V/s = mV / ms

        Components of the activating function:
        - Ri (axial resistance, Ohm, seg.ri())
            - Ri = 4*Ra*L/(pi*d^2)
        """
        # Activation function according to eq. 4
        seg_id = seg_mid.node_index()
        V_mid = nodes_V_ext[seg_id]
        R_mid = nodes_R_mid[seg_id]
        fn = 0.0
        if parent_seg is not None:
            # use seg_mid.ri() between itself and parent
            V_prev = nodes_V_ext[parent_seg.node_index()]
            fn += (V_prev - V_mid) / R_mid
        for seg in child_segs:
            # use child_seg.ri() between child and parent
            next_id = seg.node_index()
            V_next = nodes_V_ext[next_id]
            R_next = nodes_R_mid[next_id]
            fn += (V_next - V_mid) / R_next
        fn /= nodes_C_mid[seg_id]
        return fn

    act_values = {}             # seclist_name -> list[float]
    dist_values = {}            # seclist_name -> list[float]

    for sl_name in seclist_names:

        seclist = getattr(icell, sl_name, None)
        if seclist is None:
            continue

        act_values[sl_name] = []
        dist_values[sl_name] = []

        for sec in seclist:
            for seg in sec:
                # FIXME: mid or parent or child can have no mechanism 'xtra'

                # Find all compartments physically connected to this one
                # NOTE: do not select the zero-area 0-end and 1-end segments
                parent_seg = prev_seg(seg, x_loc='mid')
                child_segs = next_segs(seg, x_loc='mid')

                # At least two connected segments must exist for curvature
                # (second derivative) to be defined
                if ((parent_seg is not None and len(child_segs) > 0)
                        or len(child_segs) > 1):

                    act = activating_function(seg, parent_seg, child_segs)
                    dist = h.distance(seg.x, sec=sec)

                    act_values[sl_name].append(act)
                    dist_values[sl_name].append(dist)

                    if warn_high_act != False and act > act_warn_thresh[sl_name]:
                        logger.warning("ACT_FUN_HIGH: {}".format(seg))
                        if raise_on_warn:
                            raise Exception('breakpoint')

        assert len(act_values[sl_name]) == len(dist_values[sl_name])

    return act_values, dist_values
