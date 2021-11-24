"""
Geometrical operations on neuron morphologies.

@author     Lucas Koelman
@date       26/11/2018


Notes
-----

Additional libraries to deal with NEURON morphologies:

- PyNeuron-Toolbox:https://github.com/ahwillia/PyNeuron-Toolbox
- btmorph : https://btmorph.readthedocs.io/en/latest/index.html#
- NeuroM : https://github.com/BlueBrain/NeuroM
- AllenSDK : https://github.com/AllenInstitute/AllenSDK
- Hoc2Swc : https://github.com/JustasB/hoc2swc/
"""

from neuron import h
import numpy as np


def get_section_samples(section_lists, include_diam=True):
    """
    Return 3D sample points for all sections as Nx4 matrix.

    @return     N x 4 numpy array with samples (x, y, z, diam) as rows.
    """
    all_samples = []
    sections_numsample = []
    seclists_numsec = []

    for sections in section_lists:
        num_sec = 0
        for sec in sections:
            num_sec += 1
            num_samples = int(h.n3d(sec=sec))
            if include_diam:
                sec_samples = [
                    (h.x3d(i, sec=sec),
                     h.y3d(i, sec=sec),
                     h.z3d(i, sec=sec),
                     h.diam3d(i, sec=sec))
                        for i in range(num_samples)
                ]
            else:
                sec_samples = [
                    (h.x3d(i, sec=sec),
                     h.y3d(i, sec=sec),
                     h.z3d(i, sec=sec))
                        for i in range(num_samples)
                ]
            all_samples.extend(sec_samples)
            sections_numsample.append(num_samples)

        seclists_numsec.append(num_sec)

    return all_samples, sections_numsample, seclists_numsec


def get_segment_centers(section_lists, samples_as_rows=False):
    """
    Calculate segment center coordinates for Section that has 3d sample points 
    assigned to it.

    @param  section_lists
            Individual iterables of neuron.Section (can be list, SectionList, etc.)

    @pre    All section in sectionlist have 3d sample points assigned, either
            using Hoc.pt3dadd() or Hoc.define_shape()
    """
    x_allsec = []
    y_allsec = []
    z_allsec = []
    sections_numsample = []
    seclists_numsec = []

    for sections in section_lists:
        num_sec = 0
        for sec in sections:
            num_sec += 1
            num_samples = int(h.n3d(sec=sec))
            nseg = sec.nseg
            sections_numsample.append(nseg + 2)

            # Get 3D sample points for section
            xx = h.Vector([h.x3d(i, sec=sec) for i in range(num_samples)])
            yy = h.Vector([h.y3d(i, sec=sec) for i in range(num_samples)])
            zz = h.Vector([h.z3d(i, sec=sec) for i in range(num_samples)])

            # Length in micron from start of section to sample i
            pt_locs = h.Vector([h.arc3d(i, sec=sec) for i in range(num_samples)])
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

            x_allsec.extend(list(node_xlocs))
            y_allsec.extend(list(node_ylocs))
            z_allsec.extend(list(node_zlocs))

        seclists_numsec.append(num_sec)

    if samples_as_rows:
        return zip(x_allsec, y_allsec, z_allsec), sections_numsample, seclists_numsec
    else:
        return (x_allsec, y_allsec, z_allsec), sections_numsample, seclists_numsec


def find_closest_section(point3d, sections, measure_from='segment_centers'):
    """
    Find section closest to target point

    @param  measure_from : str
            Either 'segment_centers' or 'pt3d'

    @return i_sec, i_pt3d : tuple[int, int]
            Index of closest section in section list 'sections' and index of the
            node or 3d sample point in the section that is closest to the given
            point (depending on argument 'measure_from')
    """
    if measure_from == 'segment_centers':
        node_pt3d, node_n3d, _ = get_segment_centers(
                                        [sections], samples_as_rows=True)
    elif measure_from == 'pt3d':
        node_pt3d, node_n3d, _ = get_section_samples(
                                        [sections], include_diam=False)
    else:
        raise ValueError(measure_from)
    node_pt3d = np.array(node_pt3d)
    node_pt_idx_upper = np.cumsum(node_n3d)
    node_pt_dists = np.linalg.norm(node_pt3d - point3d, axis=1)
    i_pt3d = np.argmin(node_pt_dists)
    i_sec = next((i for i, hi in enumerate(node_pt_idx_upper) if (i_pt3d < hi)))

    # Index of the segment (including 0/1) or 3d sample point in the section
    parent_pt0_idx = node_pt_idx_upper[i_sec-1] if i_sec > 0 else 0
    i_pt3d_in_sec = i_pt3d - parent_pt0_idx
    return i_sec, i_pt3d_in_sec


def transform_sections(secs, A):
    """
    Apply transformation to sections

    @param  secs : list[neuron.Section]
            List of NEURON sections.

    @param  A : np.array
            4x4 transformation matrix in column-major layout
    """
    for sec in secs:
        # Construct matrix with vertices as rows
        num_verts = int(h.n3d(sec=sec))
        src_verts = np.array([
            [h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec), 1.0]
            for i in range(num_verts)])

        # Transform vertex matrix
        new_verts = np.dot(src_verts, A.T)

        # Update 3D info
        for i in range(num_verts):
            diam = h.diam3d(i, sec=sec)
            h.pt3dchange(i, new_verts[i, 0], new_verts[
                         i, 1], new_verts[i, 2], diam, sec=sec)


def perturb_sample_angles(samples, max_angle, rng=None, seed=None):
    """
    At each sample, rotate the subtree around the sample coordinates

    @pre    requires module 'transforms3d' installed
            e.g. pip install transforms3d

    @param  samples : list[list[float/int]]
            List of SWC samples, each represented by a list of 7 numbers.

    @param  max_angles : tuple[float, float]


    @return perturbed_sampled : list[list[float/int]]
            Samples with rotations applied to subtrees
    """
    from transforms3d import axangles

    # Each time there is a 'break', rotate all children around the break-point.
    # You can keep track of the accumulated rotation matrices
    if rng is None:
        rng = np.random
    if seed:
        rng.seed(seed)

    # Copy samples so we don't modify originals
    samples = [list(sample) for sample in samples]
    root_sample = next((s for s in samples if s[6] == -1))

    subtree_affines = {} # transformation matrices for subtree of sample i
    stack = [root_sample]

    while len(stack) > 0:
        # NOTE: list is pass-by-reference so modification modies samples
        sample = stack.pop()

        # Collect transformation matrices of all samples on path to root
        ancestor_affines = []
        next_parent_idx = sample[6]
        while next_parent_idx > 0:
            if next_parent_idx in subtree_affines:
                ancestor_affines.append(subtree_affines[next_parent_idx])
            next_ancestor = next((s for s in samples if s[0] == next_parent_idx))
            next_parent_idx = next_ancestor[6]

        # Multiply matrices from left side, in order from root to leaves
        total_xform = reduce(lambda old, new: np.dot(new, old),
                             reversed(ancestor_affines), np.eye(4))

        sample[2:5] = np.dot(total_xform, sample[2:5] + [1.0])[0:3]
        if any(np.isnan(sample)):
            raise Exception('NaN resulting from multiplication')


        # handle sample
        parent_sample = next((s for s in samples if s[0] == sample[6]), None)
        if parent_sample is None:
            # Use first child sample to construct axis
            parent_sample = next((s for s in samples if s[6] == sample[0]))

        # Construct local coordinate system
        sample_coords = np.array(sample[2:5])
        parent_coords = parent_sample[2:5]
        if not np.allclose(sample_coords, parent_coords):

            # Find the sample axis
            sample_axis = sample_coords - parent_coords
            sample_axis_unit = sample_axis / np.sqrt(np.dot(sample_axis, sample_axis))

            # Choose random axis perpendicular to sample axis
            vec = rng.random(3) - 0.5
            vec_sample_component = np.dot(vec, sample_axis_unit) * sample_axis_unit
            perp_axis = vec - vec_sample_component

            # Make rotation matrix for all samples in subtree of current sample
            x = rng.random(1) - 0.5
            angle = 2 * x[0] * max_angle
            subtree_affines[sample[0]] = axangles.axangle2aff(perp_axis, angle, point=sample_coords)

        # Find child samples
        children = (s for s in samples if s[6] == sample[0])
        for child_sample in children:
            stack.append(child_sample)

    return samples


if __name__ == '__main__':
    swc_file = '/home/luye/workspace/bgcellmodels/bgcellmodels/models/STN/Miocinovic2006/morphologies/gillies_original.swc'

    import os.path, morph_io
    samples = morph_io.read_SWC_samples(swc_file)
    max_angle = np.pi/6
    max_angle_deg = np.rad2deg(max_angle)

    for seed in range(10):
        perturbed_samples = perturb_sample_angles(samples,
                                max_angle=max_angle, seed=seed)
        
        out_file = swc_file[:-4] + '_alpha-{:.0f}deg_seed-{}.swc'.format(
                        max_angle_deg, seed)

        morph_io.write_SWC_samples(perturbed_samples, out_file,
            comment="Original file {} with angle perturbations of {} degrees".format(
                     os.path.basename(swc_file), max_angle_deg))