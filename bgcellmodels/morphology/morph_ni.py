"""
Transforming neuron morphologies based on neuroimaging data.

@author     Lucas Koelman
@date       26/11/2018
"""

import numpy as np
import nibabel as nib
import transforms3d

from . import morph_3d as morph3d


def densities_to_positions(masks, densities, rng=None, rotate=False,
                           erode='end'):
    """
    Generate point cloud for cell positions based on density masks.
    Also generate a random rotation for each cell.

    The mask images specify the cell density in each voxel.

    Arguments
    ---------

    @param  masks : list[Nibabel.image]
            List of 3D voxel-based images.

    @param  densities : dict[int, float]
            Cell density in [cells/mm^3] for each mask value

    @param  erode : str
            How to erode uniform grid of cells to achieve actual density
    """
    if rng is None:
        rng = np.random

    # Combine masks into one image, preserving highest value in case of overlap
    mask_sum = masks[0].get_fdata()
    for mask in masks[1:]:
        mask_data = mask.get_fdata()
        mask_sum += mask_data
        mask_sum[mask_sum > mask_data.max()] = mask_data.max()

    all_soma_locs = []
    vox_dimensions = masks[0].header.get('pixdim')[1:4] # nib.affines.voxel_sizes(mask.affine)
    vox_units = masks[0].header.get_xyzt_units()[0]
    vox_vol = np.prod(vox_dimensions)
    if vox_units != 'mm':
        raise ValueError('Voxel units must be mm and cell density must be in cells/mm^3')

    # For each mask value and corresponding density: position neurons
    for mask_val, density in densities.items():
        # Get voxels that should get this density
        vox_inds = np.array(np.where(mask_sum == mask_val)).T # coordinates as (Nx3)
        vox_centers = vox_inds # indices = center coordinates in voxel space
        num_vox = len(vox_centers)

        # Sample voxels with this density
        dx_cell_isotropic = density ** (1./3.)
        dx_vox_isotropic = vox_vol ** (1./3.) # dimensions of isotropic voxel with same volume
        if all(vox_dimensions == vox_dimensions[0]):
            xyz_dens = [dx_cell_isotropic] * 3
        else:
            # voxels are not isotropic
            xyz_dens = [vox_dimensions[i] / dx_vox_isotropic * dx_cell_isotropic for i in range(3)]
        xyz_num_cell = [int(np.ceil(dens)) for dens in xyz_dens]
        vox_num_cell = np.prod(xyz_num_cell)
        
        print('{} voxels of size {} with density of {} cells/mm^3.'.format(
              num_vox, vox_dimensions, density))
        print('Density per grid dimension: {} cells/(voxel dimension)'.format(xyz_dens))
        print('Actual cells per grid dimension: {} cells/(voxel dimension)'.format(xyz_num_cell))
        print('Total = {} voxels * {} cells/vox = {} cells'.format(
              num_vox, vox_num_cell, num_vox*vox_num_cell))

        # make 3D grid in 0-centered voxel (in voxel space)
        vox_grid_samples = [] # x,y,z grid coordinates, one dimension per row
        for i in range(3):
            dw = vox_dimensions[i]
            grid_ncell = xyz_num_cell[i]
            offset = dw / grid_ncell / 2. # move grid away from edges of voxel
            grid_pts = np.linspace(-dw/2, +dw/2, num=grid_ncell, endpoint=False) + offset
            vox_grid_samples.append(grid_pts)

        # We need each dimension repeated over the two others
        vox_grid = np.array([co.ravel() for co in np.meshgrid(*vox_grid_samples)])
        num_remove = int(vox_num_cell - np.ceil(density))
        if num_remove > 0:
            if erode == 'end':
                vox_grid = vox_grid[:, :-num_remove]
            elif erode == 'random':
                inds = rng.choice(vox_grid.shape[1], np.ceil(density), replace=False)
                vox_grid = vox_grid[:, inds]
            else:
                raise ValueError("Unknown erosion method '{}'".format(erode))

        print("Eroded grid for total of {} cells/vox using method '{}'".format(
              vox_grid.shape[1], erode))
        print('New total = {} voxels * {} cells/vox = {} cells'.format(
              num_vox, vox_grid.shape[1], num_vox*vox_grid.shape[1]))
        
        for vox_center in vox_centers:
            vox_locs = [vox_grid[i] + vox_center[i] for i in range(3)]
            
            # Transform from voxel to world coordinates
            vox_coords = np.array(vox_locs).T # one coordinate (x,y,z) per row
            world_locs = nib.affines.apply_affine(masks[0].affine, vox_coords)
            all_soma_locs.append(world_locs)

    all_soma_coords = np.concatenate(all_soma_locs, axis=0)
    num_cells = all_soma_coords.shape[0]

    # Random rotation in [0, 360] around z-axis + [-90, 90] around x-axis (intrinsic)
    # http://matthew-brett.github.io/transforms3d/reference/transforms3d.euler.html
    all_affines = []
    xz_rots = rng.random_sample((num_cells, 2)) * np.pi
    R_identity = np.eye(3)
    print('Generating random rotation for each cell.')
    for i, xz_rot in enumerate(xz_rots):
        z_angle = xz_rot[0]
        x_angle = xz_rot[1] * 2
        if rotate:
            R = transforms3d.euler.euler2mat(z_angle, x_angle, 0, axes='rzxy')
        else:
            R = R_identity
        T = all_soma_coords[i,:]
        Z = [1.0, 1.0, 1.0]
        all_affines.append(transforms3d.affines.compose(T, R, Z))

    return all_affines


def distribute_morphologies(cell_gen, **kwargs):
    """
    Distribute NEURON cell models in 3D space using function
    `densities_to_positions`.

    @see    densities_to_positions
            See argument description for keyword arguments that will
            be passed to this function.

    @param  cell_gen : callable
            Generator of new NEURON cell model with somatic compartment
            centered at (0, 0, 0) in 3D space.
    """
    max_num = kwargs.pop('max_num', 100)
    cell_affines = densities_to_positions(**kwargs)
    cells = []
    for i, A in enumerate(cell_affines):
        if i >= max_num:
            print('Stopping at {}/{} cell positions'.format(max_num, len(cell_affines)))
            break
        icell = cell_gen()
        print('Positioning cell {}'.format(i))
        morph3d.transform_sections(icell.all, A)
        cells.append(icell)
    return cells


def load_streamlines(file_path, max_num=1e12, min_length=0.0, indices=None):
    """
    Load streamlines from file

    Arguments
    ---------

    @param  indices : list[int]
            Indices of streamlines in tractogram file that should be used for
            aligning model axons.
    """
    tck_file = nib.streamlines.load(file_path, lazy_load=True)

    # Make sure tracts are defined in RAS+ world coordinate system
    tractogram = tck_file.tractogram.to_world(lazy=True)

    # Manual transformation to RAS+ world coordinate system
    # vox2ras = tck_file.tractogram.affine_to_rasmm
    # tck_ras_coords = nib.affines.apply_affine(vox2ras, streamline)

    streamlines_filtered = []
    for i, streamline in enumerate(tractogram.streamlines): # lazy-loading generator
        # streamline is (N x 3) matrix
        if indices and (i not in indices):
            continue
        if len(streamlines_filtered) >= max_num:
            break
        # check length
        if min_length > 0:
            tck_len = np.sum(np.linalg.norm(np.diff(streamline, axis=0), axis=1))
        else:
            tck_len = 1.0
        if tck_len >= min_length:
            streamlines_filtered.append(streamline)

    return streamlines_filtered
        


if __name__ == '__main__':

    # Test voxel positioning
    # mask_paths = [
    #     '/home/luye/workspace/fem-neuron-interface/test_data/hcp-subj_118932_ROI_box.nii.gz'
    # ]
    # masks = [nib.load(path) for path in mask_paths]
    # densities = {1 : 50.386 }

    # # Make cell generator function
    # import bgcellmodels.models.STN.GilliesWillshaw.gillies_model as gillies
    # from neuron import h
    # def centered_cell():
    #     icell = gillies.stn_cell_standardized()
    #     soma = icell.soma[0]
    #     # Automatically generate 3D layout for cell without morphology.
    #     # Using this method first point of soma is centered at (0, 0, 0)
    #     h.define_shape(sec=soma)
    #     return icell

    # cells = distribute_morphologies(centered_cell, masks=masks, densities=densities, max_num=100)

    # # Write cells to STEP file
    # import morphio
    # all_sec = sum((list(icell.all) for icell in cells), [])
    # tot_compartments = sum((sec.nseg for sec in all_sec))
    # print('Writing STEP file for {} cells ({} compartments)'.format(
    #       len(cells), tot_compartments))
    # step_outfile = '/home/luye/workspace/fem-neuron-interface/test_data/cells_in_mask.stp'
    # morphio.write_STEP(all_sec, step_outfile)
    pass
