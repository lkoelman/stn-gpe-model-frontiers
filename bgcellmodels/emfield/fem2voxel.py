"""
Convert FEM data to voxel data supprted by Blender.

@author     Lucas Koelman
"""

from __future__ import division # always float division
import os.path
import numpy as np
# import scipy.spatial


def voxels_by_nearest(file_path, dxyz, out_path=None,
        fill=0.0, scale=None, translation=None, max_filesize_bytes=5e9):
    """

    @param  dxyz : list[float]
            Resolution in x,y,z dimensions : [dx, dy, dz]
    """

    # Load FEM values (x, y, z, V)
    xyzv = np.loadtxt(file_path, delimiter=',', dtype=float, unpack=False)

    if scale:
        xyzv[:, :3] *= scale

    if translation:
        xyzv[:, :3] += translation

    # Construct voxel grid
    xyz_lims = [(xyzv[:, i].min(), xyzv[:, i].max()) for i in range(3)]
    dx, dy, dz = dxyz
    nx, ny, nz = [int((xyz_lims[i][1] - xyz_lims[i][0]) / dxyz[i]) for i in range(3)]


    # Prepare binary data for voxel file
    vox_dtype = np.dtype('<i4')
    header =  np.array([nx, ny, nz, 1])
    if nx * ny * nz * vox_dtype.itemsize > max_filesize_bytes:
        raise ValueError("Binary file will exceel file size limit of {} bytes".format(max_filesize_bytes))
    vox_vals = np.full((nx * ny * nz, ), fill, dtype='f4')

    # For each sample in x,y,z,v find voxel that it maps to
    vox_numsamples = {}
    for x, y, z, v in xyzv:
        i_x = min(int((x - xyz_lims[0][0]) / dx), nx-1) # max value to voxel below
        i_y = min(int((y - xyz_lims[1][0]) / dy), ny-1)
        i_z = min(int((z - xyz_lims[2][0]) / dz), nz-1)
        vox_index = (i_z * nx * ny) + (i_y * nx) + i_x
        vox_vals[vox_index] += v
        vox_numsamples[vox_index] = vox_numsamples.get(vox_index, 0) + 1


    # Average voxels with multiple samples
    for vox_index, num_samples in vox_numsamples.iteritems():
        vox_vals[vox_index] /= num_samples

    # Rescale to [0, 1]
    vox_vals = vox_vals - np.min(vox_vals)
    vox_vals = vox_vals / np.max(vox_vals)
    
    # Write to binary file
    if out_path is None:
        out_path = os.path.splitext(file_path)[0] + '.bvox'
    with open(out_path, 'wb') as file:
        header.astype('<i4').tofile(file)
        vox_vals.astype('<i4').tofile(file)
    print("Wrote voxel data to {}".format(out_path))


def voxels_by_interp(file_path, out_path, method='linear', scale=None, translation=None):
    # TODO: use NDInterpolator : iterate over voxel centers and interpolate
    # Coordinates vary most rapidly along x, then y, then z
    # for i_z in xrange(nz):
    #     z = xyz_lims[2][0] + i_z * dz
    #     for i_y in range(ny):
    #         y = xyz_lims[1][0] + i_y * dy
    #         for i_x in range(nx):
    #             x = xyz_lims[0][0] + i_x * dx
    pass

if __name__ == '__main__':
    scale = 1e3
    fem2blend = [19935.4375, 9846.599609375, 12519.98046875]
    file_path = '/home/luye/cloudstore_m/simdata/DBS_3D_data/fem_voltages/STN_potential.txt'
    dxyz = [5.0, 5.0, 5.0]

    voxels_by_nearest(file_path, dxyz, fill=0.0, scale=scale, translation=fem2blend)