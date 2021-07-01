import bpy
import numpy as np
from mathutils import Vector

C = bpy.context

scale = 25

file_path = '/home/luye/cloudstore_m/simdata/DBS_3D_data/fem_voltages/STN_potential.txt'
scale = 1e3
translation = [19935.4375, 9846.599609375, 12519.98046875]

xyzv = np.loadtxt(file_path, delimiter=',', dtype=float, unpack=False)
if scale:
    xyzv[:, :3] *= scale
if translation:
    xyzv[:, :3] += translation

# Get vertices
data = xyzv[:, :3]
print('vertcount = {}'.format(xyzv.shape[0]))

# Create and arrange mesh data
verts = [ Vector(data[i, :3]) for i in range(data.shape[0]) ]
m = bpy.data.meshes.new('FEM_point_cloud')
m.from_pydata(verts, [], [])

# Create mesh object and link to scene collection
mesh_obj = bpy.data.objects.new('FEM_point_cloud', m)
C.scene.objects.link(mesh_obj)

# Blender 2.8 : instancing of objects
# https://docs.blender.org/manual/en/latest/scene_layout/object/properties/instancing/index.html
# # Add minimal icosphere
# bpy.ops.mesh.primitive_ico_sphere_add(subdivisions = 1, size=0.05)
# sph_obj = bpy.data.objects[ C.object.name ]

# # Set instancing props
# for ob in sph_obj, mesh_obj:
#     ob.instance_type               = 'VERTS'
#     ob.show_instancer_for_viewport = False
#     ob.show_instancer_for_render   = False

# # Set instance parenting (parent icosphere to verts)
# mesh_obj.select_set(True)
# C.view_layer.objects.active = mesh_obj

# bpy.ops.object.parent_set(type='VERTEX', keep_transform=True)