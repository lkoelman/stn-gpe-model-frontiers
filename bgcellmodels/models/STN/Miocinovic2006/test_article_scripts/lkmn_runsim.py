"""
Run Miocinovic STN neuron with extracellular stimulation

@author Lucas Koelman
@date   06-04-2017
@note   must be run from script directory or .hoc files not found

"""

# Make sure other modules are on Python path
import sys, os.path
scriptdir, scriptfile = os.path.split(__file__)

# Load NEURON
import neuron
from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# Load own NEURON mechanisms
NRN_MECH_PATH = os.path.normpath(os.path.join(scriptdir, 'nrn_mechs'))
neuron.load_mechanisms(NRN_MECH_PATH)

# Third party modules
import numpy as np

# library:      pypi.python.org/pypi/transforms3d
from transforms3d import utils as tfutils
from transforms3d import affines


# Global variables for Miocinovic file
stn_groups_nsec = {
    'soma': 'somaelements',
    'dend': 'dendelements',
    'initseg': 'iselements',
    'node': 'axonnodes',
    'MYSA': 'paranodes1',
    'FLUT': 'paranodes2',
    'STIN': 'axoninter',
}

def transform_sections(secs, TR_hom):
    """
    Apply transformation to sections
    
    @param secs         list of Section
    @param TR_hom       4x4 transformation matrix in column-major layout
    """
    for sec in secs:
        sec.push() # make section CAS

        # Construct matrix with vertices as rows
        num_verts = int(h.n3d())
        src_verts = np.array([[h.x3d(i), h.y3d(i), h.z3d(i), 1.0] for i in xrange(num_verts)])

        # Transform vertex matrix
        new_verts = np.dot(src_verts, TR_hom.T)

        # Update 3D info
        for i in xrange(num_verts):
            h.pt3dchange(i, new_verts[i,0], new_verts[i,1], new_verts[i,2], h.diam3d(i))

        h.pop_section()

def read_V_raw(x, y, z):
    """
    Query voltage attenuation factor relative to 1V peak-to-peak at location (x,y,z)
    """
    # TODO: replace this by Diego's function
    return 1.0

def read_transform():
    """ Query transformation matrix to be applied to cells """
    # TODO: read matrix file or hard-code it here
    T = [2.0, 2.0, 2.0]
    R = [[0, -1, 0], [1, 0, 0], [0, 0, 1]] # rotation matrix
    Z = [1.0, 1.0, 1.0]
    TR = affines.compose(T, R, Z, S=None)
    return TR


def write_geometry_STEP(secs, filepath):
    """
    Write write_geometry as sections as polyline

    @param secs     SectionList containing sections whose geometry will
                    be exported. It is important that this is a SectionList
                    so that the CAS will be set correctly during iterations.
    """
    # library:      github.com/tpaviot/pythonocc-core
    # wrapper API:  api.pythonocc.org/py-modindex.html
    # c++ API:      www.opencascade.com/doc/occt-7.1.0/refman/html/toolkit_tkmath.html
    from OCC.gp import gp_Pnt, gp_Lin, gp_Ax1, gp_Dir # core geometry types
    from OCC.TColgp import TColgp_Array1OfPnt # collections
    from OCC.GeomAPI import GeomAPI_PointsToBSpline # geometry types
    from OCC.GeomAbs import GeomAbs_C0
    from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge # geometry to shape
    from OCC.TopoDS import TopoDS_Compound, TopoDS_Builder # Compound shapes
    from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCC.IFSelect import IFSelect_RetDone # Return codes
    from OCC.Interface import Interface_Static_SetCVal

    # TODO: give mesh to Diego, ask him to position it, then give me matrix
    # NOTE: for mesh generation, see https://github.com/MetaCell/NEURON-UI/blob/master/neuron_ui/neuron_geometries_utils.py
    step_writer = STEPControl_Writer()
    Interface_Static_SetCVal("write.step.schema", "AP203")
    c_edges = TopoDS_Compound() # group edges in a compound shape
    c_builder = TopoDS_Builder()
    c_builder.MakeCompound(c_edges)

    # Build polyline/curve using 3D location info
    for sec in secs:
        # Copy points
        num_verts = int(h.n3d())
        pts = TColgp_Array1OfPnt(1, num_verts)
        for i in xrange(num_verts):
            pts.SetValue(i+1, gp_Pnt(h.x3d(i), h.y3d(i), h.z3d(i)))

        # Build curve
        crv_builder = GeomAPI_PointsToBSpline(pts, 1, 1, GeomAbs_C0, 1e-3) # No continuity, max degree 1
        crv = crv_builder.Curve() # this is a Handle/reference
        edge = BRepBuilderAPI_MakeEdge(crv).Edge()

        # Add to compound
        # step_writer.Transfer(edge, STEPControl_AsIs) # accepts TopoDS_Shape
        c_builder.Add(c_edges, edge)

    # Write CAD geometry to STEP file
    step_writer.Transfer(c_edges, STEPControl_AsIs) # accepts TopoDS_Shape
    status = step_writer.Write("stn_tree_polylines.stp")
    if status != IFSelect_RetDone:
        raise Exception("Failed to write STEP file")


def export_cell():
    """ Export cell morphology used in Miocinovic simulation """
    # Set up model and functions
    h.xopen("lkmn_run_setup.hoc")

    # Get STN tree
    somasecs = getattr(h, 'soma')
    sl = h.SectionList(sec=somasecs[0])
    sl.wholetree()

    # Export geometry
    filepath = "stn_tree.stp" # write as STEP file (ISO standard for CAD exhange)
    write_geometry_STEP(sl, filepath)


def inspect_model():
    """ Export cell morphology used in Miocinovic simulation """
    # Set up model and functions
    h.xopen("lkmn_run_setup.hoc")
    from neuron import gui


def run_example():
    """ Run example with single STN cell """

    # Set up model and functions using Miocinovic code
    h.xopen("lkmn_run_setup.hoc")

    # Make rigid transformation matrix
    TR = read_transform()

    # Apply transformation to cell
    for group_name, nsec in stn_groups_nsec:
        secs = getattr(h, group_name) # this is a Section[]
        transform_sections(secs, TR)

    # fill V_raw by evaluating function at each updated compartment location
    nsec_stn = getattr(h, 'total') # sum of all values in above dict (see n17*.hoc)
    V_raw = getattr(h, 'V_raw')
    stn_secrefs = getattr(h, 's') # objref s[] = SectionRef[total] (defined in initcell() in n17*.hoc)
    jj_cell = 1 # cell we wish to modify
    for kk_sec in xrange(nsec_stn):
        # Get centroid of section coordinates (all nseg=1)
        secref = stn_secrefs[kk_sec]
        secref.sec.push()

        # Get centroid of section coordinates (all nseg=1)
        num_verts = int(h.n3d())
        sec_verts = np.array([[h.x3d(i), h.y3d(i), h.z3d(i)] for i in xrange(num_verts)])
        centroid = np.mean(sec_verts, axis=0)

        # Query voltage attenuation at this location
        V_raw.x[jj_cell*nsec_stn + kk_sec] = read_V_raw(*centroid)
        h.pop_section()

    # TODO: simulate model

    # Make local variables global
    globals().update(locals())

if __name__ == '__main__':
    # export_cell()
    inspect_model()