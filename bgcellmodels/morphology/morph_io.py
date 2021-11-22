"""
Reading and writing of neuron morphologies.


Notes
-----

Additional libraries for working with NEURON morphologies:

- PyNeuron-Toolbox: https://github.com/ahwillia/PyNeuron-Toolbox
- btmorph :         https://btmorph.readthedocs.io/en/latest/index.html#
- NeuroM :          https://github.com/BlueBrain/NeuroM
- AllenSDK :        https://github.com/AllenInstitute/AllenSDK
- Hoc2Swc :         https://github.com/JustasB/hoc2swc/
- NeuroMorphoVis:   https://github.com/BlueBrain/NeuroMorphoVis
"""
# Standard library
import json, io, re, os

# Third party libraries
from neuron import h
import numpy as np
from numpy.lib import recfunctions as rfn

# Custom libraries
from bgcellmodels.common import nrnutil
from bgcellmodels.common.treeutils import parent, parent_loc
from . import morph_3d


def morphology_to_dict(sections):
    """
    Extract morphology info from given sections.

    @see        save_json() and load_json()

    @return     list(dict()) containing one dict for each section:
                the dict contains its morphological & topological information

    @note       code modified from R.A. McDougal's post at
                https://www.neuron.yale.edu/phpBB/viewtopic.php?f=2&t=3478&p=14758
    """

    # Assign index to each sections (for expressing topology relations)
    section_map = {sec: i for i, sec in enumerate(sections)}
    
    # adds 3D info using simple algorithm if not present
    h.define_shape()
    
    result = []
    for sec in sections:
        parent_sec = parent(sec)

        parent_x = -1 if parent_sec is None else parent_loc(sec, parent_sec)
        parent_id = -1 if parent_sec is None else section_map[parent_sec] # get parent index
        
        n3d = int(h.n3d(sec=sec))
        
        result.append({
            'section_orientation':  h.section_orientation(sec=sec),
            'parent':               parent_id,
            'parent_loc':           parent_x,
            'x':                    [h.x3d(i, sec=sec) for i in range(n3d)],
            'y':                    [h.y3d(i, sec=sec) for i in range(n3d)],
            'z':                    [h.z3d(i, sec=sec) for i in range(n3d)],
            'diam':                 [h.diam3d(i, sec=sec) for i in range(n3d)],
            'name':                 sec.hname()           
        })
    
    return result


def cell_to_dict(section, descr=None, metadata=None, icell=None):
    """
    Save data required to reconstruct the cell in NEURON.

    @pre    requires NEURON >= 7.6 for function Hoc.Section.psection()
    """
    if icell is None:
        icell = section.cell() # set by NEURON

    cell_data = {
        'description'       : descr,
        'metadata'          : metadata,
        'section_data'      : [],
        'cell_sectionlists' : {},
    }

    sl = h.SectionList()
    sl.wholetree(sec=section)
    seclist_bfs = list(sl)
    section_map = {sec: i for i, sec in enumerate(seclist_bfs)}
    
    for sec in seclist_bfs:

        # Connection to parent section
        parent_sec = parent(sec)
        parent_x = -1 if parent_sec is None else parent_loc(sec, parent_sec)
        parent_idx = -1 if parent_sec is None else section_map[parent_sec] # get parent index

        # Get section data
        sec_data = sec.psection() # requires NEURON >= 7.6

        # Delete raw Hoc objects
        del sec_data['cell']
        del sec_data['morphology']['parent']
        del sec_data['morphology']['trueparent']
        for k in sec_data['point_processes'].keys():
            sec_data['point_processes'][k] = [
                str(pp) for pp in sec_data['point_processes'][k]
            ]

        # Data not included in Section.psection() dict
        sec_data['morphology'].update({
            'parent_idx'    : parent_idx,
            'parent_loc'    : parent_x,
            'section_orientation':  h.section_orientation(sec=sec),
        })

        for ion_name, ion_data in sec_data['ions'].items():
            ion_data['ion_style'] = h.ion_style(ion_name+'_ion', sec=sec)

        cell_data['section_data'].append(sec_data)

    # Default SectionList membership (reconstruct using section names)
    if icell is not None:
        named_seclists = cell_data['cell_sectionlists']
        proto_seclist_names = list(nrnutil.nrn_proto_seclists_arrays.keys())
        proto_seclist_names.remove('all')
        for seclist_name in proto_seclist_names:
            if hasattr(icell, seclist_name):
                named_seclists[seclist_name] = [
                    sec.name() for sec in getattr(icell, seclist_name)
                ]

    return cell_data



def cell_from_dict(
        cell_data,
        hoc_sections=False,
        name_prefix=None,
        name_substitutions=None):
    """
    Recreate cell from saved data

    @see    cell_to_dict

    @param  name_substitutions : list[tuple[str, str/function]]
            Substitutions to make in original name, as first and second
            arguments to function re.sub. For example:
            [
                ('.', '_'),
                (r"\\[(\\d+)\\]", lambda m: m.groups()[0])
            ]

    @return cell_seclists : dict[str, list[h.Section]]
            Section lists containing instantiated sections
    """
    if name_prefix is None:
        name_prefix = ''
    if name_substitutions is None:
        name_substitutions = []

    # Sections by section list
    seclists = {}
    seclists['all'] = sections = []
    for sl_name in cell_data['cell_sectionlists'].keys():
        seclists[sl_name] = []

    # Create all sections
    for sec_data in cell_data['section_data']:
        # make section
        secname = name_prefix + sec_data['name']
        for pattern, repl in name_substitutions:
            secname = re.sub(r"[\[\]\.]", "", secname)
        if hoc_sections:
            created = h("create %s" % secname)
            if created != 1:
                raise Exception("Could not create section with name '{}'".format(secname))
            sec = getattr(h, secname)
        else:
            sec = h.Section(name=secname)
        sections.append(sec)

        # Geometry
        if len(sec_data['morphology']['pts3d']) > 0:
            for x, y, z, d in sec_data['morphology']['pts3d']:
                h.pt3dadd(x, y, z, d, sec=sec)
            sec.nseg    = sec_data['nseg']
        else:
            sec.L       = sec_data['morphology']['L']
            sec.nseg    = sec_data['nseg']
            for i, seg in enumerate(sec):
                seg.diam = sec_data['morphology']['diam'][i]
        
        # Passive biophysical properties
        sec.Ra = sec_data['Ra']
        for i, seg in enumerate(sec):
            seg.cm = sec_data['cm'][i]

        # RANGE properties
        for mech_name, mech_params in sec_data['density_mechs'].items():
            sec.insert(mech_name)
            for pname, pvals in mech_params.items():
                for i, seg in enumerate(sec):
                    setattr(seg, pname+'_'+mech_name, pvals[i])

        # Ions
        for ion_name, ion_data in sec_data['ions'].items():
            flags = nrnutil.make_ion_style_flags(ion_data['ion_style'])
            h.ion_style(ion_name+'_ion', *flags, sec=sec)
            ion_erev = 'e' + ion_name
            # if flags[1] == 1 or flags[1] == 2:
            for i, seg in enumerate(sec):
                setattr(seg, ion_erev, ion_data[ion_erev][i])

        # Section lists
        for sl_name, sl_members in cell_data['cell_sectionlists'].items():
            if sec_data['name'] in sl_members:
                seclists[sl_name].append(sec)

    # Connect sections following topology
    # NOTE: reversed(...) will make sure that order of children() is same
    for sec, sec_data in reversed(zip(sections, cell_data['section_data'])):
        if sec_data['morphology']['parent_loc'] >= 0:
            parent_sec  = sections[sec_data['morphology']['parent_idx']]
            parent_loc  = sec_data['morphology']['parent_loc']
            orient_loc  = sec_data['morphology']['section_orientation']
            sec.connect(parent_sec(parent_loc), orient_loc)

    return seclists



def morphology_to_SWC(sections, filename):
    """
    Convert all instantiated NEURON cells to SWC. 

    Each isolated tree is written to a separate .swc file.

    Instead use:
    - built in exporters in ModelView (to NeuroML)
    - https://github.com/JustasB/hoc2swc
    - http://neuronland.org/NLMorphologyConverter/NLMorphologyConverter.html
    """
    from hoc2swc import neuron2swc # this function exports all loaded cells to SWC
    neuron2swc(filename)


def read_SWC_samples(file_path, structured=True):
    """
    Read samples in SWC file.

    @return     samples : list[list[<7 sample elements>]]

                A sample consists of 7 elements, i.e. [sample_number (int), 
                structure_identifier (int), x (float), y (float), z (float),
                radius (float), parent_sample_number (int)]

    """
    # sample_data_types = [int, int, float, float, float, float, int]
    # samples = []
    # with open(file_path, 'r') as swc_file:
    #     for line in swc_file:
    #         if line.startswith('#'):
    #             continue

    #         sample_data = line.split() # default using whitespace
    #         sample_parsed = [
    #             value_type(sample_data[i]) for i, value_type in enumerate(
    #                 sample_data_types)
    #         ]
    #         samples.append(sample_parsed)
    #         # sample_idx, sample_type, x, y, z, radius, parent_idx = sample_parsed

    swc_dtype = [('sample_id', 'i4'), ('structure_id', 'i4'),
                 ('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('radius', 'f4'),
                 ('parent_id', 'i4')]

    if structured:
        target_dtype = swc_dtype
    else:
        target_dtype = float

    return np.loadtxt(file_path, dtype=target_dtype)


def write_SWC_samples(samples, file_path, comment=None):
    """
    Write samples to SWC file.
    """
    sample_format = "{:.0f} {:.0f} {:f} {:f} {:f} {:f} {:.0f}\n"
    with open(file_path, 'w') as file:
        if comment:
            file.writelines(
                ['# ' + comment + '\n' for line in comment.splitlines()])
        for sample in samples:
            file.write(sample_format.format(*sample))


def save_json(sections, filename, encoding='utf-8'):
    """
    Save morphology to JSON
    """
    morph_dicts = morphology_to_dict(sections)
    
    # json_string = json.dumps(morph_dict, indent=2)

    if encoding == 'ascii':
        with open(filename, 'w') as f:
            json.dump(morph_dicts, f)

    else:
        with io.open(filename, 'w', encoding=encoding) as f:
            f.write(json.dumps(morph_dicts, ensure_ascii=False))


def load_json(morphfile):
    """
    Load morphology from JSON.
    """

    with open(morphfile, 'r') as f:
        secdata = json.load(f)

    seclist = []
    for sd in secdata:
        # make section
        sec = h.Section(name=sd['name'])
        seclist.append(sec)

        # make 3d morphology
        for x,y,z,d in zip(sd['x'], sd['y'], sd['z'], sd('diam')):
            h.pt3dadd(x, y, z, d, sec=sec)

    # connect children to parent compartments
    for sec,sd in zip(seclist,secdata):
        if sd['parent_loc'] >= 0:
            parent_sec = seclist[sd['parent']] # not parent_loc, grab from sec_list not sec 
            sec.connect(parent_sec(sd['parent_loc']), sd['section_orientation'])

    return seclist


def test_json_export():
    """
    Text exporting a simple example morphology to JSON format
    """

    s = [h.Section(name='s[%d]' % i) for i in range(13)]

    """
        Create the tree
       
              s0
        s1    s2         s3
        s4           s5      s6
        s7         s8 s9       s10
    """
    for p, c in [[0, 1], [0, 2], [0, 3], [1, 4], [4, 7], [3, 5], [3, 6], [5, 8], [5, 9], [6, 10]]:
        s[c].connect(s[p])
   
    print json.dumps(morphology_to_dict([s[3], s[5], s[8], s[0], s[1], s[4], s[7]]), indent=2)


def prepare_hoc_morphology_for_SWC_export(hoc_script_path, hoc_out_path):
    """
    Replace expressions in brackets (not containing variables) by their
    evaluated values.

    Usage
    -----

    - first replace all variable names in 'connect' 'access' and 'pt3dadd' statements
    - ensure there is only one statement per line
        - split 'create A, B, C' into separate statements
    """
    hoc_clean_morph = ''
    def match_evaluator(match):
        """ process a regex match and return the corrected line """
        assert match.lastindex == 1 # only one match
        expr = match.group(0)
        repl = str(eval(expr))
        print("{}\n>>>>>>\n{}".format(expr, repl))
        return repl

    pattern = r"\[([\d\*\+]+)\]" # expression within brackets e.g. [i*j+2]

    with open(hoc_script_path, 'r') as hoc_script:
        hoc_dirty_morph = hoc_script.read()

    hoc_clean_morph = re.sub(pattern, match_evaluator, hoc_dirty_morph)

    with open(hoc_out_path, 'w') as out_script:
        out_script.write(hoc_clean_morph)


def morphology_to_STEP_1D(secs, filepath):
    """
    Write morphology as one-dimensional STEP file.

    @param  secs : neuron.SectionList

            SectionList containing sections whose geometry will be exported. 
            It is important that this is a SectionList so that the CAS will be
            set correctly during iterations.
    """
    if not filepath.endswith('.stp'):
        filepath += '.stp'

    # library:      github.com/tpaviot/pythonocc-core
    # wrapper API:  api.pythonocc.org/py-modindex.html
    # c++ API:      www.opencascade.com/doc/occt-7.1.0/refman/html/toolkit_tkmath.html
    from OCC.gp import gp_Pnt # gp_Lin, gp_Ax1, gp_Dir # core geometry types
    from OCC.TColgp import TColgp_Array1OfPnt # collections
    from OCC.GeomAPI import GeomAPI_PointsToBSpline # geometry types
    from OCC.GeomAbs import GeomAbs_C0
    from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge # geometry to shape
    from OCC.TopoDS import TopoDS_Compound, TopoDS_Builder # Compound shapes
    from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCC.IFSelect import IFSelect_RetDone # Return codes
    from OCC.Interface import Interface_Static_SetCVal
    
    step_writer = STEPControl_Writer()
    Interface_Static_SetCVal("write.step.schema", "AP203")
    c_edges = TopoDS_Compound() # group edges in a compound shape
    c_builder = TopoDS_Builder()
    c_builder.MakeCompound(c_edges)

    # Build polyline/curve using 3D location info
    for sec in secs:
        # Copy points
        num_verts = int(h.n3d(sec=sec))
        pts = TColgp_Array1OfPnt(1, num_verts)
        for i in range(num_verts):
            pts.SetValue(i+1,
                gp_Pnt(h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)))

        # Build curve
        crv_builder = GeomAPI_PointsToBSpline(pts, 1, 1, GeomAbs_C0, 1e-3) # No continuity, max degree 1
        crv = crv_builder.Curve() # this is a Handle/reference
        edge = BRepBuilderAPI_MakeEdge(crv).Edge()

        # Add to compound
        # step_writer.Transfer(edge, STEPControl_AsIs) # accepts TopoDS_Shape
        c_builder.Add(c_edges, edge)

    # Write CAD geometry to STEP file
    step_writer.Transfer(c_edges, STEPControl_AsIs) # accepts TopoDS_Shape
    status = step_writer.Write(filepath)
    if status != IFSelect_RetDone:
        raise Exception("Failed to write STEP file")


def morphologies_to_edges(section_lists, segment_centers=True,
                         scale=1.0, translation=None, transform=None,
                         flatten_cells=True):
    """
    Convert neuron morphologies to lists of vertices and edges.

    Each segment is represented by a degenerate face, i.e. a triangle
    with vertices [a b a].

    @pre    requires package plyfile, e.g. `pip install plyfile`

    @see    morphology_to_PLY()

    @return (vertices, edges) : tuple
            
            - If flatten_cells is True, vertices is an Nx3 numpy array containing
            all the compartment coordinates, and edges are pairs of indices
            into this array.
            
            - If flatten_cells is False, vertices is a list of Nx3 arrays, 
            one for each SectionList, and edges a list of edges for each vertex
            array.
    """

    # Get 3D samples
    if segment_centers:
        verts_data = morph_3d.get_segment_centers(section_lists, samples_as_rows=True)
    else:
        verts_data = morph_3d.get_section_samples(section_lists, include_diam=False)

    samples_xyz, secs_num3d, seclist_numsec = verts_data
    samples_xyz = np.array(samples_xyz)

    # Apply transformation before writing
    if (translation is not None) or (transform is not None) or (scale != 1.0):
        samples_mat = np.ones((len(samples_xyz), 4))
        samples_mat[:,:3] = samples_xyz
        A = np.array(transform) if transform else np.eye(4)
        if (translation is not None):
            A[:-1, 3] += translation
        samples_mat = np.dot(samples_mat, A.T)
        if scale != 1.0:
            samples_mat *= scale
        samples_xyz = samples_mat[:, :3]

    # Create vertices
    if flatten_cells:
        all_verts = samples_xyz  # flat vertex list
    else:
        all_verts = []           # vertex list per cell
    
    # Create edges
    all_edges = []
    sample_offset = 0
    seclist_index = 0
    seclist_bounds = np.cumsum(seclist_numsec)

    for i_sec, num_samples in enumerate(secs_num3d):

        if flatten_cells:
            edge_list = all_edges
            vert_offset = sample_offset
        else:
            # Start a new vertex and edge list if section belongs to new SectionList
            first_seclist = (i_sec == 0)
            new_seclist = (i_sec >= seclist_bounds[seclist_index])
            if first_seclist or new_seclist:
                all_edges.append([])
                all_verts.append([])
                num_verts_added = 0
                if new_seclist:
                    seclist_index += 1
            
            vert_offset = num_verts_added
            num_verts_added += num_samples
            edge_list = all_edges[-1]
            vert_list = all_verts[-1]
            sec_samples = samples_xyz[sample_offset:(sample_offset+num_samples), :]
            vert_list.extend(sec_samples)
        
        # Add one edge between each successive pair of samples in the section
        edges = [
            (i, i+1) for i in range(vert_offset, vert_offset + num_samples - 1)
        ]
        edge_list.extend(edges)

        sample_offset += num_samples

    # Correct return datatypes
    if flatten_cells:
        all_edges = np.array(all_edges)
    else:
        all_verts = [np.array(verts) for verts in all_verts]
    return all_verts, all_edges


def edges_to_PLY(vertices, edges, filepath, rgb=(0.0, 0.0, 0.0), text=False,
                 multiple=False):
    """
    Write edges defined as pairs of vertices to PLY file.

    @param  vertices : iterable[iterable[float]]
            Collection of 3D vertices.

    @param  edges : iterable[iterable[int]]
            Collection of edges index pairs into vertex collection

    @param  multiple : bool
            If True, treat arguments 'vertices' and 'edges' as multiple
            collections of vertices and edges. These will be concatenated into
            one PLY file  
    """
    from plyfile import PlyData, PlyElement

    # multiple = not isinstance(vertices[0][0], float)
    if multiple:
        sets_num_verts = [len(verts_set) for verts_set in vertices]
        sets_num_preceding = pre = [0] + list(np.cumsum(sets_num_verts))
        vertices = [tuple(v) for verts_set in vertices for v in verts_set]

    if not isinstance(vertices[0], tuple):
        vertices = [tuple(v) for v in vertices]

    # Create vertices as PLY elements
    vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
    vertices = np.array(vertices, dtype=vertex_dtype) # argument must be list of tuple
    verts_element = PlyElement.describe(vertices, 'vertex')
    
    # Create edges in PLY format
    if multiple:
        secs_edges = [
            (e[0] + pre[i_set], e[1] + pre[i_set], rgb[0], rgb[1], rgb[2]) 
                for i_set, edge_set in enumerate(edges) for e in edge_set
        ]
    else:
        secs_edges = [(e[0], e[1], rgb[0], rgb[1], rgb[2]) for e in edges]

    # Concatenate all edges
    edge_dtype = [ # edge format: http://paulbourke.net/dataformats/ply/
        ('vertex1', 'i4'),
        ('vertex2', 'i4'), 
        ('red',   'u1'),
        ('green', 'u1'),
        ('blue',  'u1')
    ]
    edges = np.array(secs_edges, dtype=edge_dtype)
    edges_element = PlyElement.describe(edges, 'edge')

    # Write vertices and faces to PLY file
    elements = [verts_element, edges_element]
    PlyData(elements, text=text).write(filepath)


def morphology_to_PLY(section_lists, filepath, segment_centers=True,
                      scale=1.0, rgb=(0.0, 0.0, 0.0), text=False,
                      translation=None, transform=None,
                      make_edges=True, make_faces=False):
    """
    Write neuron morphology to PLY file using degenerate faces.

    Each segment is represented by a degenerate face, i.e. a triangle
    with vertices [a b a].

    @pre    requires package plyfile, e.g. `pip install plyfile`

    @param  segment_centers : bool
            If true, write 3D locations of segment centers (nodes, i.e. centers
            of simulated compartments). This is useful for knowing the locations
            of compartments, and their current/voltage sources in 3D space.

    @param  scale : float
            Scale factor applied to coordinates after transform and translation
            is applied. Translation and transform not affected by scale.
    """
    from plyfile import PlyData, PlyElement

    # Get 3D samples
    if segment_centers:
        verts_data = morph_3d.get_segment_centers(section_lists, samples_as_rows=True)
    else:
        verts_data = morph_3d.get_section_samples(section_lists, include_diam=False)
    samples_xyz, secs_num3d, _ = verts_data

    # Apply transformation before writing
    if (translation is not None) or (transform is not None) or (scale != 1.0):
        samples_mat = np.ones((len(samples_xyz), 4))
        samples_mat[:,:3] = samples_xyz
        A = np.array(transform) if transform else np.eye(4)
        if (translation is not None):
            A[:-1, 3] += translation
        samples_mat = np.dot(samples_mat, A.T)
        if scale != 1.0:
            samples_mat *= scale
        samples_xyz = [tuple(row[:3]) for row in samples_mat]

    # Create vertices as PLY elements
    vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
    vertices = np.array(samples_xyz, dtype=vertex_dtype) # argument must be list of tuple
    verts_element = PlyElement.describe(vertices, 'vertex')
    
    # Create all faces and edges
    secs_faces = []
    secs_edges = []

    sample_offset = 0
    for i_sec, num_samples in enumerate(secs_num3d):

        # Face is degenerate face [a, b, a]
        faces = [
            ([i, i+1, i], rgb[0], rgb[1], rgb[2])
                for i in range(sample_offset, sample_offset + num_samples - 1)
        ]
        secs_faces.extend(faces)

        # Add edge [a, b]
        edges = [
            (i, i+1, rgb[0], rgb[1], rgb[2])
                for i in range(sample_offset, sample_offset + num_samples - 1)
        ]
        secs_edges.extend(edges)

        sample_offset += num_samples


    # Concatenate all faces
    face_dtype = [
        ('vertex_indices', 'i4', (3,)), 
        ('red',   'u1'),
        ('green', 'u1'),
        ('blue',  'u1')
    ]
    faces = np.array(secs_faces, dtype=face_dtype)
    faces_element = PlyElement.describe(faces, 'face')

    # Concatenate all edges
    edge_dtype = [ # edge format: http://paulbourke.net/dataformats/ply/
        ('vertex1', 'i4'),
        ('vertex2', 'i4'), 
        ('red',   'u1'),
        ('green', 'u1'),
        ('blue',  'u1')
    ]
    edges = np.array(secs_edges, dtype=edge_dtype)
    edges_element = PlyElement.describe(edges, 'edge')

    # Write vertices and faces to PLY file
    elements = [verts_element]
    if make_faces:
        elements.append(faces_element)
    if make_edges:
        elements.append(edges_element)
    PlyData(elements, text=text).write(filepath)


def SWC_nrn_to_PLY(morphology_path, ply_path=None, icell=None, quiet=False, 
                      **ply_kwargs):
    """
    Export SWC or ASC morphology to PLY file (easy to import in Blender
    as vertices + edges).

    Morphology is first instantiated in NEURON and then written to PLY file.

    @return     icell : object
                dummy cell object containing SectionLists of section types
                defined in the morphology file.
    """
    h.load_file('stdlib.hoc')
    h.load_file('import3d.hoc')

    extension = morphology_path.split('.')[-1]
    if extension.lower() == 'swc':
        imorphology = h.Import3d_SWC_read()
    elif extension.lower() == 'asc':
        imorphology = h.Import3d_Neurolucida3()
    else:
        raise ValueError("Unknown filetype: %s" % extension)

    if quiet:
        imorphology.quiet = 1
        h.hoc_stdout('/dev/null') # use NULL for Windows

    # Read morphology
    imorphology.input(str(morphology_path))

    # Instantiate in NEURON
    if icell is None:
        class ICell(object):
            pass
        icell = ICell
    importer = h.Import3d_GUI(imorphology, 0)
    importer.instantiate(icell)

    if quiet:
        h.hoc_stdout()

    # Convert to PLY
    if ply_path is None:
        ply_path = morphology_path[:-3] + 'ply'
    morphology_to_PLY([icell.all], ply_path, **ply_kwargs)
    print("Wrote morphology to %s" % ply_path)

    return icell


def SWC_raw_to_PLY(morphology_path, ply_path=None, single_ply_file=False,
                   scale=1.0, **ply_kwargs):
    """
    Export SWC or ASC morphology to PLY file (easy to import in Blender
    as vertices + edges).
    """
    # TODO: use structure ID to set edge color in PLY file
    if isinstance(morphology_path, (list, tuple)):
        morph_paths = morphology_path
    elif not morphology_path.endswith('.swc'):
        morph_paths = [os.path.join(morphology_path, f) for f in os.listdir(morphology_path) if f.endswith('.swc')]
    else:
        morph_paths = [morphology_path]

    all_cell_vertices = []
    all_cell_edges = []
    for morph_path in morph_paths:
        samples_list = read_SWC_samples(morph_path)
        vertices =  rfn.structured_to_unstructured(samples_list, dtype=float)[:, 2:5]
        if scale != 1.0:
            vertices *= scale
        edges = []
        for sample in samples_list:
            sample_id, structure_id, x, y, z, radius, parent_id = sample
            # sample_id and parent_id are 1-based indices
            if parent_id >= 0:
                edges.append([parent_id-1, sample_id-1])

        if not single_ply_file:
            out_filepath = os.path.join(ply_path,
                os.path.split(morph_path)[1][:-3] + 'ply')
            edges_to_PLY(vertices, edges, out_filepath, **ply_kwargs)
        else:
            all_cell_vertices.append(vertices)
            all_cell_edges.append(edges)

    if single_ply_file:
        edges_to_PLY(all_cell_vertices, all_cell_edges, ply_path, multiple=True,
                     **ply_kwargs)
            


def morphology_to_TXT(section_lists, filepath, segment_centers=True,
                      scale=1.0, rgb=(0.0, 0.0, 0.0),
                      translation=None, transform=None, precision_mm=1e-6):
    """
    Write neuron morphology to PLY file using degenerate faces.

    Each segment is represented by a degenerate face, i.e. a triangle
    with vertices [a b a].

    @pre    requires package plyfile, e.g. `pip install plyfile`

    @param  segment_centers : bool
            If true, write 3D locations of segment centers (nodes, i.e. centers
            of simulated compartments). This is useful for knowing the locations
            of compartments, and their current/voltage sources in 3D space.

    @param  scale : float
            Scale factor applied to coordinates after transform and translation
            is applied. Translation and transform not affected by scale.
    """

    # Get 3D samples
    if segment_centers:
        verts_data = morph_3d.get_segment_centers(section_lists, samples_as_rows=True)
    else:
        verts_data = morph_3d.get_section_samples(section_lists, include_diam=False)
    samples_xyz, secs_num3d, _ = verts_data

    # Apply transformation before writing
    if (translation is not None) or (transform is not None) or (scale != 1.0):
        samples_mat = np.ones((len(samples_xyz), 4))
        samples_mat[:,:3] = samples_xyz
        A = np.array(transform) if transform else np.eye(4)
        if (translation is not None):
            A[:-1, 3] += translation
        samples_mat = np.dot(samples_mat, A.T)
        if scale != 1.0:
            samples_mat *= scale
        samples_xyz = samples_mat[:, :3]
    else:
        samples_xyz = np.array(samples_xyz)

    # Set precision to 0.001 um
    num_significant = int(-np.log10(precision_mm) - 3 - np.log10(scale))
    fmt = '%.{:d}e'.format(num_significant)

    # Save to text file
    np.savetxt(filepath, samples_xyz, fmt=fmt)


def coordinates_um_to_txt(coords, filepath, precision_mm=1e-6, scale=1.0):
    """
    @param  scale : float
            Scale factor applied to coordinates after transform and translation
            is applied. Translation and transform not affected by scale.
    """
    # Set precision to 0.001 um
    num_significant = int(-np.log10(precision_mm) - 3 - np.log10(scale))
    fmt = '%.{:d}e'.format(num_significant)

    # Save to text file
    np.savetxt(filepath, coords * scale, fmt=fmt)


if __name__ == '__main__':
    test_json_export()