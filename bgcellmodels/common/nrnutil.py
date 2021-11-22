"""
Utilities for dealing with NEURON

@author Lucas Koelman
"""

import os

import neuron
h = neuron.h

from .nrnmodutil import get_mod_name # for backwards compatibility

# default NEURON section arrays and SectionList names for cell prototypes
# created by Hoc Import3d tool (for loading 3D morphologies)
nrn_proto_seclists_arrays = {
    'all'       : None,
    'somatic'   : 'soma',
    'basal'     : 'dend',
    'apical'    : 'apic',
    'axonal'    : 'axon',
    'myelinated': 'myelin',
}


class ExtSection(neuron.hclass(h.Section)):
    """ Extension of Section to allow modifying attributes """
    pass


class ExtSecRef(neuron.hclass(h.SectionRef)):
    """ Extension of SectionRef to allow modifying attributes """
    pass


class ICell(object):
    """
    Cell object with section array and SectionList names
    conforming to Hoc prototype of Import3D morphology importer.
    """
    def __init__(self, **kwargs):
        """
        Make new ICell

        @param  kwargs : dict[str, <list[Section] or h.SectionList>]
                Section arrays or section lists according to Import3d prototype
        """
        # NOTE: since SectionList do not keep references alive, sections
        #       must be stored in secarray or python list
        self._all = []
        self.all = h.SectionList()
        if 'all' in kwargs:
            for sec in kwargs['all']:
                self._all.append(sec)
                self.all.append(sec=sec)

        for seclist_name, array_name in nrn_proto_seclists_arrays.items():
            if seclist_name == 'all':
                continue

            sl = h.SectionList()
            setattr(self, seclist_name, sl)
            source = None
            
            # Check if sections were provided
            if array_name in kwargs:
                source = list(kwargs[array_name])
            elif seclist_name in kwargs:
                source = list(kwargs[seclist_name])
            
            # Fill section array and SectionList
            if source is None:
                setattr(self, array_name, [])
            else:
                setattr(self, array_name, source)
                for sec in source:
                    sl.append(sec=sec)

            if 'all' not in kwargs:
                for sec in getattr(self, array_name):
                    self._all.append(sec)
                    self.all.append(sec=sec)


def hoc_load_from(working_dir, script):
    """
    Load Hoc code with given working directory.

    NOTE: alternative is to append dirrectoires to the Hoc library path
          accessible via environment variable os.environ['HOC_LIBRARY_PATH']
    """
    prev_cwd = os.getcwd()
    os.chdir(working_dir)
    h.xopen(script)
    os.chdir(prev_cwd)


# Create equivalent section
def create_hoc_section(secname):
    """
    Create Section with given name on global Hoc object

    @post   Section will be available as attribute on global Hoc object (neuron.h).
            This will ensure that the Section is not destroyed, even though
            no Python reference to it exists.

    @return         tuple(Section, SectionRef)
    """

    if secname in [sec.name() for sec in h.allsec()]:
        raise Exception('Section named {} already exists'.format(secname))

    created = h("create %s" % secname)
    if created != 1:
        raise Exception("Could not create section with name '{}'".format(secname))

    eqsec = getattr(h, secname)
    eqref = ExtSecRef(sec=eqsec)
    return eqsec, eqref


def join_seclists(*args):
    """
    Join together all Hoc SectionList into a single SectionList.
    """
    all_sections = h.SectionList()
    for seclist in args:
        for sec in seclist:
            all_sections.append(sec=sec)
    return all_sections


def getsecref(sec, refs):
    """
    Look for SectionRef pointing to Section sec in enumerable of SectionRef

    @return     <SectionRef/NoneType>
                first SectionRef in refs with same section name as sec
    """
    if sec is None: return None
    # Section names are unique, but alternatively use sec.same(ref.sec)
    return next((ref for ref in refs if (ref.exists() and ref.sec.name()==sec.name())), None)


def contains_sec(seclist, sec):
    """
    Check if enumerable contains given section
    """
    return any([sec_b.same(sec) for sec_b in seclist])


def seg_index(tar_seg):
    """
    Get index of given segment on Section
    """
    seg_dx = 1.0/tar_seg.sec.nseg
    seg_id = int(tar_seg.x // seg_dx)
    return min(seg_id, tar_seg.sec.nseg-1)

    # NOTE: == operator compares actual segments
    # if tar_seg.x == 0.0: # start node
    #     return 0
    # elif tar_seg.x == 1.0: # end node
    #     return tar_seg.sec.nseg-1
    # else: # internal node
    #     return next(i for i,seg in enumerate(tar_seg.sec) if seg==tar_seg)


def seg_at_index(sec, iseg):
    """
    Get the i-th segment of Section

    @return     nrn.Segment with x-value equal to segment midpoint
    """
    xmid = (2.*(iseg+1)-1.)/(2.*sec.nseg) # See NEURON book p. 8
    # return next(seg for i,seg in enumerate(sec) if i==iseg)
    return sec(xmid)


def seg_at_node(sec, inode):
    """
    Get the segment at given node index.

    A section's nodes are the centers of its segments and the 0 and 1 ends.
    """
    if inode == 0:
        return sec(0.0)
    elif inode == sec.nseg + 1:
        return sec(1.0)
    else:
        return sec((2.*(inode)-1.)/(2.*sec.nseg))


def all_xnode(sec):
    """
    Normalized locations of nodes, which are the segment centers
    and the nodes at the 0-end and 1-end of the section.

    @return   list[float]
    """
    nseg = sec.nseg
    node_locs = h.Vector(nseg + 2)
    node_locs.indgen(1.0 / nseg)
    node_locs.sub(1.0 / (2 * nseg))
    node_locs.x[0] = 0.0
    node_locs.x[nseg+1] = 1.0
    return list(node_locs)


def seg_xmid(seg):
    """
    x-value at segment midpoint
    """
    nseg = seg.sec.nseg
    iseg = seg_index(seg)
    xmid = (2.*(iseg+1)-1.)/(2.*nseg) # See NEURON book p. 8
    return xmid


def seg_xmin(seg, side=None):
    """
    x-value at left boundary of segment (towards 0-end)

    @param  side : str
            Relative location of return x-value to the exact segment boundary
                - 'inside' : inside given segment
                - 'outside': in previous segment or 0-end node
                - 'boundary' or None: exactly on segment boundary, no guarantee
                   whether this is inside or outside given segment
    """
    nseg = seg.sec.nseg
    iseg = seg_index(seg)
    x_lo = (1.0/nseg) * iseg

    if side=='inside':
        x_lo += 1e-12
    elif side=='outside':
        x_lo -= 1e-12
    elif (side=='boundary') or (side is None):
        return x_lo
    else:
        raise ValueError(side)

    return max(0.0, x_lo)

seg_xleft = seg_xmin


def seg_xmax(seg, side=None):
    """
    x-value at right boundary of segment (towards 1-end)

    @param  side : str
            Relative location of return x-value to the exact segment boundary
                - 'inside' : inside given segment
                - 'outside': in previous segment or 0-end node
                - 'boundary' or None: exactly on segment boundary, no guarantee
                   whether this is inside or outside given segment
    """
    nseg = seg.sec.nseg
    iseg = seg_index(seg)
    x_hi = (1.0/nseg) * (iseg + 1)

    if side=='inside':
        x_hi -= 1e-12
    elif side=='outside':
        x_hi += 1e-12
    elif (side=='boundary') or (side is None):
        return x_hi
    else:
        raise ValueError(side)

    return min(1.0, x_hi)

seg_xright = seg_xmax


def same_seg(seg_a, seg_b):
    """
    Check whether both segments are in same Section and their x location
    maps to the same segment.
    """
    if not(seg_a.sec.same(seg_b.sec)):
        return False
    # check if x locs map to same section
    return seg_index(seg_a) == seg_index(seg_b)


def get_range_var(seg, varname, default=0.0):
    """
    Get RANGE variable at segment regardless whether the mechanism exists
    or not.
    """
    try:
        return getattr(seg, varname, default)
    except NameError: # mechanisms is known but not inserted
        return default


def copy_ion_styles(src_sec, tar_sec, ions=None):
    """
    Copy ion styles from source to target section

    NOTE:

    oldstyle = ion_style("name_ion")

    oldstyle = int:
        int( 1*c_style + 4*cinit + 8*e_style + 32*einit + 64*eadvance )
        c_style:    0, 1, 2, 3  (2 bits)
        e_style:    0, 1, 2, 3  (2 bits)
        einit:      0, 1        (1 bits)
        eadvance:   0, 1        (1 bits)
        cinit:      0, 1        (1 bits)

    ion_style("name_ion", c_style, e_style, einit, eadvance, cinit)

    """
    if ions is None:
        ions = ['na', 'k', 'ca']

    # Get ion style for each ion species
    src_sec.push()
    styles = dict(((ion, h.ion_style(ion+'_ion')) for ion in ions))

    # Copy to target Section
    set_ion_styles(tar_sec, **styles)

    h.pop_section()


def get_ion_styles(src_sec, ions=None):
    """
    Get ion styles as integer for each ion.

    @return     dict[str, int] mapping ion species to integer containing
                the bit flags that indicate ion styles
    """
    if not isinstance(ions, list):
        ions = ['na', 'k', 'ca']

    # Get ion style for each ion species
    styles = {}
    for ion in ions:
        ion_var = ion + '_ion'
        if hasattr(src_sec(0.5), ion_var):
            styles[ion] = h.ion_style(ion_var, sec=src_sec)

    return styles


def set_ion_styles(tar_sec, **kwargs):
    """
    Set ion styles from integer containing bit flags.

    @param  kwargs      keyword arguments ion_name: style_int
    """
    for ion, style in kwargs.items():

        # Decompose int into bit flags
        c_style = int(style) & (1+2)
        cinit = (int(style) & 4) >> 2
        e_style = (int(style) & (8+16)) >> 3
        einit = (int(style) & 32) >> 5
        eadvance = (int(style) & 64) >> 6

        # Copy to new section
        h.ion_style(ion+'_ion', c_style, e_style, einit, eadvance, cinit, sec=tar_sec)


def make_ion_style_flags(style):
    """
    Decompose bit flags containing ion styles into individual
    flags

    EXAMPLE
    -------

    >>> styles = h.ion_style('na_ion', sec=sec)
    >>> flags = make_ion_style_flags(style)
    >>> h.ion_style('na_ion', *flags, sec=other)
    """
    c_style = int(style) & (1+2)
    cinit = (int(style) & 4) >> 2
    e_style = (int(style) & (8+16)) >> 3
    einit = (int(style) & 32) >> 5
    eadvance = (int(style) & 64) >> 6

    return c_style, e_style, einit, eadvance, cinit


def ion_styles_bits_to_dict(style):
    """
    Convert a float representing the styles of one ion to a dictionary
    containing values for all the style flags.

    @param  style : float
            Result of call to h.ion_style(ion) for the CAS

    @return styles: dict
            Names of styles flags and their values
    """
    # Decompose int into bit flags
    styles = {}
    styles['c_style'] = int(style) & (1+2)
    styles['cinit'] = (int(style) & 4) >> 2
    styles['e_style'] = (int(style) & (8+16)) >> 3
    styles['einit'] = (int(style) & 32) >> 5
    styles['eadvance'] = (int(style) & 64) >> 6

    return styles


def test_segment_boundaries():
    """
    Test functions related to segment x-values
    """
    from stdutil import isclose

    for nseg in range(1,11):
        sec = h.Section()
        sec.nseg = nseg
        dx = 1.0 / nseg

        # Assign true index to segment property
        sec.insert('hh')
        for i, seg in enumerate(sec):
            seg.gnabar_hh = float(i)

        # Test each segment
        for i, seg in enumerate(sec):

            # test midpoint of segment
            x_mid = seg_xmid(seg)
            assert isclose(seg.x, x_mid, abs_tol=1e-9)
            assert seg_index(seg) == i

            # test low boundary of segment
            x_min = seg_xmin(seg, side='inside')
            lo_seg = sec(x_min)
            assert seg_index(lo_seg) == i
            assert lo_seg.gnabar_hh == float(i)
            assert isclose(x_min, x_mid-dx/2, abs_tol=1e-9)

            if i != 0:
                x_min = seg_xmin(seg, side='outside')
                prev_seg = sec(x_min)
                assert seg_index(prev_seg) == i-1
                assert prev_seg.gnabar_hh == float(i-1)
                assert isclose(x_min, x_mid-dx/2, abs_tol=1e-9)

            # test high boundary of segment
            x_max = seg_xmax(seg, side='inside')
            hi_seg = sec(x_max)
            assert seg_index(hi_seg) == i
            assert hi_seg.gnabar_hh == float(i)
            assert isclose(x_max, x_mid+dx/2, abs_tol=1e-9)

            if i != nseg-1:
                x_max = seg_xmax(seg, side='outside')
                next_seg = sec(x_max)
                assert seg_index(next_seg) == i+1
                assert next_seg.gnabar_hh == float(i+1)
                assert isclose(x_max, x_mid+dx/2, abs_tol=1e-9)


    print("Test 'test_segment_boundaries' passed.")


def independent_random_stream(stream_len, used_indices,
                              force_low_index=None,
                              start_low_index=0):
    """
    Return a Hoc Random object that will generate a stream of N
    independent random numbers, provided a list of seeds/indices that
    are already in use.

    @see        Documentation of MCellRan4 at following page:
                https://neuron.yale.edu/neuron/static/new_doc/programming/math/random.html?#Random.MCellRan4

    @note       Note that the distribution generated by the RNG has to be set
                manually by the caller, on the returned RNG object. E.g. call
                `rng.uniform(0,1)`, `rng.negexp(1)`, `rng.random(0,1)`, ...

    Arguments
    ---------

    @param      stream_len : int
                Desired length of stream.

    @param      start_low_index : int
                Only use low_index values starting from this value.

    @param      force_low_index : int
                Force the use of this low_index value. ValueError is raised
                if there are less than stream_len increments left.

    @param      used_indices : dict[int, int]
                The first integers are equal to the low_index values
                already in use. The second integers are the highest high_index
                value in use for that low_index.

    Returns
    -------

    @return     (RNG, init_rng) : tuple[Hoc.Random, callable]
                The first object is the Hoc Random object and the second
                is an initialization function to initialize the stream.

    @post       The used_indices dict is modified to reflect the new indices
                in use.


    """
    # Look for a low_index that still has stream_len increments left in
    # its high_index
    low_index, high_index = None, None
    valid_low_indices = sorted([k for k in used_indices.keys() if k >= start_low_index])

    if len(valid_low_indices) == 0:
        # Case 1: no indices in use
        low_index = start_low_index
        high_index = stream_len
        used_indices[low_index] = high_index

    elif force_low_index is not None:
        # Case 2: caller specified the low_index
        low_index = force_low_index
        if force_low_index not in used_indices:
            high_index = stream_len
            used_indices[low_index] = high_index
        else:
            hi = used_indices[low_index]
            if (2**32 - 1 - hi) < stream_len:
                raise ValueError(
                    "Cannot generate stream of length {} using "
                    "given low_index {}: insufficient increments left.".format(
                        stream_len, force_low_index))
            else:
                high_index = hi + stream_len
                used_indices[low_index] = high_index
    else:
        # Case 3: find suitable low and high index in dict
        for lo in valid_low_indices:
            hi = used_indices[lo]
            if (2**32 - 1 - hi) >= stream_len:
                low_index = lo
                high_index = hi + stream_len
        if low_index is not None:
            # Found a low index with stream_len increments left
            used_indices[low_index] = high_index
        else:
            # No low index found that has stream_len increments left
            low_index = valid_low_indices[-1] + 1
            high_index = stream_len
            used_indices[low_index] = high_index

    rng = h.Random()
    rng.MCellRan4(high_index, low_index)

    def start_stream():
        # captures rng and high_index from local scope
        rng.seq(high_index)

    return rng, start_stream
