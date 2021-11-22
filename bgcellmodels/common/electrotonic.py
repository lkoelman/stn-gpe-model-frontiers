"""
Electrotonic analysis in NEURON

@author Lucas Koelman
"""

import math
from neuron import h

from . import treeutils

def lambda_DC(sec, gleak):
    """
    Compute electrotonic length constant of section in units of micron [um]
    """
    # Convert membrane resistance to same units as Ra
    # R_m = 1./(gleak*math.pi*sec.diam*1e-4) # r_m = R_m [Ohm*cm^2] /(pi*d) [Ohm*cm]
    R_m = 1./gleak # units [Ohm*cm^2]
    return 1e2 * math.sqrt(sec.diam*R_m/(4*sec.Ra)) # units: ([um]*[Ohm*cm^2]/(Ohm*cm))^1/2 = [um*1e2]


def lambda_AC(sec, f):
    """
    Compute electrotonic length constant (taken from stdlib.hoc)
    """
    return 1e5 * math.sqrt(sec.diam/(4*math.pi*f*sec.Ra*sec.cm))


def calc_lambda_AC(f, diam, Ra, cm):
    """
    Compute electrotonic length constant
    """
    return 1e5 * math.sqrt(diam/(4*math.pi*f*Ra*cm))


def calc_lambda(f, diam, Ra, gleak, cm):
    """
    Compute electrotonic length constant
    """
    if f <= 0:
        R_m = 1./gleak # units [Ohm*cm^2]
        return 1e2 * math.sqrt(diam*R_m/(4*Ra))
    else:
        return 1e5 * math.sqrt(diam/(4*math.pi*f*Ra*cm))


def sec_lambda(sec, gleak, f):
    """
    Compute electrotonic length constant at given frequency
    """
    if f <= 0:
        return lambda_DC(sec, gleak)
    else:
        return lambda_AC(sec, f)


def seg_lambda(seg, gleak, f):
    """
    Compute length constant of segment
    """
    Ra = seg.sec.Ra # Ra is section property
    if f <= 0:
        if isinstance(gleak, str):
            Rm = 1./getattr(seg, gleak)
        else:
            Rm = 1./gleak # units [Ohm*cm^2]
        return 1e2 * math.sqrt(seg.diam*Rm/(4*Ra)) # units: ([um]*[Ohm*cm^2]/(Ohm*cm))^1/2 = [um*1e2]
    else:
        return 1e5 * math.sqrt(seg.diam/(4*math.pi*f*Ra*seg.cm))


def seg_L_elec(seg, gleak, f):
    """
    Electrotonic length of segment
    """
    return (seg.sec.L/seg.sec.nseg)/seg_lambda(seg, gleak, f)


def min_nseg_hines(sec, f=100.):
    """
    Minimum number of segments based on electrotonic length
    """
    return int(sec.L/(0.1*lambda_AC(sec, f))) + 1


def min_nseg_marasco(sec):
    """
    Minimum number of segments based on electrotonic length
    """
    return int((sec.L/(0.1*lambda_AC(sec,100.))+0.9)/2)*2 + 1  


def calc_min_nseg_hines(f, L, diam, Ra, cm, round_up=True):
    lamb_AC = 1e5 * math.sqrt(diam/(4*math.pi*f*Ra*cm))
    nseg_trunc = int(L/(0.1*lamb_AC)) # rounded down
    if round_up:
        return nseg_trunc + 1
    else:
        return max(nseg_trunc, 1)


def set_min_nseg_hines(seclist, f_lambda, round_up=True,
                       add=True, remove=False):
    """
    Set number of segments based on Hines' rule of thumb.

    NOTE: uses diam and cm in the middle segment (x = 0.5)

    @return     The number of additionally created segments.
    """
    extra_nseg = 0
    for sec in seclist:
        min_nseg = calc_min_nseg_hines(
                        f_lambda, sec.L, sec.diam, sec.Ra, sec.cm,
                        round_up=False)
        if add and min_nseg > sec.nseg:
            extra_nseg += min_nseg - sec.nseg
            sec.nseg = min_nseg
        if remove and min_nseg < sec.nseg:
            extra_nseg += min_nseg - sec.nseg
            sec.nseg = min_nseg
    return extra_nseg


def inputresistance_inf(sec, gleak, f):
    """
    Input resistance for semi-infinite cable in units of [Ohm*1e6]
    """
    lamb = sec_lambda(sec, gleak, f)
    R_m = 1./gleak # units [Ohm*cm^2]
    return 1e2 * R_m/(math.pi*sec.diam*lamb) # units: [Ohm*cm^2]/[um^2] = [Ohm*1e8]


def inputresistance_sealed(sec, gleak, f):
    """
    Input resistance of finite cable with sealed end in units of [Ohm*1e6]
    """
    x = sec.L/sec_lambda(sec, gleak, f)
    return inputresistance_inf(sec, gleak, f) * (math.cosh(x)/math.sinh(x))


def inputresistance_leaky(sec, gleak, f, R_end):
    """
    Input resistance of finite cable with leaky end in units of [Ohm*1e6]
    
    @param R_end    input resistance of connected cables at end of section
                    in units of [Ohm*1e6]
    """
    R_inf = inputresistance_inf(sec, gleak, f)
    x = sec.L/sec_lambda(sec, gleak, f)
    return R_inf * (R_end/R_inf*math.cosh(x) + math.sinh(x)) / (R_end/R_inf*math.sinh(x) + math.cosh(x))


def inputresistance_tree(rootsec, f, glname):
    """
    Compute input resistance to branching tree
    """

    childsecs = rootsec.children()
    gleak = getattr(rootsec, glname)

    # Handle leaf sections
    if not any(childsecs):
        return inputresistance_sealed(rootsec, gleak, f)

    # Calc input conductance of children
    g_end = 0.
    for childsec in childsecs:
        g_end += 1./inputresistance_tree(childsec, f, glname)
    return inputresistance_leaky(rootsec, gleak, f, 1./g_end)


def measure_current_transfer(
        source=None,
        target=None,
        imp=None,
        freq=None,
        linearize_gating=None,
        precompute=False,
    ):
    """
    Measure current attenuation from source segment to target segment.

    @param  imp : h.Impedance
            Impedance measuring prove

    @param  freq : float
            Frequency to calculate measure at

    @param  source : nrn.Segment
            source segment

    @param  target : nrn.Segment
            target segment

    @return A_{I,source->target} : float

            A_{I,source->target} = (Itarget/Isource)
                                 = (Vsource/Varget | Isource=0)
    """
    # NOTE: Impedance.ratio(x, sec=sec) measures |v(impedance.loc)/v(sec.x)|
    #       == A_{V,sec.x->Impedance.loc}
    #       == A_{I,Impedance.loc->sec.x}
    if imp is None:
        imp = h.Impedance()
    if not precompute:
        imp.loc(source.x, sec=source.sec)
        imp.compute(freq, int(linearize_gating))
    return imp.ratio(target.x, sec=target.sec)


def measure_voltage_transfer(
        source=None,
        target=None,
        imp=None,
        freq=None,
        linearize_gating=None,
        precompute=False
    ):
    """
    Measure voltage attenuation from source segment to target segment.

    @return A_{V,source->target} : float

            A_{V,source->target} = (Vtarget/Vsource | Itarget=0)
                                 = (Isource/Itarget)
                                 = A_{I,target->source}

    @see    measure_Ai

    @note   this is the same as a call to with switched arguments:
            measure_current_transfer(source=target, target=source)
    """
    if imp is None:
        imp = h.Impedance()
    if not precompute:
        imp.loc(target.x, sec=target.sec)
        imp.compute(freq, int(linearize_gating))
    return imp.ratio(source.x, sec=source.sec)


def measure_transfer_impedance(
        source=None,
        target=None,
        imp=None,
        freq=None,
        linearize_gating=None,
        precompute=False
    ):
    """
    Measure transfer impedance from source segment to target segment,
    
    @return Zc : float

            Zc = (Vtarget/Isource | Itarget=0) = (Vsource/Itarget | Isource=0)

    @see    measure_Ai
    """
    if imp is None:
        imp = h.Impedance()
    if not precompute:
        imp.loc(target.x, sec=target.sec)
        imp.compute(freq, int(linearize_gating))
    return imp.transfer(source.x, sec=source.sec)


def measure_input_impedance(
        source=None,
        target=None,
        imp=None,
        freq=None,
        linearize_gating=None,
        precompute=False
    ):
    """
    Measure input impedance at source segment.
    
    @return Zin : float

            Zin = Vsource/Isource

    @see    measure_Ai

    @note   argument 'target' is unused, but present to provide same function
            signature as other electrotonic measurement functions
    """
    if imp is None:
        imp = h.Impedance()
    if not precompute:
        imp.loc(source.x, sec=source.sec)
        imp.compute(freq, int(linearize_gating))
    return imp.input(source.x, sec=source.sec)


def measure_along_paths(root, leaves, measure_funcs, freq):
    """
    Measure electrotonic properties along paths.
    
    @param  root : nrn.Section
            Root section for measurement

    @param  leaves : iterable[nrn.Section]
            Endpoints for measurements

    @param  measure_funcs : callable
            Any of the functions measure_X from this module or function
            with the same signature.
    """
    # make paths to leaves and measure each measure
    probe = h.Impedance()
    leaf_distance_measures = [] # one dict for each leaf

    # Get path lengths
    h.distance(0, 0.5, sec=root)
    pathlen_func = lambda seg: h.distance(seg.x, sec=seg.sec)

    # Get electrotonic measurements
    for leaf_sec in leaves:

        dist_measures_vecs = {}

        path_sections = treeutils.path_sections(leaf_sec)
        path_segments = [seg for sec in path_sections for seg in sec]

        # Get path lengths
        dist_measures_vecs['pathlen_micron'] = map(pathlen_func, path_segments)

        # do each measurement
        for measure, dist_func in measure_funcs.items():
            measure_func = lambda seg: dist_func(
                                        source=seg,
                                        target=root(0.5),
                                        imp=probe,
                                        freq=freq,
                                        linearize_gating=0)
            distance_vec = map(measure_func, path_segments)
            dist_measures_vecs[measure] = distance_vec

        leaf_distance_measures.append(dist_measures_vecs)
        
    return leaf_distance_measures


def segs_at_dX(cur_seg, dX, f, gleak):
    """
    Get next segments (in direction 0 to 1 end) at electrotonic distance dX
    from current segment (dX = L/lambda).

    I.e. get the segment at X(cur_seg) + dX

    Assumes entire section has same diameter
    """
    nseg = cur_seg.sec.nseg
    sec_L = cur_seg.sec.L
    sec_dL = sec_L / nseg
    sec_lamb = [seg_lambda(seg, gleak, f) for seg in cur_seg.sec]
    sec_Xi = [sec_dL / sec_lamb[i] for i in range(nseg)]
    sec_Xtot = sum(sec_Xi) # total L/lambda
    sec_Xacc = [sum((sec_Xi[j] for j in range(i)), 0.0) for i in range(nseg)] # accumulated X to left border
    sec_dx = 1.0/nseg

    # Get electrotonic length from start of Section to current x
    cur_i = min(int(cur_seg.x/sec_dx), nseg-1)
    l_x = cur_i*sec_dx      # x of left boundary
    l_X = sec_Xacc[cur_i]   # L/lambda of left boundary

    frac = (cur_seg.x - l_x) / sec_dx
    assert 0.0 <= frac <= 1.01
    frac = min(1.0, frac)
    cur_X = l_X + frac*sec_Xi[i]

    # Get length to next discretization step
    bX = cur_X + dX
    next_segs = []

    if bX <= sec_Xtot:

        # Find segment that X falls in
        i_b = next((i for i, Xsum in enumerate(sec_Xacc) if Xsum+sec_Xi[i] >= bX))
        
        # Interpolate
        frac = (bX - sec_Xacc[i_b]) / sec_Xi[i_b]
        assert 0.0 <= frac <= 1.01
        frac = min(1.0, frac)
        
        b_x = i_b*sec_dx + frac*sec_dx
        next_segs.append(cur_seg.sec(b_x))

    else:
        ddX = bX - sec_Xtot # how far to advance in next section
        next_segs = []
        for child_sec in cur_seg.sec.children():
            next_segs.extend(segs_at_dX(child_sec(0.0), ddX, f, gleak))

    return next_segs