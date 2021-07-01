"""
Algorithms for mapping synaptic input locations from detailed compartmental neuron model
to reduced morphology model.
"""

import re
import functools
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
logger = logging.getLogger('redops') # create logger for this module

from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

from .. import redutils
from bgcellmodels.common import configutil
from bgcellmodels.common.nrnutil import getsecref, seg_index
from bgcellmodels.common.treeutils import subtree_has_node
from bgcellmodels.common.stdutil import isclose
from bgcellmodels.common.electrotonic import (
    measure_voltage_transfer, measure_current_transfer, 
    measure_transfer_impedance, measure_input_impedance
)

class SynInfo(object):
    """
    Struct-like object containing synapse properties.
    

    It may contain following attributes:

        mod_name            str
        
        sec_name            str
        
        sec_hname           str
        
        sec_loc             float

        'mech_attr_i'       float
                            one attribute for each synaptic mechanism parameter

        afferent_netcons    list(NetCon)
        
        afferent_weights    list(list(float))
                            weight vector for each incoming NetCon

        'secref_attr_i'     object
                            saved SectionRef attributes


        path_ri             float
                            summed seg.ri() up to synapse segment
        
        max_path_ri         float
                            = max(syn_sec.pathri_seg) # max path resistance in Section
        
        min_path_ri         float
                            = min(syn_sec.pathri_seg) # min path resistance in Section

        Zc                  float
                            Ztransfer, i.e. |v(soma)/i(syn)| or |v(syn)/i(soma)|
        
        Zin                 float
                            Zinput, i.e. v(x)/i(x) measured at synapse

        Av_syn_to_soma          float
                            voltage transfer ratio, i.e. |v(soma)/v(syn|

        mapped_syn          HocObject
                            the new synapse that the original was mapped to
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    # def __repr__(self):
    #   return "{} at {}({})".format(self.mod_name, self.sec_hname, self.sec_loc)


# Parameter names for synaptic mechanisms in .MOD files
synmech_parnames = {
    'ExpSyn': ['tau', 'e'], # gmax stored in NetCon weight vector
    
    'Exp2Syn': ['tau1', 'tau2', 'e'], # gmax stored in NetCon weight vector
    
    'AlphaSynapse': ['onset', 'tau', 'gmax', 'e'],
    
    'GABAsyn': ['tau_r_GABAA', 'tau_d_GABAA', 'tau_r_GABAB', 'tau_d_GABAB', 
                'gmax_GABAA', 'gmax_GABAB', 'Erev_GABAA', 'Erev_GABAB', 
                'tau_rec', 'tau_facil', 'U1', 'use_stdp_A', 'use_stdp_B'],
    
    'GLUsyn': ['tau_r_AMPA', 'tau_d_AMPA', 'tau_r_NMDA', 'tau_d_NMDA', 'mg', 
                'gmax_AMPA', 'gmax_NMDA', 'tau_rec', 'tau_facil', 'U1', 'e'],
}




def measure_root_distance_micron(
        source=None,
        target=None,
        imp=None,
        freq=None,
        linearize_gating=None,
        precompute=False
    ):
    """
    Measure path distance to the root Section in units of micron [m*1e-6].
    """
    return redutils.seg_path_L(source, endpoint='segment_x')


electrotonic_measurement_funcs = {
    'Av_syn_to_soma' : measure_voltage_transfer,
    'Ai_syn_to_soma' : measure_current_transfer,
    'Ztransfer' : measure_transfer_impedance,
    'Zc' : measure_transfer_impedance,
    'root_distance_micron' : measure_root_distance_micron,
}


def get_syn_mapping_info(
        rootsec,
        allsecrefs,
        syn_mod_pars=None,
        Z_freq=25.,
        gleak_name=None,
        linearize_gating=False,
        init_cell=None,
        save_ref_attrs=None,
        sever_netcons=True,
        attr_mappers=None,
        nc_tomap=None,
        syn_tomap=None,
        electrotonic_measures=None,
    ):
    """
    For each synapse on the neuron, calculate and save information for placing an equivalent
    synaptic input on a morphologically simplified neuron.

    
    ARGUMENTS

        @param rootsec          any section of the cell

        @param syn_mod_pars     dict of <synaptic mechanism name> : list(<attribute names>)
                                containing attributes that need to be stored

        @param init_cell        function to bring the cell to the desired state to measure transfer
                                impedances, e.g. simulating under a particular input

        @param sever_netcons    Set NetCon target to None for each connection, this is useful
                                when the cell or synapses are destroyed


    RETURN VALUES

        @return                 list(SynInfo) containing information about each synapse
                                to be mapped


    USAGE

        - either provide syn_tomap (list(SynInfo)), or nc_tomap, or neither.
          In the latter case, all NetCon targetting the cell are used.
    """
    if syn_mod_pars is None:
        syn_mod_pars = synmech_parnames

    if save_ref_attrs is None:
        save_ref_attrs = []

    if attr_mappers is None:
        attr_mappers = {}

    # Calculate section path properties for entire tree
    for secref in allsecrefs:
        redutils.sec_path_props(secref, 100., gleak_name)

    # Measure transfer impedance and filter parameters
    logger.debug("Initializing cell for electrotonic measurements...")
    init_cell()
    
    imp = h.Impedance() # imp = h.zz # previously created
    last_freq = None


    # Find all Synapses on cell (all Sections in whole tree)
    if syn_tomap is not None:
        cell_syns = syn_tomap

    else:
        # Get references to all synapses ourselves
        if nc_tomap is None:
            # Get all NetCon targetting the cell
            dummy_syn = h.Exp2Syn(rootsec(0.5))
            dummy_nc = h.NetCon(None, dummy_syn)

            # unique synapses targeting the same cell
            cell_ncs = [nc for nc in list(dummy_nc.postcelllist()) if not nc.same(dummy_nc)] # all NetCon targeting same tree as dummy
            syns = set([nc.syn() for nc in cell_ncs]) # NOTE: this works, since set uses == for uniqueness, which compares using syn.hocobjptr()
        
        else:
            cell_ncs = nc_tomap
            syns = set([nc.syn() for nc in cell_ncs]) # NOTE: this works, since set uses == for comparison, which uses syn.hocobjptr()

        # Make a list of SynInfo ourselves
        cell_syns = [SynInfo(orig_syn=syn) for syn in syns]
        for syn in cell_syns:
            syn.PSP_median_frequency = Z_freq
            syn.afferent_netcons = [nc for nc in cell_ncs if nc.syn().same(syn.orig_syn)]

    if len(cell_syns) == 0:
        logger.warn("No synapses found on tree of Section {}".format(rootsec))

    # Save synapse properties
    logger.debug("Getting synapse properties...")
    for syn_info in cell_syns:

        # get synapse HocObject
        syn = syn_info.orig_syn
        
        # Get synaptic mechanism name
        match_mechname = re.search(r'^[a-zA-Z0-9]+', syn.hname())
        synmech = match_mechname.group()
        if synmech not in syn_mod_pars:
            raise Exception("Synaptic mechanism '{}' not in given mechanism list".format(synmech))

        # Get its segment and location
        syn_seg = syn.get_segment()
        syn_sec = syn_seg.sec
        syn_secref = getsecref(syn_sec, allsecrefs)
        syn_loc = syn.get_loc() # changes CAS
        h.pop_section()

        # Create struct to save synapse information
        syn_info.mod_name = synmech
        syn_info.sec_name = syn_sec.name()
        syn_info.sec_hname = syn_sec.hname()
        syn_info.sec_loc = syn_loc # can also use nc.postcell() and nc.postloc()

        # Save synaptic mechanism parameters
        mech_params = syn_mod_pars[synmech]
        for par in mech_params:
            setattr(syn_info, par, getattr(syn, par))

        # Save all NetCon objects targetting this synapse
        syn_info.afferent_weights = [[nc.weight[i] for i in xrange(int(nc.wcnt()))] for nc in syn_info.afferent_netcons]

        # Save requested properties of synapse SectionRef
        syn_info.saved_ref_attrs = save_ref_attrs
        for attr in save_ref_attrs:
            if hasattr(syn_secref, attr):
                setattr(syn_info, attr, getattr(syn_secref, attr))

        # Save other computed properties
        for attr, mapper in attr_mappers.iteritems():
            setattr(syn_info, attr, mapper(syn))

        # Get axial path resistance to synapse
        syn_info.path_ri = syn_secref.pathri_seg[seg_index(syn_seg)] # summed seg.ri() up to synapse segment
        syn_info.max_path_ri = max(syn_secref.pathri_seg) # max path resistance in Section
        syn_info.min_path_ri = min(syn_secref.pathri_seg) # min path resistance in Section
        syn_info.root_distance_micron = redutils.seg_path_L(syn_seg, endpoint='segment_x')
        syn_info.path_L = syn_info.root_distance_micron

        # Get electrotonic measures

        ## Get measures for sec.x (synapse) -> impedance.loc (soma)
        imp.loc(0.5, sec=rootsec) # injection site
        syn_freq = syn_info.PSP_median_frequency
        if (last_freq is None) or (syn_freq != last_freq):
            imp.compute(syn_freq, int(linearize_gating))

        syn_info.Zc = imp.transfer(syn_loc, sec=syn_sec) # query transfer impedanc,e i.e.  |v(loc)/i(x)| or |v(x)/i(loc)|
        syn_info.Ztransfer = syn_info.Zc
        syn_info.Zin = imp.input(syn_loc, sec=syn_sec) # query input impedance, i.e. v(x)/i(x)
        syn_info.Av_syn_to_soma = imp.ratio(syn_loc, sec=syn_sec) # A_{V,syn->soma} = |v(impedance.loc)/v(sec_x)|
        syn_info.Ai_soma_to_syn = syn_info.Av_syn_to_soma # A_{I,soma->syn}

        ## Get measures for impedance.loc (soma) -> sec.x (synapse)
        imp.loc(syn_loc, sec=syn_sec) # injection site
        imp.compute(syn_freq, int(linearize_gating))
        syn_info.Av_soma_to_syn = imp.ratio(0.5, sec=rootsec) # A_{V,soma->syn}
        syn_info.Ai_syn_to_soma = syn_info.Av_soma_to_syn # A_{I,syn->oma}

        # Point all afferent NetCon connections to nothing
        if sever_netcons:
            for nc in syn_info.afferent_netcons:
                nc.setpost(None)

        # Delete reference to synapse about to be deleted
        syn_info.orig_syn = None

    return cell_syns


def syn_was_in_absorbed(syn_info, secref):
    """
    Check if the synapse represented by SynInfo object was in one of the Sections
    that was absrbed into the given SectionRef.
    """
    return ((secref.gid==syn_info.gid)
              or (hasattr(secref, 'merged_sec_gids') and (syn_info.gid in secref.merged_sec_gids))
              or (hasattr(secref, 'orig_props') and (syn_info.gid in secref.orig_props.merged_sec_gids)))

def map_syn_to_loc(
        noderef,
        allsecrefs,
        syn_info,
        distance_func,
        distance_attr,
        did_ascend_original=False
    ):
    """
    Find an equivalent location for the synapse in the subtree of the given
    Section, using the measure of electrotonic attenuation dependent on 'method'.

    @param  method : str
            Measure of electrotonic attenuation to use for finding equivalent location.

    @param  imp : Hoc.Impedance
            Probe for measuring electrotonic attenuation

    @return map_loc_dist : tuple(nrn.Segment, float)
            Tuple containing segment with associated x-value that synapse was
            mapped to, and the value of the distance measure at that location.
    """
    cur_sec = noderef.sec
    freq = syn_info.PSP_median_frequency

    Asyn = getattr(syn_info, distance_attr)
    A_0 = distance_func(source=cur_sec(0.0), freq=freq)
    A_1 = distance_func(source=cur_sec(1.0), freq=freq)

    logger.anal("Entering section with A(0.0)={} , A(1.0)={} (target A={})".format(A_0, A_1, Asyn))

    # Assume monotonically decreasing A away from root section
    if (A_0 <= Asyn <= A_1) or (A_1 <= Asyn <= A_0):

        # Calculate A at midpoint of each internal segment
        locs_As = [(seg.x, distance_func(source=seg, freq=freq)) for seg in cur_sec]
        A_diffs = [abs(Asyn-pts[1]) for pts in locs_As]

        # Map synapse with closest A at midpoint
        seg_index = A_diffs.index(min(A_diffs))
        x_map, A_map = locs_As[seg_index]
        # map_seg = cur_sec(x_map)
        # seg_x0, seg_x1 = nrnutil.seg_xmin(map_seg), nrnutil.seg_xmax(map_seg)
        # seg_A0, seg_A1 = distance_func(seg_x0), distance_func(seg_x1)
        # x_map = interp(Asyn, [seg_A0, seg_A1], [seg_x0, seg_x1])

        err_A = (A_map - Asyn) / Asyn
        if not isclose(A_map, Asyn, rel_tol=.1):
            logger.warning("Warning: error in electrotonic distance measure {} "
                            "is > 10\%: error = {}".format(distance_attr, err_A))
        
        return cur_sec(x_map), A_map

    else:
        logger.anal("No map: descending further...")
        
        # Get child branches
        childrefs = [getsecref(sec, allsecrefs) for sec in noderef.sec.children()]

        # If we are in correct tree but A smaller than endpoints, return endpoint
        if not any(childrefs):
            logger.warning("Arrived at branch terminal without encountering segment with matching value of electrotonic measure")
            return cur_sec(1.0), A_1

        # Else, recursively search child nodes
        did_contain_syn = functools.partial(syn_was_in_absorbed, syn_info)

        # Did we pass the synapse's original section?
        did_ascend_original = did_ascend_original or syn_was_in_absorbed(syn_info, noderef)
        
        # child_mapsto_synsec = [subtree_has_node(mapsto_synsec, ref, allsecrefs) for ref in childrefs]
        # did_ascend_original = not any(child_mapsto_synsec)
        # if did_ascend_original:
        #   assert mapsto_synsec(noderef) # assume that synapse was positioned on this Section

        for childref in childrefs:

            # Only descend subtree if original synapse section already passed, or in subtree
            if did_ascend_original or subtree_has_node(did_contain_syn, childref, allsecrefs):
                return map_syn_to_loc(
                        childref, allsecrefs, syn_info, 
                        distance_func, distance_attr, did_ascend_original)
        
        # TODO: plot Av_soma_to_syn before and after reduction, give option to position synapses based on path length rather than 'equivalent' location. This is more logical since equivalent location may not exist, and may change after optimization of parameters.
        assert did_ascend_original
        raise Exception("The synapse did not map onto any segment in this subtree.")



def map_synapses(
        rootref,
        allsecrefs,
        orig_syn_info,
        init_cell,
        Z_freq, 
        syn_mod_pars=None,
        method=None,
        pos_method=None,
        pos_attr=None,
        linearize_gating=False):
    """
    Map synapses to equivalent synaptic inputs on given morphologically
    reduced cell.

    @param  rootref : SectionRef
            root section of reduced cell

    @param  orig_syn_info : SynInfo
            Synapse info for each synapse on original cell

    @param  method : str

            Method for positioning synapses and scaling synaptic conductances,
            specified as a measure of electrotonic extent:

                - 'Ztransfer' positions synapses at loc with ~= transfer impedance

                - 'Ai_syn_to_soma': place synapse at location with approx. 
                  the same current attenuation A_{I,syn->soma}

                - 'Av_syn_to_soma': place synapse at location with approx. 
                  the same voltage attenuation A_{V,syn->soma}

    @param  pos_method : str
            
            Optional method for measuring distance used to position synapses
            at an equivalent location. If given, synapses will be placed in a 
            segment where distance_func(segment) ~= syn_info.distance_attr. 
            Accepted values are the same as for argument 'method', as well as
            'root_distance_micron'. If not given, the measure of electrotonic 
            extent given in 'method' will be used.

    @param  pos_attr : str
            See argument 'distance_func' for meaning.

    @effect Create one synapse for each original synapse in 
            orig_syn_info. A reference to this synapse is saved 
            as as attribute 'mapped_syn' on the SynInfo object.
    """
    # Synaptic mechanisms
    if syn_mod_pars is None:
        syn_mod_pars = synmech_parnames

    # Ready cell and make Impedance probe
    logger.debug("Initializing cell for electrotonic measurements...")
    init_cell()
    logger.debug("Placing impedance measuring electrode...")
    imp = h.Impedance() # imp = h.zz # previously created

    # Get measurement frequency
    for syn in orig_syn_info:
        if not hasattr(syn, 'PSP_median_frequency'):
            syn.PSP_median_frequency = Z_freq

    # Make electrotonic distance measurement function
    if pos_method is None:
        pos_method = method
    if pos_attr is None:
        pos_attr = pos_method
    distance_func = functools.partial(
                            electrotonic_measurement_funcs[pos_method],
                            target=rootref.sec(0.5),
                            imp=imp,
                            linearize_gating=int(linearize_gating),
                            precompute=False)

    # Loop over all synapses
    logger.debug("Mapping synapses to reduced cell...")
    for syn_info in orig_syn_info:

        # Find the segment with closest match in electrotonic measure
        map_seg, map_pos = map_syn_to_loc(
                            rootref, allsecrefs, syn_info, 
                            distance_func, pos_attr)

        logger.anal("Synapse was in {}({}) -> mapped to {}\n".format(
                        syn_info.sec_name, syn_info.sec_loc, map_seg))

        # Make the synapse
        syn_mod = syn_info.mod_name
        synmech_ctor = getattr(h, syn_mod) # constructor function for synaptic mechanism
        mapped_syn = synmech_ctor(map_seg)
        syn_info.mapped_syn = mapped_syn

        # Copy synapse properties
        mech_params = syn_mod_pars[syn_mod]
        for par in mech_params:
            val = getattr(syn_info, par)
            setattr(mapped_syn, par, val)

        # Change target of afferent connections
        for aff_nc in syn_info.afferent_netcons:
            aff_nc.setpost(mapped_syn)
            aff_nc.active(1) # NetCons are turned off if target was set to None!

        # Ensure synaptic conductances scaled correctly to produce same effect at soma
        # using relationship `Vsoma = Zc*gsyn*(V-Esyn) = k_{syn->soma}*Zin*gsyn*(V-Esyn)`

        # Calculate scale factor
        #   `Vsoma = k_{syn->soma} * Vsyn`          (Vsyn is local EPSP)
        #   `Vsoma = k_{syn->soma} * Zin * Isyn`    (Isyn is synaptic current source)
        #   `Vsoma = (Zc / Zin) * Zin * Isyn`
        #   `Vsoma = Zc * Isyn`
        #   `Vsoma = Zc * gsyn * (V-Esyn)`
        if method == 'Ztransfer':
            # method 1: (don't conserve local Vsyn):
            #           - place at x with +/- same Zc 
            #           - ensure Isyn is the same
            #               - correct gsyn using exact Zc measurement (see below)

            # Vsoma = Zc(new) * gsyn * (V-Esyn)
            # ... but we want:
            # Vsoma = Zc(old) * gsyn * (V-Esyn)
            map_Zc = measure_transfer_impedance(
                            source=map_seg,
                            target=rootref.sec(0.5),
                            imp=imp,
                            freq=syn_info.PSP_median_frequency,
                            linearize_gating=int(linearize_gating),
                            precompute=False)
            scale_g = syn_info.Zc/map_Zc
        
        elif method == 'Av_syn_to_soma':
            # Vsyn = Zin(new) * gsyn * (V-Esyn)
            # ... but we want:
            # Vsyn = Zin(old) * gsyn * (V-Esyn)
            map_Zin = measure_input_impedance(
                            source=map_seg,
                            target=None,
                            imp=imp,
                            freq=syn_info.PSP_median_frequency,
                            linearize_gating=int(linearize_gating),
                            precompute=False)
            scale_g = syn_info.Zin/map_Zin

        elif method == 'Ai_syn_to_soma':
            # Isoma = Ai(map_x) * gsyn * (V-Esyn)
            # ... but we want:
            # Isoma = Ai(orig_x) * gsyn * (V-Esyn)
            map_Ai = measure_current_transfer(
                            source=map_seg,
                            target=rootref.sec(0.5),
                            imp=imp,
                            freq=syn_info.PSP_median_frequency,
                            linearize_gating=int(linearize_gating),
                            precompute=False)
            scale_g = syn_info.Ai_syn_to_soma / map_Ai
        
        else:
            raise Exception("Unknown synapse placement method '{}'".format(method))

        # Scale conductances (weights) of all incoming connections
        scale_netcon_weights = False

        # Scale parameter of synaptic mechanism
        if len(getattr(syn_info, 'gbar_param_specs', [])) > 0:
            for gbar_spec in syn_info.gbar_param_specs:

                mechtype, mech_parname, mech_paridx = configutil.interpretParamSpec(gbar_spec)
                
                if mechtype == 'syn':
                    target = mapped_syn
                    if mech_paridx is None:
                        old_val = getattr(target, mech_parname)
                        setattr(target, mech_parname, old_val*scale_g)
                    else:
                        old_val = getattr(target, mech_parname)[int(mech_paridx)]
                        getattr(target, mech_parname)[int(mech_paridx)] = old_val*scale_g
                
                elif mechtype == 'netcon':
                    scale_netcon_weights = True

        else:
            scale_netcon_weights = True
                    
        # Scale weights of associated NetCon
        if scale_netcon_weights:
            # Default action: scale NetCon weights
            for i_nc, nc in enumerate(syn_info.afferent_netcons):
                
                # Get original weights (weights are reset when disconnecting)
                orig_weights = syn_info.afferent_weights[i_nc]
                assert len(orig_weights) == int(nc.wcnt())
                
                # Scale original weights
                for i_w in xrange(int(nc.wcnt())):
                    nc.weight[i_w] = orig_weights[i_w] * scale_g
                    logger.anal("Scaled weight {} by factor {}".format(orig_weights[i_w], scale_g))
