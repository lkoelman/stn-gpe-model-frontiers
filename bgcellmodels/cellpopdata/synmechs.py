"""
Utilities for dealing with synaptic mechanisms in NEURON.

@author     Lucas Koelman
@date       12/07/2017
"""

import types
import re

import numpy as np
import neuron
nrn = neuron.nrn # types nrn.Section and nrn.Segment
h = neuron.h

from bgcellmodels.mechanisms import synapses
from bgcellmodels.cellpopdata.physiotypes import Populations, NTReceptors, ParameterSource
Pop = Populations
Rec = NTReceptors
Src = ParameterSource


#######################################################################
# Parameter mapping for synaptic mechanisms
#######################################################################

# mapping synapse type -> parameter names
syn_par_maps = {
    'GABAsyn' : {}, # see below
    'GLUsyn' : {}, # see below
    'Exp2Syn' : {}, # see below
}

# Wrapper clases allow setting arbitrary attributes
class Exp2SynWrapper(neuron.hclass(h.Exp2Syn)):
    pass
class GABAsynWrapper(neuron.hclass(h.GABAsyn)):
    pass
class GLUsynWrapper(neuron.hclass(h.GLUsyn)):
    pass
class NetConWrapper(neuron.hclass(h.NetCon)):
    pass

syn_wrapper_classes = {
    'GABAsyn' : GABAsynWrapper,
    'GLUsyn' : GLUsynWrapper,
    'Exp2Syn' : Exp2SynWrapper,
}
# Dynamic creation:
# syn_wrapper_classes = {
#   modname: type(modname+'Wrapper',(neuron.hclass(getattr(h, modname)),), {})
#       for modname in syn_par_maps.keys()
# }


# GABAsyn.mod cam be used for GABA-A and GABA-B receptor
syn_par_maps['GABAsyn'] = {
    Rec.GABAA : {
        'Erev': 'syn:Erev_GABAA',
        'tau_rise_g': 'syn:tau_r_GABAA',
        'tau_decay_g': 'syn:tau_d_GABAA',
        'gbar': 'syn:gmax_GABAA',
        'delay': 'netcon:delay',
        'Vpre_threshold': 'netcon:threshold',
        'tau_rec_STD': 'syn:tau_rec', # NOTE: GABAA & GABAB use shared vars for depression/facilitation
        'tau_rec_STP': 'syn:tau_facil',
        'P_release_base': 'syn:U1', # initial release probability (fraction of vesicles in RRP released initially)
    },
    Rec.GABAB : {
        'Erev': 'syn:Erev_GABAB',
        'tau_rise_NT': 'syn:tau_r_GABAB', # in GABAsyn.mod, tau_r represents rise time of NT concentration that kicks off signaling cascade
        'tau_decay_NT': 'syn:tau_d_GABAB', # in GABAsyn.mod, tau_d represents decay time of NT concentration that kicks off signaling cascade
        'gbar': 'syn:gmax_GABAB',
        'delay': 'netcon:delay',
        'Vpre_threshold': 'netcon:threshold',
    },
}

# GABAsyn.mod cam be used for GABA-A and GABA-B receptor
syn_par_maps['GLUsyn'] = {
    Rec.AMPA : {
        'Erev': 'syn:e',
        'tau_rise_g': 'syn:tau_r_AMPA',
        'tau_decay_g': 'syn:tau_d_AMPA',
        'gbar': 'syn:gmax_AMPA',
        'delay': 'netcon:delay',
        'Vpre_threshold': 'netcon:threshold',
        'tau_rec_STD': 'syn:tau_rec', # NOTE: AMPA & NMDA use shared vars for depression/facilitation
        'tau_rec_STP': 'syn:tau_facil',
        'P_release_base': 'syn:U1',
    },
    Rec.NMDA : {
        'Erev': 'syn:e', # AMPA,NMDA have same reversal potential
        'tau_rise_g': 'syn:tau_r_NMDA',
        'tau_decay_g': 'syn:tau_d_NMDA',
        'gbar': 'syn:gmax_NMDA',
        'delay': 'netcon:delay',
        'Vpre_threshold': 'netcon:threshold',
    },
}

# Exp2Syn.mod can be used for any receptor
exp2syn_parmap = {
    'Erev': 'syn:e',
    'tau_rise_g': 'syn:tau1',
    'tau_decay_g': 'syn:tau2',
    'gbar': 'netcon:weight[0]',
    'delay': 'netcon:delay',
    'Vpre_threshold': 'netcon:threshold',
    
}
for rec in list(NTReceptors):
    syn_par_maps['Exp2Syn'][rec] = dict(exp2syn_parmap)


def getNrnConParamMap(mech_name):
    """
    For given synaptic mechanism (POINT_PROCESS defined in .mod file),
    get mapping from parameter name in dict getPhysioConParams()
    to parameters of the synaptic mechanism and NetCon.

    In other words: how each key in the parameter dictionary
    should be interpreted.
    """
    return dict(syn_par_maps[mech_name]) # return a copy


def get_mod_name(syn):
    """
    Get NEURON mechanism name of given synapse object

    @param  syn     HocObject: synapse POINT_PROCESS
    """
    if hasattr(syn, 'htype'):
        hoc_name = syn.htype.hname() # for wrapped HocObject, mechanism name is in htype attribute
    else:
        hoc_name = syn.hname()
    match_mod = re.search(r'^[a-zA-Z0-9]+', hoc_name)
    modname = match_mod.group()
    return modname


def getSynMechReceptors(mech_name):
    """
    Get receptor types implemented by the given synaptic mechanism.

    @param  mech_name   str: name of the NEURON mechanism

    @return             list(NTReceptors)
    """
    return syn_par_maps[mech_name].keys()


def getSynMechParamNames(mech_name):
    """
    Get parameter names for synaptic mechanism.

    I.e. the parameter names defined in the .mod file that can be changed.

    @return     list(str): list of parameter names declared in mod file
    """
    ntr_params_names = syn_par_maps[mech_name]
    mech_parnames = []

    # Get all parameter names prefixed by 'syn:'
    for ntr, params_names in ntr_params_names.items():
        for pname in params_names.values():
            matches = re.search(r'^(?P<mech>\w+):(?P<parname>\w+)(\[(?P<idx>\d+)\])?', pname)
            mechtype = matches.group('mech')
            if mechtype == 'syn':
                parname = matches.group('parname')
                mech_parnames.append(parname)

    return mech_parnames


# Aliases
get_mechanism_receptors = getSynMechReceptors
get_mech_param_names = getSynMechParamNames
get_physio2mech_map = getNrnConParamMap


def evalValueSpec(value_spec, rng=None):
    """
    Evaluate specification of a parameter value in a range of formats.

    @param value_spec   numeric value, function, or dict with parameters of distribution

    @param rng          numpy random object
    
    @return             float
    """
    if rng is None:
        rng = np.random

    if isinstance(value_spec, (float, int)):
        # parameter is numerical value
        value = value_spec

    elif isinstance(value_spec, types.FunctionType):
        # parameter is a function
        value = value_spec()

    elif isinstance(value_spec, dict):
        # parameter is described by other parameters (e.g. distribution)

        if ('min' in value_spec and 'max' in value_spec):
            lower = value_spec['min']
            upper = value_spec['max']
            value = lower + rng.rand()*(upper-lower)

        elif ('mean' in value_spec and 'deviation' in value_spec):
            lower = value_spec['mean'] - value_spec['deviation']
            upper = value_spec['mean'] + value_spec['deviation']
            value = lower + rng.rand()*(upper-lower)

        elif ('mean' in value_spec and 'stddev' in value_spec):
            value = rng.normal(value_spec['mean'], value_spec['stddev'])

        else:
            raise ValueError('Could not infer distribution from parameters in {}'.format(value_spec))
    return value


def correct_GABAsyn(syn):
    """
    Correct parameters of GABAsyn.mod so peak synaptic conductance
    is equal to value of 'gmax_' parameters
    """

    # Compensate for effect max value Hill factor and U1 on gmax_GABAA and gmax_GABAB
    if isinstance(syn, dict):
        syn['pointprocess']['gmax_GABAA'] /= syn['pointprocess']['U1']
        syn['pointprocess']['gmax_GABAB'] /= 0.21
    else:
        syn.gmax_GABAA = syn.gmax_GABAA / syn.U1
        syn.gmax_GABAB = syn.gmax_GABAB / 0.21


def correct_GLUsyn(syn):
    """
    Correct parameters of GLUsyn.mod so peak synaptic conductance
    is equal to value of 'gmax_' parameters
    """

    # Compensate for effect max value Hill factor and U1 on gmax_GABAA and gmax_GABAB
    if isinstance(syn, dict):
        syn['pointprocess']['gmax_AMPA'] /= syn['pointprocess']['U1']
        syn['pointprocess']['gmax_NMDA'] /= syn['pointprocess']['U1']
    else:
        syn.gmax_AMPA = syn.gmax_AMPA / syn.U1
        syn.gmax_NMDA = syn.gmax_NMDA / syn.U1


# Parameter correction functions for synaptic mechanisms
syn_mech_correctors = {
    'GABAsyn' : correct_GABAsyn,
    'GLUsyn' : correct_GLUsyn,
}


# @unique
class MSRCorrection(object):
    """
    Correction method to take into account multi-synaptic contacts
    (Multi Synapse Rule).
    """
    SCALE_GSYN_MSR = 0,     # Divide synaptic conductance by average number of contacts

    SCALE_NUMSYN_MSR = 1,   # Number of SYNAPSE objects = number of synapses 
                            # (observed) divided by average number of contacts
    
    SCALE_NUMSYN_GSYN = 2,  # Number of SYNAPSE objects = number needed to get
                            # total observed synaptic conductance


class SynInfo(object):
    """ Synapse info bunch/struct """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def get_synapse_data(connector, synapses, netcons):
    """
    Get list of SynInfo properties containing a reference
    to each synapse, its NetCon and pre-synaptic population.

    @param  synapses : list(Synapse)
            List of synapses where Synapse is a wrapped HocObject with
            attributes 'pre_pop' and 'post_pop'

    @param  netcons : list(Hoc.NetCon)
            List of NetCon objects containing afferent connections that connect
            to one of the synapse objects in 'synapses'

    @return             list(namedtuple)
    """
    
    syn_list = []

    # Save properties
    for syn in synapses:

        # Get connection parameters
        pre = syn.pre_pop
        post = syn.post_pop
        if not isinstance(pre, Populations):
            pre = Populations.from_descr(pre)
        if not isinstance(post, Populations):
            post = Populations.from_descr(post)

        con_par = connector.getPhysioConParams(pre, post, [Src.Default]) 
        syn_info = SynInfo()

        # HocObjects
        syn_info.orig_syn = syn
        syn_info.afferent_netcons = [nc for nc in netcons if nc.syn().same(syn)]
        
        # meta-information
        syn_info.pre_pop = pre

        # For PSP frequency: need to find the receptor types that this synaptic mechanism implements
        modname = get_mod_name(syn)
        syn_info.mod_name = modname

        syn_receptors = connector.getSynMechReceptors(modname)
        freqs = [con_par[receptor]['f_med_PSP_burst'] for receptor in syn_receptors]
        syn_info.PSP_median_frequency = max(freqs)

        # gbar parameters that need to be scaled
        syn_info.gbar_param_specs = []
        ntrs_params = getNrnConParamMap(modname)
        for ntr, syn_param_specs in ntrs_params.items():
            if 'gbar' in syn_param_specs:
                syn_info.gbar_param_specs.append(syn_param_specs['gbar'])

        # Save SynInfo object
        syn_list.append(syn_info)

    return syn_list
