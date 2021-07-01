"""
Creation of EFEL parameters for STN model optimization.

@author Lucas Koelman

@date   3/10/2017

"""

import bluepyopt.ephys as ephys

from bgcellmodels.extensions.bluepyopt.bpop_parameters import (
    NrnScaleRangeParameter, DistanceParamaterizedRangeVar)

# Gillies & Willshaw model mechanisms
from bgcellmodels.models.STN.GilliesWillshaw import gillies_model
gleak_name = gillies_model.gleak_name

################################################################################
# MODEL REGIONS
################################################################################

# seclist_name are names of SectionList declared in the cell model we optimize

somatic_region = ephys.locations.NrnSeclistLocation('somatic', seclist_name='somatic')

dendritic_region = ephys.locations.NrnSeclistLocation('dendritic', seclist_name='dendritic')

################################################################################
# PASSIVE PARAMETERS
################################################################################

# NOTE: use parameter.value attribute for default value

# SOMATIC ----------------------------------------------------------------------

soma_gl_param = ephys.parameters.NrnSectionParameter(                                    
                        name='gleak_soma',      # assigned name
                        param_name=gleak_name,  # NEURON name
                        locations=[somatic_region],
                        value = 7.84112e-05, # default value
                        bounds=[7.84112e-7, 7.84112e-3],
                        frozen=False)

soma_cm_param = ephys.parameters.NrnSectionParameter(
                        name='cm_soma',
                        param_name='cm',
                        value = 1.0,
                        bounds=[0.05, 10.0],
                        locations=[somatic_region],
                        frozen=False)


# DENDRITIC --------------------------------------------------------------------

# Cell reductions: gbar and cm have been scaled in distance-dependent manner.
# To keep distance-dependent profile, use parameter scalers instead of direct
# parameter assignment.


dend_gl_factor = NrnScaleRangeParameter(
                    name='gleak_dend_scale',
                    param_name=gleak_name,
                    value = 1.0,
                    bounds=[0.05, 10.0],
                    locations=[dendritic_region],
                    frozen=False)

dend_cm_factor = NrnScaleRangeParameter(
                    name='cm_dend_scale',
                    param_name='cm',
                    value = 1.0,
                    bounds=[0.05, 10.0],
                    locations=[dendritic_region],
                    frozen=False)

dend_ra_param = ephys.parameters.NrnSectionParameter(
                    name='Ra_dend',
                    param_name='Ra',
                    value = 150.224, # default in Gillies model
                    bounds=[50, 500.0],
                    locations=[dendritic_region],
                    frozen=False)

# Groups of parameters to be used in optimizations
soma_passive_params = [soma_gl_param, soma_cm_param]
dend_passive_params = [dend_gl_factor, dend_cm_factor, dend_ra_param]
all_passive_params = soma_passive_params + dend_passive_params

################################################################################
# ACTIVE PARAMETERS
################################################################################

# Scalers ----------------------------------------------------------------------

# for set of most important active conductance: scale factor
scaled_gbar = ['gna_NaL', 'gk_Ih', 'gk_sKCa', 'gcaT_CaT', 'gcaL_HVA', 'gcaN_HVA']

# Make parameters to scale channel conductances
MIN_SCALE_GBAR = 0.1
MAX_SCALE_GBAR = 10.0
dend_gbar_scale_params = []
for gbar_name in scaled_gbar:

    gbar_scale_param = NrnScaleRangeParameter(
                        name        = gbar_name + '_dend_scale',
                        param_name  = gbar_name,
                        value       = 1.0,
                        bounds      = [MIN_SCALE_GBAR, MAX_SCALE_GBAR],
                        locations   = [dendritic_region],
                        frozen      = False)

    dend_gbar_scale_params.append(gbar_scale_param)

# Distributions ----------------------------------------------------------------

# NOTE: 
# - the default way to solve distributions in BluePyOpt is to use 
#   `ephys.parameters.NrnRangeParameter`
#
# - its argument 'value_scaler' can be set set to one of the scale functions in 
#   `ephys.parameterscalers`.
#
# - Parameters of the scale function can be controlled by using
#   `ephys.parameters.MetaParameter`, as shown in the example
#   `bluepyopt/tests/test_ephys/test_parameters.py`
#
# - HOWEVER: we use a more concise and flexible solution, that is implemented
#   in `bgcellmodels.extensions.bluepyopt.bpop_parameters`

# TODO: linked parameters for gbar distribution
# Sub-parameters (MetaParameter)
# - assigning value stores it
# - instantitating does nothing
# DistributionParameter
# - assigning value does nothing
# - instantiating reads values from sub-parameters

def gillies_density_tapered(dist_soma, B=None, C=None, D=None, max_dist=None):
    """
    Parameterized ion channel density according to Gillies & Willshaw (2005)

    @see    models/STN/Miocinovic2006/stn_proto_arcdist.hoc
    """
    dist = dist_soma / max_dist # normalized distance
    prox = D

    if prox == 0:
        f = 1
    elif prox > 0:
        f = (1 - prox - dist) / (1 - prox)
    else:
        f = (prox + dist) / (1 + prox)

    if f < 0:
        f = 0.0

    return C + (f * B)


def gillies_density_step(dist_soma, g1=0, g2=0, g3=0, d1=0, d2=0):
    """
    Parameterized ion channel used in Gillies & Willshaw (2005)

    Step-function with three levels:
    - g1 in 0 <= x < d1
    - g2 in d1 <= x < d2
    - g3 in d2 <= x < inf

    Used for gk_KDR and gk_sKCa
    """
    if 0.0 <= dist_soma < d1:
        return g1
    elif d1 <= dist_soma < d2:
        return g2
    else:
        assert d2 <= dist_soma
        return g3


# A,B,C,D parameters from Gillies 2005, Fig. 2., Values in Table A1)
# A: density at the soma (not used in distribution)
# B: overall density to be distributed across the dendritic trees.
# C: proportion of the density that is uniformly distributed across the trees
# D: proximity, how the remaining density is distributed
# gbar_dist_article_params = {
#     # 'gk_KDR'  : {'B': 9.32e-5, 'C': 4.22e-5, 'D': -0.05},
#     # 'gk_sKCa'     : {'B': 3.92e-5 * 1.8, 'C': 0.0, 'D': -0.52},
#     'gk_Kv31'   : {'B': 1.0e-3, 'C': 8.91e-4, 'D': 0.8},
#     'gk_Ih'     : {'B': 5.1e-4, 'C': 0.0, 'D': -0.39},
#     'gcaT_CaT'  : {'B': 1.67e-3, 'C': 1.17e-3, 'D': -0.01},
#     'gcaL_HVA'  : {'B': 1.87e-3, 'C': 1.21e-4, 'D': -0.57},
#     'gcaN_HVA'  : {'B': 4.79e-4, 'C': 0.0, 'D': 0.5},   
# }

# To fit actual gbar distributions of 2005 model code:
# NOTE: sKCa and KDR are step-wise
gbar_dist_code_params = {
    # 'gk_KDR'  : {'B': 9.32e-5, 'C': 4.22e-5, 'D': -0.05},
    # 'gk_sKCa'     : {'B': 3.92e-5 * 5.3, 'C': 0.0, 'D': -0.52},
    'gk_Kv31'   : {'B': 1.0e-3 * 4.5, 'C': 8.91e-4, 'D': 0.8},
    'gk_Ih'     : {'B': 5.1e-4 * 4.9, 'C': 0.0, 'D': -0.39},
    'gcaT_CaT'  : {'B': 1.67e-3 * 2.3, 'C': 1.17e-3, 'D': -0.01},
    'gcaL_HVA'  : {'B': 1.87e-3 * 7.1, 'C': 1.21e-4, 'D': -0.57},
    'gcaN_HVA'  : {'B': 4.79e-4 * 3.8, 'C': 0.0, 'D': 0.5}, 
}

gbar_dist_step_params = {
    'gk_KDR'    : {'g1': 0, 'g2': 1e-4, 'g3': 2e-4, 'd1': 50.0, 'd2': 210.0},
    'gk_sKCa'   : {'g1': 0, 'g2': 1e-4, 'g3': 2e-4, 'd1': 240.0, 'd2': 350.0},
}

gbar_dist_uniform_params = {
    'gna_NaL'   : {'gbar': 8e-6 * 1.3 } # original model: 8e-6 (uniform)}
}

gNaL_base = gbar_dist_uniform_params['gna_NaL']['gbar']

def tapered_param_bounds(X, dist_params):
    """ Get bounds for B,C,D parameters """
    if X == 'B':
        return [dist_params['B']*0.1, dist_params['B']*10.0]
    elif X == 'C':
        return [0.0, dist_params['B']]
    elif X == 'D':
        return [-1.0, 1.0]
    else:
        raise ValueError(X)

gbar_NaL_param = ephys.parameters.NrnSectionParameter(
                    name        = 'gna_NaL_dend',
                    param_name  = 'gna_NaL',
                    locations   = [dendritic_region],
                    value       = gNaL_base,
                    bounds      = [gNaL_base*0.1, gNaL_base*10.0],
                    frozen      = False)

dend_gbar_dist_params = [gbar_NaL_param]

# Make distribution with A, B, C als metaparameters for each conductance
for gbar_name, dist_params in gbar_dist_code_params.items():
    
    gbar_dist_param = DistanceParamaterizedRangeVar(
                    name        = gbar_name + '_dend',
                    param_name  = gbar_name,
                    locations   = [dendritic_region],
                    dist_func   = gillies_density_tapered,
                    B           = dist_params['B'],
                    C           = dist_params['C'],
                    D           = dist_params['D'],
                    max_dist    = 410.0) # 369 in longest, 330 in second-longest dendrites

    dend_gbar_dist_params.append(gbar_dist_param)

    # Metaparameters
    for X in 'B', 'C', 'D':
        meta_param = ephys.parameters.MetaParameter(
                    name        = 'dist_{}_{}'.format(gbar_name, X),
                    obj         = gbar_dist_param,
                    attr_name   = X,
                    value       = dist_params[X],
                    bounds      = tapered_param_bounds(X, dist_params),
                    frozen      = False)

        dend_gbar_dist_params.append(meta_param)


# Make step-distributions for KDR and sKCa
for gbar_name, dist_params in gbar_dist_step_params.items():
    
    gbar_dist_param = DistanceParamaterizedRangeVar(
                    name        = gbar_name + '_dend',
                    param_name  = gbar_name,
                    locations   = [dendritic_region],
                    dist_func   = gillies_density_step,
                    **dist_params)

    dend_gbar_dist_params.append(gbar_dist_param)

    # Metaparameters
    for X,V in dist_params.items():
        if X == 'g1':
            continue
        meta_param = ephys.parameters.MetaParameter(
                    name        = 'dist_{}_{}'.format(gbar_name, X),
                    obj         = gbar_dist_param,
                    attr_name   = X,
                    value       = dist_params[X],
                    bounds      = [0.1*V, 10.0*V],
                    frozen      = False)

        dend_gbar_dist_params.append(meta_param)

# Default values for parameters
# default_params = {p.name: p.value for p in all_params}

