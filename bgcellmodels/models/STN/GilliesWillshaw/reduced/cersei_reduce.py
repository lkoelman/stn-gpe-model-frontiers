"""
Test suite for reduce_cell.py : reductions of specific cell models.
"""

import re
import logging

# import sys
# pkg_root = ".." # root dir for our packages
# sys.path.append(pkg_root)

from bgcellmodels.common import logutils
from bgcellmodels.cersei.collapse.fold_algorithm import ReductionMethod
from bgcellmodels.cersei.collapse.fold_reduction import FoldReduction
from bgcellmodels.models.STN.GilliesWillshaw import gillies_model
from neuron import h

logger = logging.getLogger('gillies')
logutils.setLogLevel('verbose', ['gillies', 'marasco', 'folding'])

################################################################################
# Cell model-specific implementations of reduction functions
################################################################################

class StnCellReduction(FoldReduction):
    """
    Model-specific functions for folding reduction of Gillies & Willshaw (2005)
    STN cell model.
    """

    def __init__(self, **kwargs):
        """
        Make new Gillies model reduction.

        @param  balbi_motocell_id   (int) morphology file identifier
        """

        gillies_model.stn_cell_gillies()

        # Make sections accesible by name and index
        soma = h.SThcell[0].soma
        somatic = [soma]

        dendL = list(h.SThcell[0].dend0) # 0 is left tree
        dendR = list(h.SThcell[0].dend1) # 1 is right tree
        dendritic = dendL + dendR

        # Get references to root sections of the 3 identical trees
        dendR_root          = h.SThcell[0].dend1[0]
        dendL_upper_root    = h.SThcell[0].dend0[1] # root section of upper left dendrite
        dendL_lower_root    = h.SThcell[0].dend0[2] # root section of lower left dendrite
        fold_roots = [dendR_root, dendL_upper_root, dendL_lower_root]

        # Set reduction parameters
        kwargs['soma_secs'] = somatic
        kwargs['dend_secs'] = dendritic
        kwargs['fold_root_secs'] = fold_roots

        # Set all parameters
        kwargs['gleak_name'] = gillies_model.gleak_name
        kwargs['mechs_gbars_dict'] = gillies_model.gbar_dict
        kwargs['mechs_params_nogbar'] = gillies_model.mechs_params_nogbar
        
        super(StnCellReduction, self).__init__(**kwargs)


    def assign_region_label(reduction, secref):
        """
        Assign region labels to sections.
        """
        arrayname = re.sub(r"\[?\d+\]?", "", secref.sec.name())

        # Original sections
        if secref.is_original:
            if arrayname.endswith('soma'):
                secref.region_label = 'somatic'
            elif arrayname.endswith('dend'):
                secref.region_label = 'dendritic'
            else:
                raise Exception("Unrecognized original section {}".format(secref.sec))
        
        # Substituted / equivalent sections
        elif hasattr(secref, 'merged_region_labels'):
            secref.region_label = '-'.join(sorted(secref.merged_region_labels))

        elif hasattr(secref, 'orig_props'):
            secref.region_label = '-'.join(sorted(secref.orig_props.merged_region_labels))


    def fix_section_properties(self, new_sec_refs):
        """
        Fix properties of newly created sections.

        @override   abstract method FoldReduction.fix_section_properties
        """
        for ref in new_sec_refs:
            StnCellReduction.set_ion_styles(ref)


    @staticmethod
    def init_cell_steadystate():
        """
        Initialize cell for analyzing electrical properties.
        
        @note   Ideally, we should restore the following, like SaveState.restore():
                - all mechanism STATE variables
                - voltage for all segments (seg.v)
                - ion concentrations (nao,nai,ko,ki, ...)
                - reversal potentials (ena,ek, ...)
        """
        h.celsius = 35
        h.v_init = -68.0
        h.set_aCSF(4)
        h.init()


    @staticmethod
    def set_ion_styles(secref):
        """
        Set correct ion styles for each Section.

        @note   assigned to key 'set_ion_styles_func'
        """
        # Set ion styles
        secref.sec.push()
        # NOTE: h.ion_style(ion, c_style, e_style, einit, eadvance, cinit)
        h.ion_style("na_ion",1,2,1,0,1)
        h.ion_style("k_ion",1,2,1,0,1)
        h.ion_style("ca_ion",3,2,1,1,1)
        h.pop_section()


################################################################################
# Cell model-specific tweaks
################################################################################

def adjust_gbar_spontaneous(reduction):
    """
    Adjust gbar (NaL) to fix spontaneous firing rate.

    @note   put in list assigned to key 'post_tweak_funcs'

    @note   this is a manual model parameter tweak and should not be considered 
            part of the reduction
    """
    # Apply correction (TODO: remove this after fixing reduction)
    for ref in reduction.all_sec_refs:
        sec = ref.sec
        
        if sec.name().endswith('soma'):
            print("Skipping soma")
            continue
        
        logger.anal("Scaled gna_NaL in section {}".format(sec))
        for seg in sec:
            # seg.gna_NaL = 1.075 * seg.gna_NaL
            seg.gna_NaL = 8e-6 * 1.3 # full model value = uniform 8e-6


################################################################################
# Gillies Model Reduction Experiments
################################################################################


def make_reduction(method, reduction_params=None, tweak=False):
    """
    Make FoldReduction object using given collasping method

    @param  method : ReductionMethod
            Accepted values are BushSejnowski and Marasco

    @param  reduction_params : dict[str, object]
            Parameters that define the reduction (optional).
            These parameters are determined by the reduction method
            but can be overridden.
    """
    if not isinstance(method, ReductionMethod):
        method = ReductionMethod.from_str(str(method))
    reduction = StnCellReduction(method=method)

    # Common reduction parameters
    reduction.set_reduction_params({
            'Z_freq' :              25.,
            'Z_init_func' :         reduction.init_cell_steadystate,
            'Z_linearize_gating' :  False,
            'f_lambda':             100.0,
            'syn_scale_method' :    'Ai_syn_to_soma',
            'syn_position_method':  'root_distance_micron',
            })

    # Method-specific parameters
    if method == ReductionMethod.BushSejnowski:
        reduction.set_reduction_params({
            'gbar_init_method':     'area_weighted_average',
            'gbar_scale_method':    'surface_area_ratio',
            'passive_scale_method': 'surface_area_ratio',
            # 'gbar_scale_method':    'match_input_impedance_subtrees',
            # 'passive_scale_method': 'match_input_impedance_subtrees',
            ### Splitting cylinders based on L/lambda ##########################
            # 'split_criterion':      'eq_electrotonic_distance',
            # 'split_dX':             3.0,
            # 'lookahead_units':      'lambda',
            # 'lookahead_dX':         3.0,
            ### Splitting cylinders based on dL in micron ######################
            'split_criterion':      'micron_distance',
            'split_dX':             50.0,
        })
    
    elif method == ReductionMethod.Marasco:
        reduction.set_reduction_params({
            'gbar_scaling' :        'area',
            'set_ion_styles_func':  reduction.set_ion_styles,
            'post_tweak_funcs' :    [adjust_gbar_spontaneous] if tweak else [],
        })
    
    else:
        raise ValueError("Reduction method {} not supported".format(method))

    # Apply addition parameters (override)
    if reduction_params is not None:
        reduction.set_reduction_params(reduction_params)

    return reduction


def fold_bush(export_locals=False):
    """
    Fold Gillies STN model using Bush & Sejnowski method
    
    @param  export_locals : bool

            if True, local variables will be exported to the global
            namespace for easy inspection
    """
    # Make reduction object
    reduction = make_reduction(method=ReductionMethod.BushSejnowski)
    
    # Do reduction
    reduction.reduce_model(num_passes=1)

    if export_locals:
        globals().update(locals())
    return reduction._soma_refs, reduction._dend_refs


def fold_marasco(export_locals=False):
    """
    Fold Gillies STN model using Marasco method
    
    @param  export_locals : bool
    
            if True, local variables will be exported to the global
            namespace for easy inspection
    """
    # Make reduction object
    reduction = make_reduction(
                    method=ReductionMethod.Marasco,
                    tweak=True)
    
    # Do reduction
    reduction.reduce_model(num_passes=7)

    if export_locals:
        globals().update(locals())
    return reduction._soma_refs, reduction._dend_refs


if __name__ == '__main__':

    # Reduce cell
    reduction_method = ReductionMethod.BushSejnowski
    # reduction = make_reduction(reduction_method, reduction_params=None)
    # reduction.reduce_model(num_passes=1, map_synapses=False)

    cell_pkl_file = 'stn-cell_Gillies2005_reduced-{}.pkl'.format(
                        str(reduction_method)[16:])

    # Save cell
    # reduction.pickle_reduced_cell(cell_pkl_file)

    # Load cell
    # from bgcellmodels.morphology import morph_io
    # import pickle
    # with open(cell_pkl_file, 'rb') as file:
    #     cell_data = pickle.load(file)
    # seclists = morph_io.cell_from_dict(cell_data)