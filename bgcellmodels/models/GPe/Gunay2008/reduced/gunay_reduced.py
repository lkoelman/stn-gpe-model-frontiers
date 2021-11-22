"""
Morphology reduction of Gunay (2008) GPe cell model
"""
import logging

from bgcellmodels.common import logutils
logging.basicConfig(level=logging.WARNING, format=logutils.DEFAULT_FORMAT)

from bgcellmodels.cersei.collapse.fold_algorithm import ReductionMethod
from bgcellmodels.cersei.collapse.fold_reduction import FoldReduction
from bgcellmodels.models.GPe.Gunay2008 import gunay_model
from neuron import h

logger = logutils.getBasicLogger('gunay')
logutils.setLogLevel('warning', ['gunay', 'bpop_ext', 'matplotlib',
    'bluepyopt.ephys.parameters', 'bluepyopt.ephys.mechanisms'])

################################################################################
# Cell model-specific implementations of reduction functions
################################################################################

class GpeCellReduction(FoldReduction):
    """
    Model-specific functions for folding reduction of Gillies & Willshaw (2005)
    STN cell model.
    """

    def __init__(self, **kwargs):
        """
        Make new Gillies model reduction.

        @param  balbi_motocell_id   (int) morphology file identifier
        """

        ephys_cell = kwargs.pop('ephys_model', None)
        if ephys_cell is None:
            ephys_cell, nrnsim = gunay_model.create_cell(
                                    model=gunay_model.MODEL_GUNAY2008_AXONLESS)
        icell = ephys_cell.icell
        self.ephys_cell = ephys_cell # save ref.

        # Dendrites are connected

        # Set reduction parameters
        kwargs['soma_secs'] = list(icell.somatic)
        kwargs['dend_secs'] = list(icell.basal)
        kwargs['axon_secs'] = list(icell.axonal)
        kwargs['fold_root_secs'] = [
            # sec for sec in soma.children() if sec not in icell.axonal
            icell.dend[0], icell.dend[11], icell.dend[22]
        ]

        # Set all parameters
        kwargs['gleak_name']            = gunay_model.gleak_name
        kwargs['mechs_gbars_dict']      = gunay_model.gbar_dict
        kwargs['mechs_params_nogbar']   = gunay_model.mechs_params_nogbar
        
        super(GpeCellReduction, self).__init__(**kwargs)


    def assign_region_label(reduction, secref):
        """
        Assign region labels to sections.

        @override   FoldReduction.assign_region_label
        """
        # Original sections
        if secref.is_original:
            if 'soma' in secref.sec.name():
                secref.region_label = 'somatic'
            elif 'dend' in secref.sec.name():
                secref.region_label = 'dendritic'
            elif 'axon' in secref.sec.name():
                secref.region_label = 'axonal'
            else:
                raise Exception("Unrecognized original section {}".format(secref.sec))

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
            GpeCellReduction.set_ion_styles(ref)
            ref.sec.ena = 50.0
            ref.sec.ek = -90.0


    @staticmethod
    def init_cell_steadystate():
        """
        Initialize cell for analyzing electrical properties.
        """
        h.v_init = -68.0
        h.init()


    @staticmethod
    def set_ion_styles(secref):
        """
        Set correct ion styles for each Section.

        @note   assigned to key 'set_ion_styles_func'
        """
        # Should be done automatically based on inserted mechanisms. If not:
        # h.ion_style('na_ion', 0, 1, 0, 0, 0, sec=secref.sec)
        # h.ion_style('k_ion', 0, 1, 0, 0, 0, sec=secref.sec)
        # if 'axon' not in secref.sec.name():
        #     h.ion_style('ca_ion', 3, 2, 1, 1, 1, sec=secref.sec)
        pass



def make_reduction(method, reduction_params=None, ephys_model=None, tweak=False):
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
    
    reduction = GpeCellReduction(method=method, ephys_model=ephys_model)

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
            'split_criterion':      'micron_distance',
            'split_dX':             2.0,
        })
    
    elif method == ReductionMethod.Marasco:
        reduction.set_reduction_params({
            'gbar_scaling' :        'area',
            'set_ion_styles_func':  reduction.set_ion_styles,
            'post_tweak_funcs' :    [],
        })
    
    else:
        print("WARNING: ensure all parameters for reduction method {} "
              "are included in dict: {}".format(method, reduction_params))

    # Apply addition parameters (override)
    if reduction_params is not None:
        reduction.set_reduction_params(reduction_params)

    return reduction


if __name__ == '__main__':
    # Reduce cell
    reduction_method = ReductionMethod.BushSejnowski
    reduction_params = {'split_dX': 1.0}
    reduction = make_reduction(reduction_method, reduction_params)
    reduction.reduce_model(num_passes=1, map_synapses=False)

    cell_pkl_file = 'gpe-cell_gunay2008_reduce-{}_dL-{:.0f}.pkl'.format(
                        str(reduction_method)[16:],
                        reduction_params['split_dX'])

    # # Save cell
    from bgcellmodels.common import electrotonic
    icell = reduction.make_icell()
    electrotonic.set_min_nseg_hines(icell.all, 100.0, remove=True)
    reduction.pickle_reduced_cell(cell_pkl_file, icell=icell)

    # # Load cell
    # from bgcellmodels.morphology import morph_io
    # import pickle
    # with open(cell_pkl_file, 'rb') as file:
    #     cell_data = pickle.load(file)
    # seclists = morph_io.cell_from_dict(cell_data)
