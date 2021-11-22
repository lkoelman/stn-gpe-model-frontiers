import logging

import numpy as np
import neuron
import bluepyopt.ephys as ephys

# PyNN imports
from pyNN.parameters import ArrayParameter

from bgcellmodels.common import logutils
from bgcellmodels.morphology import morph_3d, morph_io
from bgcellmodels.models.STN import GilliesWillshaw as gillies
from bgcellmodels.models.STN import Miocinovic2006 as miocinovic
# from bgcellmodels.models.axon.mcintyre2002 import AxonMcintyre2002
from bgcellmodels.models.axon.foust2011 import AxonFoust2011
from bgcellmodels.extensions.pynn import cell_base


h = neuron.h
gillies.load_mechanisms()

logger = logging.getLogger('miocinovic_model')
logger.setLevel(logging.DEBUG)
logging.basicConfig(format=logutils.DEFAULT_FORMAT)



class GilliesSwcModel(cell_base.MorphModelBase):
    """
    Morphological SWC model that uses the channel density distributions
    and morphology from the original Gillies (2005) model.
    """

    # Combined with celltype.receptors in MorphCellType constructor
    # to make celltype.receptor_types in format 'region.receptor'
    regions = ['proximal', 'distal']

    # Related to PyNN properties
    _mechs_params_dict = {
        'STh':  ['gpas'],
        'Na':   ['gna'],
        'NaL':  ['gna'],
        'KDR':  ['gk'],
        'Kv31': ['gk'],
        'sKCa': ['gk'],
        'Ih':   ['gk'],
        'CaT':  ['gcaT'],
        'HVA':  ['gcaL', 'gcaN'],
        'Cacum':[],
    }
    rangevar_names = [p+'_'+m for m,params in _mechs_params_dict.items() for p in params]
    gleak_name = 'gpas_STh'

    # Scaling of NEURON RANGE variables
    rangevar_scaled_seclists = {}
    
    spike_threshold = {
        'soma': -10.0,
        'AIS': -10.0,
        'axon_terminal': -10.0,
    }

    # Regions for synapses
    regions = ['proximal', 'distal']
    synapse_spacing = 0.25 # (um)
    region_boundaries = {
        'proximal': (0.0, 120.0),   # (um)
        'distal':   (100.0, 1e12),  # (um)
    }

    def __init__(self, *args, **kwargs):
        # Define parameter names before calling superclass constructor
        self.parameter_names = StnMorphType.default_parameters.keys() + \
                               StnMorphType.extra_parameters.keys()
        for rangevar in self.rangevar_names:
            self.parameter_names.append(rangevar + '_scale')
        
        super(GilliesSwcModel, self).__init__(*args, **kwargs)


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator.
        """
        if sim is None:
            sim = cell_base.ephys_sim_from_pynn()

        # Get the Hoc template
        miocinovic.load_template(self.template_name) # xopen -> loads once
        template_constructor = getattr(h, self.template_name)
        
        # Instantiate template
        self.icell = icell = template_constructor()
        icell.with_extracellular = (
            self.with_extracellular_rec or self.with_extracellular_stim)

        # Load morphology into template
        morphology = ephys.morphologies.NrnFileMorphology(
                        str(self.morphology_path), do_replace_axon=False)
        morphology.instantiate(sim=sim, icell=icell)

        # Setup biophysical properties
        ais_diam = 1.2 # nodal diam from Foust axon
        ais_relative_length = 0.2
        icell.create_hillock(ais_diam, ais_relative_length, 0, 0.0)
        icell.del_unused_sections()
        icell.insert_biophys()

        # Instead of cel.set_biophys_spatial(), load channel density
        # distributions from file.
        self._init_gbar()

        # Transform morphology
        if len(self.transform) > 0 and not np.allclose(self.transform, np.eye(4)):
            morph_3d.transform_sections(icell.all, self.transform)

        # Create and append axon
        if len(self.streamline_coordinates_mm) > 0:
            self._init_axon(self.axon_class, with_ais_compartment=False)

        # Init extracellular stimulation & recording
        self._init_emfield()


    def _init_gbar(self):
        """
        Load channel conductances from Gillies & Wilshaw data files.
        """
        sec_arrays_lists = {'soma': 'somatic', 'dend': 'basal'}
        gbar_names = ['gk_KDR', 'gcaT_CaT', 'gcaL_HVA', 'gcaN_HVA', 
                        'gk_sKCa', 'gk_Kv31', 'gk_Ih']
        
        for gbar_name in gbar_names:
            gbar_mat = gillies.load_gbar_dist(gbar_name)
            for secarray_name, seclist_name in sec_arrays_lists.items():
                for i, sec in enumerate(getattr(self.icell, seclist_name)):

                    # Get location of compartment in original Gillies
                    # schematic morphology
                    tree_index, array_index = miocinovic.swc_to_gillies_index(
                                                i, secarray_name=secarray_name)

                    # Get samples for section
                    sample_mask = (gbar_mat[:,0] == tree_index) & (gbar_mat[:,1] == array_index)
                    gbar_samples = gbar_mat[sample_mask]
                    if gbar_samples.ndim == 1:
                        gbar_samples = gbar_samples[np.newaxis, :]
                    xvals_gvals = gbar_samples[:, [2, 3]]

                    # Choose closest sample to assign gbar (nseg discretization mismatch)
                    for seg in sec:
                        i_close = np.abs(xvals_gvals[:,0] - seg.x).argmin()
                        setattr(seg, gbar_name, xvals_gvals[i_close, 1])
                    # for seg_x, gbar_val in xvals_gvals:
                    #     setattr(sec(seg_x), gbar_name, gbar_val)

        # Increase persistent Na current to compensate for axon Zin
        for sec in list(self.icell.somatic) + [self.icell.axon[0]]:
            for seg in sec:
                seg.gna_NaL = self.somatic_gNaP_factor * seg.gna_NaL



class StnMorphType(cell_base.MorphCellType):
    """
    Cell type associated with a PyNN population.
    """
    model = GilliesSwcModel

    # NOTE: default_parameters is used to make 'schema' for checking and 
    #       converting datatypes. It supports only basic numpy-compatible types.
    default_parameters = {
        # Cell biophysics
        'somatic_gNaP_factor': 1.0,
        'membrane_noise_std': 0.1,
        # Morphology & 3D specification
        'morphology_path': np.array('placeholder/path'), # workaround for strings
        'transform': ArrayParameter([]),
        # Inputs
        'max_num_gpe_syn': 19,
        'max_num_ctx_syn': 30,
        'max_num_stn_syn': 10,
        
    }
    default_parameters.update(cell_base.MorphCellType._axon_parameters)
    default_parameters.update(cell_base.MorphCellType._emf_parameters)

    # NOTE: extra_parameters supports non-numpy types. 
    #       - They are are passed to model.__init__()
    #       - using fix below, they can be passed as argument to cell type 
    #       - You still cannot pass arrays of values (one value per cell) 
    extra_parameters = {
        'template_name': 'STN_morph_arcdist',
        'axon_class': AxonFoust2011,
        'default_GABA_mechanism': 'GABAsyn2',
        'default_GLU_mechanism': 'GLUsyn',
    }

    default_initial_values = {'v': -65.0}
    # recordable = ['spikes', 'v', 'lfp']

    # Combined with self.model.regions by MorphCellType constructor
    receptor_types = ['AMPA', 'NMDA', 'AMPA+NMDA',
                      'GABAA', 'GABAB', 'GABAA+GABAB']

    def __init__(self, **parameters):
        """
        Trick for allowing extra parameters as kwargs.
        """
        self.extra_parameters = {
            k: parameters.pop(k, v) for k,v in StnMorphType.extra_parameters.items()
        }
        cell_base.MorphCellType.__init__(self, **parameters)


    def can_record(self, variable):
        """
        Override or it uses pynn.neuron.record.recordable_pattern.match(variable)
        """
        return super(StnMorphType, self).can_record(variable)


################################################################################
# Testing
################################################################################


def test_simulate_population(export_locals=False):
    """
    Test PyNN model creation, running, and recording.

    @see    Based on test in:
            https://github.com/NeuralEnsemble/PyNN/blob/master/test/system/test_neuron.py
    """
    from pyNN.random import RandomDistribution
    from pyNN.utility import init_logging
    import pyNN.neuron as nrn

    init_logging(logfile=None, debug=True)
    nrn.setup()

    # GPe cell population
    parameters = {
        'template_name': 'STN_morph_arcdist',
        'morphology_path': '/home/luye/workspace/bgcellmodels/bgcellmodels/models/STN/Miocinovic2006/morphologies/type1RD_axonless-with-AIS.swc',
        'streamlines_path': '/home/luye/Documents/mri_data/Waxholm_rat_brain_atlas/WHS_DTI_v1_ALS/S56280_track_filter-ROI-STN.tck',
    }
    cell_type = StnMorphType(**parameters)
    p1 = nrn.Population(5, cell_type)

    # Modify population
    # p1.rset('Ra', RandomDistribution('uniform', low=100., high=300.))
    p1.initialize(v=-63.0)

    # Stimulation electrode
    current_source = nrn.StepCurrentSource(times=[50.0, 110.0, 150.0, 210.0],
                                           amplitudes=[0.4, 0.6, -0.2, 0.2])
    p1.inject(current_source)

    # Stimulation spike source
    p2 = nrn.Population(1, nrn.SpikeSourcePoisson(rate=100.0))
    connector = nrn.AllToAllConnector()
    syn = nrn.StaticSynapse(weight=0.1, delay=2.0)
    prj_alpha = nrn.Projection(p2, p1, connector, syn, 
                               receptor_type='distal.AMPA')

    # Recording
    # p1.record(['apical(1.0).v', 'soma(0.5).ina'])
    nrn.run(250.0)

    if export_locals:
        print("Adding to global namespace: {}".format(locals().keys()))
        globals().update(locals())


def test_instantiate_cell(export_locals=False):
    """
    Test cell instantiation without PyNN.
    """
    cell = StnMorphModel()
    sim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)
    cell.instantiate(sim=sim)

    if export_locals:
        print("Adding to global namespace: {}".format(locals().keys()))
        globals().update(locals())


if __name__ == '__main__':
    test_simulate_population(export_locals=True)