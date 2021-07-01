"""
PyNN compatible cell models for GPe cell model.

@author     Lucas Koelman

@date       7/03/2018

"""

import cPickle as pickle
import os.path

from neuron import h
import numpy as np

# PyNN imports
from pyNN.parameters import ArrayParameter

from bgcellmodels.common import nrnutil
from bgcellmodels.extensions.pynn import cell_base, ephys_models as ephys_pynn
from bgcellmodels.extensions.pynn.ephys_locations import SomaDistanceRangeLocation
from bgcellmodels.morphology import morph_3d
from bgcellmodels.models.axon.foust2011 import AxonFoust2011

import gunay_model

# logutils.setLogLevel('quiet', [
#     'bluepyopt.ephys.parameters',
#     'bluepyopt.ephys.mechanisms'])


def define_synapse_locations():
    """
    Define locations / regions on the cell that will function as the target
    of synaptic connections.

    @return     list(SomaDistanceRangeLocation)
                List of location / region definitions.
    """

    # NOTE: distances are not realistic since Gunay model adjusts
    #       compartment dimensions to 1 micron after loading morphology
    proximal_dend = SomaDistanceRangeLocation(
        name='proximal_dend',
        seclist_name='basal',
        min_distance=0.01,
        max_distance=1.5)

    distal_dend = SomaDistanceRangeLocation(
        name='distal_dend',
        seclist_name='basal',
        min_distance=1.5,
        max_distance=1000.0)

    return [proximal_dend, distal_dend]


class GPeCellModel(ephys_pynn.EphysModelWrapper):
    """
    Gunay (2008) multi-compartmental GPe cell model.
    """
    _ephys_morphology = gunay_model.define_morphology(
                        'morphology/bg0121b_axonless_GENESIS_import.swc',
                        replace_axon=False)
    
    _ephys_mechanisms = gunay_model.define_mechanisms(
                        'config/mechanisms.json')

    _param_locations = gunay_model.define_locations(
                        'config/locations.json')

    # Ephys parameters will automaticall by converted to PyNN parameters
    _ephys_parameters = gunay_model.define_parameters(
                            'config/params_gunay2008_GENESIS.json',
                            'config/map_params_gunay2008_v2.json',
                            _param_locations)

    # _ephys_locations = define_synapse_locations()
    regions = ['proximal', 'distal']

    # Related to PyNN properties
    _mechs_params_dict = {
        # Nonspecific channels
        'HCN':      ['gmax'],
        'HCN2':     ['gmax'],
        'pas':      ['g'],
        # Na channels
        'NaF':      ['gmax'],
        'NaP':      ['gmax'],
        # K-channels
        'Kv2':      ['gmax'],
        'Kv3':      ['gmax'],
        'Kv4f':     ['gmax'],
        'Kv4s':     ['gmax'],
        'KCNQ':     ['gmax'],
        'SK':       ['gmax'],
        # Calcium channels / buffering
        'CaHVA':    ['gmax'],
        # 'Calcium':  [''],
    }
    rangevar_names = [p+'_'+m for m,params in _mechs_params_dict.iteritems() for p in params]
    gleak_name = 'g_pas'

    # map rangevar to cell region where it should be scaled, default is 'all'
    tau_m_scaled_regions = ['somatic', 'basal', 'apical', 'axonal']
    rangevar_scaled_seclists = {} 


    # Properties defined for synapse position
    synapse_spacing = 0.2
    region_boundaries = {
        'proximal': (0.0, 2.0),   # (um)
        'distal':   (1.0, 1e12),  # (um)
    }

    # Spike threshold (mV)
    spike_threshold = {
        'soma': -10.0,
        'AIS': -10.0,
        'axon_terminal': -10.0,
    }

    def __init__(self, *args, **kwargs):
        # Define parameter names before calling superclass constructor
        self.parameter_names = GPeCellType.default_parameters.keys() + \
                               GPeCellType.extra_parameters.keys()
        for rangevar in self.rangevar_names:
            self.parameter_names.append(rangevar + '_scale')
        
        super(GPeCellModel, self).__init__(*args, **kwargs)
    

    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator

        @note   The call order is:

                - GpeCellModel.__init__()
                `- EphysModelWrapper.__init__()
                    `- ephys.models.CellModel.__init__()
                    `- cell_base.MorphModelBase.__init__()
                     `- GpeCellModel.instantiate()
                      `- EphysModelWrapper.instantiate()

        @override       ephys.models.CellModel.instantiate()
        """
        # Call instantiate method from Ephys model class
        ephys_pynn.EphysModelWrapper.instantiate(self, sim)

        # Transform morphology
        if len(self.transform) > 0 and not np.allclose(self.transform, np.eye(4)):
            morph_3d.transform_sections(self.icell.all, self.transform)

        # Create and append axon
        self.with_axon = (len(self.streamline_coordinates_mm) > 0)
        if self.with_axon:
            self._init_axon(self.axon_class)

        # Fix conductances if axon is present (compensate loading)
        # self._init_gbar()

        # Init extracellular stimulation & recording
        self._init_emfield() # BEFORE adjusting comp. dimensions.

        # Adjust compartment dimensions like in GENESIS code
        self._fix_compartment_dimensions()


    def _init_axon(self, axon_class):
        """
        Initialize axon and update source sections for spike connections.
        """
        super(GPeCellModel, self)._init_axon(axon_class,
                                             with_ais_compartment=False)


    def _fix_compartment_dimensions(self):
        """
        Normalizes all compartment dimensions to 1 um as in original
        Gunay (2008) and Hendrickson (2010) model code (see files
        names 'GPcomps.g'.

        Note that this is a huge headache since dimensions are not
        biphysically realistic anymore.
        """
        def normalize_dimensions(secs):
            for sec in secs:
                sec.L = 1.0
                for seg in sec:
                    seg.diam = 1.0

        soma = self.icell.soma[0]
        axonal_sections = list(self.icell.axonal)
        ais = axonal_sections[0]
        
        # Change dimensions depending on whether axon is attached
        if self.with_axon:
            # With axon: compensate for input impedance decrease
            myis = axonal_sections[1]
            soma.diam = 1.0
            ais.diam = 1.0
            myis.diam = 1.0
            # NOTE: section.L not scaled as opposed to normalize_dimensions()

            for seg in soma:
                seg.gmax_NaP = 6.0 * 0.000102 # last factor is default
            for seg in ais:
                seg.gmax_NaP = 6.0 * 0.004000 # last factor is default
        else:
            # Without axon: normalize like in article
            normalize_dimensions(self.icell.somatic)
            normalize_dimensions([ais])

        # Only fix section defined in SWC morphology
        normalize_dimensions(self.icell.basal)
        

    def _update_position(self, xyz):
        pass

    
class GpeCellReduced(GPeCellModel):
    """
    Morphology reduction of Gunay (2008) GPe cell model, for use with PyNN.
    """

    # NOTE: we only subclass to inherit data members, but we override
    #       al methods. Could also subclass EphysModelWrapper and
    #       copy relevant data members.
    
    # Override/hide ephys parameters
    _ephys_morphology   = None
    _ephys_mechanisms   = []
    _param_locations    = []
    _ephys_parameters   = []
        
    def __init__(self, *args, **kwargs):
        # Define parameter names before calling superclass constructor
        self.parameter_names = GpeReducedType.default_parameters.keys() + \
                               GpeReducedType.extra_parameters.keys()
        for rangevar in self.rangevar_names:
            self.parameter_names.append(rangevar + '_scale')
            
        self._pickle_file = kwargs.pop('cell_pickle_file')
        
        # skip superclass __init__
        super(GPeCellModel, self).__init__(*args, **kwargs)


    class NrnGpeProto(object):
        """ Container for GPe cell, like NEURON template/prototype """
        pass


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator

        @override       ephys.models.CellModel.instantiate()
        """

        # Load cell from file
        from bgcellmodels.morphology import morph_io
        with open(os.path.expanduser(self._pickle_file), 'rb') as file:
            cell_data = pickle.load(file)

        name_subs_regex = [(r"[\[\]]", ""), [r"\.", "_"]]

        saved_seclists = morph_io.cell_from_dict(
                            cell_data,
                            name_substitutions=name_subs_regex)

        # default NEURON secarrays and seclists for 3D morphology
        seclists_arrays = nrnutil.nrn_proto_seclists_arrays

        # Instantiate empty NEURON template
        self.icell = icell = self.NrnGpeProto() # stdutil.Bunch()

        # Copy to template
        # NOTE: SecList does not keep refs alive -> need to store lists
        for sl_name, array_name in seclists_arrays.items():
            if array_name is not None:
                setattr(icell, array_name, list(saved_seclists.get(sl_name, [])))
            seclist = h.SectionList()
            setattr(icell, sl_name, seclist)
            for sec in saved_seclists.get(sl_name, []):
                seclist.append(sec=sec)

        # Transform morphology
        if len(self.transform) > 0 and not np.allclose(self.transform, np.eye(4)):
            morph_3d.transform_sections(self.icell.all, self.transform)

        # Create and append axon
        self.with_axon = (len(self.streamline_coordinates_mm) > 0)
        if self.with_axon:
            self._init_axon(self.axon_class)

        # Init extracellular stimulation & recording
        self._init_emfield() # BEFORE adjusting comp. dimensions.

        return icell


class GPeCellType(cell_base.MorphCellType):
    """
    Encapsulates a GPe model described as a BluePyOpt Ephys model 
    for interoperability with PyNN.

    @see    Based on definition of SimpleNeuronType and standardized cell types in:
                - https://github.com/NeuralEnsemble/PyNN/blob/master/test/system/test_neuron.py
                - https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/neuron/standardmodels/cells.py

            And on documentation at:
                - http://neuralensemble.org/docs/PyNN/backends/NEURON.html#using-native-cell-models

    IMPLEMENTATION NOTES
    --------------------

    @note   The class must be callable with parameters passed as keyword arguments.
            These are passed up the inheritance chain to BaseModelType, which
            fills the parameter space.

    @note   The encapsulated model (class attribute 'model') can be any object
            as long as the __call__ method instantiates it in NEURON, and accepts
            keyword arguments containing parameters and values
    """

    # The encapsualted model available as class attribute 'model'
    model = GPeCellModel

    # Defaults for our custom PyNN parameters
    default_parameters = {
        # Morphology & 3D specification
        'transform': ArrayParameter([]),
        # Inputs
        'default_GABA_mechanism': np.array('GABAsyn'),
        'default_GLU_mechanism': np.array('GLUsyn'),
        'membrane_noise_std': 0.0,
        # Biophysical properties
        'tau_m_scale': 1.0,
    }
    default_parameters.update(cell_base.MorphCellType._axon_parameters)
    default_parameters.update(cell_base.MorphCellType._emf_parameters)

    # Defaults for Ephys parameters
    default_parameters.update({
        # ephys_model.params.values are ephys.Parameter objects
        # ephys_param.name is same as key in ephys_model.params
        p.name.replace(".", "_"): p.value for p in model._ephys_parameters
    })

    # NOTE: extra_parameters supports non-numpy types. 
    extra_parameters = {
        'axon_class': AxonFoust2011,
    }


    # extra_parameters = {}
    default_initial_values = {'v': -68.0}
    
    # Combined with self.model.regions by MorphCellType constructor
    receptor_types = ['AMPA', 'NMDA', 'AMPA+NMDA',
                      'GABAA', 'GABAB', 'GABAA+GABAB']


    def can_record(self, variable):
        """
        Override or it uses pynn.neuron.record.recordable_pattern.match(variable)
        """
        return super(GPeCellType, self).can_record(variable)


GpeProtoCellType = GPeCellType


class GpeReducedType(GPeCellType):
    """
    PyNN cell type for reduced GPe model.
    """
    # All parameters the same as superclass
    model = GpeCellReduced

    # Defaults for our custom PyNN parameters
    # NOTE: same as GpeCellType but leave out ephys parameter values
    default_parameters = {
        # Morphology & 3D specification
        'transform': ArrayParameter([]),
        # Inputs
        'default_GABA_mechanism': np.array('GABAsyn'),
        'default_GLU_mechanism': np.array('GLUsyn'),
        'membrane_noise_std': 0.0,
        # Biophysical properties
        'tau_m_scale': 1.0,
    }
    default_parameters.update(cell_base.MorphCellType._axon_parameters)
    default_parameters.update(cell_base.MorphCellType._emf_parameters)
    
    default_parameters['cell_pickle_file'] = np.array(
        '~/workspace/bgcellmodels/bgcellmodels/models/GPe/Gunay2008/reduced/gpe-cell_gunay2008_reduce-BushSejnowski_dL-1.pkl')


class GpeArkyCellType(GPeCellType):
    """
    GPe ArkyPallidal cell. It uses the same Gunay (2008) GPe cell model with
    modified parameters to reduce sponaneous firing rate and rebound firing.

    Sources
    -------

    Abdi, Mallet et al (2015) - Prototypical and Arkypallidal Neurons ...
        - See abstract: weaker persistent Na current and rebound firing
        - Fig. 7

    Bogacz, Moraud, et al (2016) - Properties of Neurons in External ...
        - Fig. 3
    """

    # Defaults for our custom PyNN parameters
    default_parameters = dict(GPeCellType.default_parameters)

    default_parameters.update({
        'gmax_NaP_scale': 0.45,
    })
    

def test_record_gpe_model(export_locals=False):
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
        'Ra_basal': 200.0
    }
    cell_type = GPeCellType(**parameters)
    p1 = nrn.Population(5, cell_type)
    
    print(p1.get('Ra_basal'))
    p1.rset('Ra_basal', RandomDistribution('uniform', low=100., high=300.))
    print(p1.get('Ra_basal'))
    
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
        receptor_type='distal_dend.Exp2Syn')

    # Recording
    # p1.record(['apical(1.0).v', 'soma(0.5).ina'])
    nrn.run(250.0)

    if export_locals:
        print("Adding to global namespace: {}".format(locals().keys()))
        globals().update(locals())


if __name__ == '__main__':
    test_record_gpe_model(export_locals=True)
