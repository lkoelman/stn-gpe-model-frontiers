"""
PyNN compatible version of Fujita, Kitano et al. (2011) GPe cell model.

@author     Lucas Koelman

@date       14/09/2018
"""

# Python standard library
import os, os.path, logging

# Third party libraries
import numpy as np
import neuron
from pyNN.parameters import ArrayParameter

# Our own modules
from bgcellmodels.common import logutils
from bgcellmodels.morphology import morph_3d, morph_io
from bgcellmodels.models.axon.foust2011 import AxonFoust2011
from bgcellmodels.extensions.pynn import cell_base


# Load mechanisms and Hoc functions for cell model
h = neuron.h
script_dir = os.path.dirname(__file__)
neuron.load_mechanisms(os.path.join(script_dir, 'mechanisms'))
prev_cwd = os.getcwd()
os.chdir(script_dir)
h.xopen("fujita_createcell.hoc") # instantiates all functions & data structures on Hoc object
os.chdir(prev_cwd)

# Set up logging
logger = logging.getLogger('fujita2011')
logger.setLevel(logging.DEBUG)
logging.basicConfig(format=logutils.DEFAULT_FORMAT)


class GpeCellModel(cell_base.MorphModelBase):
    """
    Model class for Mahon/Corbit MSN cell.


    EXAMPLE
    -------

    >>> from fujita_pynn_model import GpeCellModel
    >>> GpeCellModel.default_GABA_mechanism = 'MyMechanism'
    >>> cell = GpeCellModel()
    >>> nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)
    >>> icell = cell.instantiate(sim=nrnsim)
    
    """

    # Related to PyNN properties
    _mechs_params_dict = {
        # Nonspecific channels
        'HCN':      ['gmax'],
        'leak':     ['gmax'],
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
        'CaH':    ['gmax'],
        # 'Calcium':  [''],
    }
    rangevar_names = [
        rvar + '_' + mech for mech, params in _mechs_params_dict.items() 
            for rvar in params
    ]
    gleak_name = 'gmax_leak'

    # map rangevar to cell region where it should be scaled, default is 'all'
    tau_m_scaled_regions = ['somatic']

    # Combined with celltype.receptors in MorphCellType constructor
    # to make celltype.receptor_types in format 'region.receptor'
    regions = ['proximal']
    synapse_spacing = 0.1 # single compartment so doesn't matter

    spike_threshold = {
        'soma': -10.0,
    }

    allow_synapse_reuse = False

    spike_threshold = {'soma': 0.0}


    def __init__(self, *args, **kwargs):
        """
        Populate self.parameter names with eligible parameters.
        """
        self.parameter_names = FujitaGpePrototypic.default_parameters.keys() + \
                               FujitaGpePrototypic.extra_parameters.keys()
        for rangevar in self.rangevar_names:
            self.parameter_names.append(rangevar + '_scale')
        
        super(GpeCellModel, self).__init__(*args, **kwargs)


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator

        @override       ephys.models.CellModel.instantiate()

                        Since the wrapped model is a pure Hoc model completely
                        defined by its Hoc template, i.e. without ephys
                        morphology, parameters, or mechanisms definitions,
                        we have to override instantiate().
        """
        self.icell = h.FujitaGPE()
        self.icell.setparams_corbit_2016()
        
        with_transform = len(self.transform) > 0 and not np.allclose(self.transform, np.eye(4))
        with_axon = len(self.streamline_coordinates_mm) > 0

        # Assign 3D shape
        if with_transform or with_axon:
            h.define_shape(sec=self.icell.soma[0])
            logger.debug("Assigned 3D shape to cell {}".format(self))

        # Transform morphology
        if with_transform:
            morph_3d.transform_sections(self.icell.all, self.transform)

        # Create and append axon
        if with_axon:
            self._init_axon(self.axon_class, axonmodel_use_ais=False)

        # Init extracellular stimulation & recording
        if self.with_extracellular:
            self._init_emfield()

        self.region_boundaries = {
            'proximal': (0.0, self.icell.soma[0].L),
        }


class FujitaGpePrototypic(cell_base.MorphCellType):
    """
    Encapsulates an MSN model described as a BluePyOpt Ephys model 
    for interoperability with PyNN.
    """

    # The encapsulated model available as class attribute 'model'
    model = GpeCellModel

    
    # Defaults for our custom PyNN parameters
    default_parameters = {
        'default_GABA_mechanism': np.array('GABAsyn'),
        'default_GLU_mechanism': np.array('GLUsyn'),
        'membrane_noise_std': 0.0,
        # Biophysical properties
        'tau_m_scale': 1.0,
        # Extracellular stim & rec
        'with_extracellular': False,
        'electrode_coordinates_um' : ArrayParameter([]),
        'rho_extracellular_ohm_cm' : 0.03, 
        'transfer_impedance_matrix': ArrayParameter([]),
        # 3D specification
        'transform': ArrayParameter([]),
        'streamline_coordinates_mm': ArrayParameter([]), # Sequence([])
    }

    # NOTE: extra_parameters supports non-numpy types. 
    extra_parameters = {
        'axon_class': AxonFoust2011,
    }

    default_initial_values = {'v': -65.0}
    # recordable = ['spikes', 'v']

    # Combined with self.model.regions by MorphCellType constructor
    receptor_types = ['AMPA', 'NMDA', 'AMPA+NMDA',
                      'GABAA', 'GABAB', 'GABAA+GABAB']


    def can_record(self, variable):
        """
        Override or it uses pynn.neuron.record.recordable_pattern.match(variable)
        """
        return super(FujitaGpePrototypic, self).can_record(variable)


class FujitaGpeArkypallidal(cell_base.MorphCellType):
    """
    Encapsulates an MSN model described as a BluePyOpt Ephys model 
    for interoperability with PyNN.
    """
    # The encapsulated model available as class attribute 'model'
    model = GpeCellModel

    # NOTE: default_parameters is used to make 'schema' for checking & converting datatypes
    default_parameters = {}
    # extra_parameters = {}
    default_initial_values = {'v': -65.0}
    # TODO: decrease NaP for arky type
    # recordable = ['spikes', 'v']

    # Combined with self.model.regions by MorphCellType constructor
    receptor_types = ['AMPA', 'NMDA', 'AMPA+NMDA',
                      'GABAA', 'GABAB', 'GABAA+GABAB']


    def can_record(self, variable):
        """
        Override or it uses pynn.neuron.record.recordable_pattern.match(variable)
        """
        return super(FujitaGpeArkypallidal, self).can_record(variable)


def test_gpe_population(export_locals=True):
    """
    Test creation of PyNN population of MSN cells.
    """
    from pyNN.utility import init_logging
    import pyNN.neuron as nrn

    init_logging(logfile=None, debug=True)
    nrn.setup()

    # STN cell population
    cell_type = FujitaGpePrototypic()
    p1 = nrn.Population(5, cell_type)

    # Stimulation electrode
    # Fujita et al. (2011): "Change the IClamp delay to 500, duration to
    # 1000 and current to either -.003, 0, .003, .006 nA to correspond to
    # -1, 0, 1, 2 uA/cm2."
    current_source = nrn.StepCurrentSource(times=[50.0, 110.0, 150.0, 210.0],
                                           amplitudes=[0.4, 0.6, -0.2, 0.2])
    p1.inject(current_source)

    # Recording
    # p1.record(['apical(1.0).v', 'soma(0.5).ina'])
    nrn.run(500.0)

    if export_locals:
        print("Adding to global namespace: {}".format(locals().keys()))
        globals().update(locals())


if __name__ == '__main__':
    test_gpe_population()