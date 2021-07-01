"""
PyNN compatible cell models for Gillies STN cell model.

@author     Lucas Koelman

@date       15/03/2018

"""

import os, os.path
import cPickle as pickle

import neuron
import numpy as np

# Load NEURON libraries, mechanisms
script_dir = os.path.dirname(__file__)
neuron.load_mechanisms(os.path.join(script_dir, 'mechanisms'))

# Load STN cell Hoc libraries
h = neuron.h
prev_cwd = os.getcwd()
os.chdir(script_dir)
h.load_file("gillies_create_factory.hoc") 
os.chdir(prev_cwd)

from bgcellmodels.common import nrnutil
from bgcellmodels.extensions.pynn import cell_base, ephys_models as ephys_pynn
from bgcellmodels.extensions.pynn.ephys_locations import SomaDistanceRangeLocation
# import lfpsim # loads Hoc functions

# Debug messages
# logutils.setLogLevel('quiet', [
#     'bluepyopt.ephys.parameters', 
#     'bluepyopt.ephys.mechanisms', 
#     'bluepyopt.ephys.morphologies'])


def define_locations():
    """
    Define locations / regions on the cell that will function as the target
    of synaptic connections.

    @return     list(SomaDistanceRangeLocation)
                List of location / region definitions.
    """

    proximal_dend = SomaDistanceRangeLocation(
        name='proximal_dend',
        seclist_name='basal',
        min_distance=5.0,
        max_distance=100.0)

    distal_dend = SomaDistanceRangeLocation(
        name='distal_dend',
        seclist_name='basal',
        min_distance=100.0,
        max_distance=1000.0)

    return [proximal_dend, distal_dend]


class StnCellModel(ephys_pynn.EphysModelWrapper):
    """
    Model class for Gillies STN cell.

    NOTES
    ------

    - instantiated using Population.cell_type.model(**parameters) 
      and assigned to ID._cell

        - parameters are those passed to the CellType on creation

    - !!! Don't forget to set initial ion concentrations globally.


    EXAMPLE
    -------

    >>> cell = StnCellModel()
    >>> nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)
    >>> icell = cell.instantiate(sim=nrnsim)
    
    """

    # _ephys_locations = define_locations()


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
    rangevar_names = [p+'_'+m for m,params in _mechs_params_dict.iteritems() for p in params]
    gleak_name = 'gpas_STh'

    tau_m_scaled_regions = ['somatic', 'basal', 'apical', 'axonal']
    rangevar_scaled_seclists = {}

    spike_threshold = {'soma': -10.0}

    # Combined with celltype.receptors in MorphCellType constructor
    # to make celltype.receptor_types in format 'region.receptor'
    regions = ['proximal', 'distal']
    synapse_spacing = 0.25 # (um)
    region_boundaries = {
        'proximal': (0.0, 120.0),   # (um)
        'distal':   (100.0, 1e12),  # (um)
    }

    def __init__(self, *args, **kwargs):
        # Define parameter names before calling superclass constructor
        self.parameter_names = StnCellType.default_parameters.keys()
        for rangevar in self.rangevar_names:
            self.parameter_names.append(rangevar + '_scale')
        
        super(StnCellModel, self).__init__(*args, **kwargs)


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator

        @override       ephys.models.CellModel.instantiate()
                        
                        Since the wrapped model is a pure Hoc model completely
                        defined by its Hoc template, i.e. without ephys
                        morphology, parameters, or mechanisms definitions,
                        we have to override instantiate().
        """
        cell_idx = h.make_stn_cell_global()
        cell_idx = int(cell_idx)
        self.icell = h.SThcells[cell_idx]


    def memb_init(self):
        """
        Initialization function required by PyNN.

        @override     EphysModelWrapper.memb_init()
        """
        super(StnCellModel, self).memb_init()

        for sec in self.icell.all:
            h.ion_style("na_ion",1,2,1,0,1, sec=sec)
            h.ion_style("k_ion",1,2,1,0,1, sec=sec)
            h.ion_style("ca_ion",3,2,1,1,1, sec=sec)


    def _update_position(self, xyz):
        """
        Called when the cell's position is changed, e.g. when changing 
        the space/structure of the parent Population.

        @effect     Adds xyz to all coordinates of the root sections and then
                    calls h.define_shape() so that whole tree is translated.
        """
        if self.calculate_lfp:
            # translate the root section and re-define shape to translate entire cell
            # source_ref = h.SectionRef(sec=self.source_section)
            # root_sec = source_ref.root
            root_sec = self.icell.soma[0]
            # root_sec.push() # 3D functions operate on CAS

            def get_coordinate(i, sec):
                return [h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)]

            # initial define shape to make sure 3D info is present
            h.define_shape(sec=root_sec)
            # root_origin = get_coordinate(0, root_sec)

            # FIXME: uncomment cell position update after test
            # # Translate each point of root_sec so that first point is at xyz
            # for i in range(int(h.n3d(sec=root_sec))):
            #     old_xyz = get_coordinate(i, root_sec)
            #     old_diam = h.diam3d(i, sec=root_sec)
            #     h.pt3dchange(i,
            #         xyz[0]-root_origin[0]+old_xyz[0],
            #         xyz[1]-root_origin[1]+old_xyz[0],
            #         xyz[2]-root_origin[2]+old_xyz[0],
            #         old_diam)

            # # redefine shape to translate tree based on updated root position
            # h.define_shape(sec=root_sec)
            # # h.pop_section()


    def _init_lfp(self):
        """
        Initialize LFP sources for this cell.

        @pre        Named parameter 'lfp_electrode_coords' must be passed
                    to cell

        @return     lfp_tracker : nrn.POINT_PROCESS
                    Object with recordable variable 'summed' that represents 
                    the cell's summed LFP contributions
        """
        if not self.calculate_lfp:
            self.lfp = None
            return
        
        # define_shape() called in _update_position()
        # h.define_shape(sec=self.icell.soma[0])
        coords = h.Vector([self.lfp_electrode_x, 
                           self.lfp_electrode_y,
                           self.lfp_electrode_z])
        sigma = self.lfp_sigma_extracellular

        # Method A: using LFP summator object
        self.lfp_tracker = h.LfpTracker(
                                self.icell.soma[0], True, "PSA", sigma, 
                                coords, self.icell.somatic, self.icell.basal)


class StnCellReduced(StnCellModel):
    """
    Model class for Gillies STN cell.

    DEVNOTES
    --------

    - instantiated using following call stack:
        - ID._cell = Population.cell_type.model(**parameters)
            - StnCellReduced.__init__()
                - EphysModelWrapper.__init__()
                    - MorphModelBase.__init__()
                        - StnCellReduced.instantiate()

        - parameters are those passed to the CellType on creation
    
    """

    def __init__(self, *args, **kwargs):
        # Define parameter names before calling superclass constructor
        self.parameter_names = StnReducedType.default_parameters.keys()
        for rangevar in self.rangevar_names:
            self.parameter_names.append(rangevar + '_scale')

        self._pickle_file = kwargs.pop('cell_pickle_file')
        
        super(StnCellModel, self).__init__(*args, **kwargs) # skip superclass


    class NrnStnProto(object):
        """ Container for STN cell, like NEURON template/prototype """
        pass


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator

        NOTE: called automatically by base class __init__

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
        icell = self.NrnStnProto() # stdutil.Bunch()

        # Copy to template
        # NOTE: SecList does not keep refs alive -> need to store lists
        for sl_name, array_name in seclists_arrays.items():
            if array_name is not None:
                setattr(icell, array_name, list(saved_seclists.get(sl_name, [])))
            seclist = h.SectionList()
            setattr(icell, sl_name, seclist)
            for sec in saved_seclists.get(sl_name, []):
                seclist.append(sec=sec)

        self.icell = icell
        return icell


class StnCellType(cell_base.MorphCellType):
    """
    Encapsulates an STN model described as a BluePyOpt Ephys model 
    for interoperability with PyNN.
    """

    # The encapsualted model available as class attribute 'model'
    model = StnCellModel

    # NOTE: default_parameters is used to make 'schema' for checking & converting datatypes
    default_parameters = {
        # 'default_GABA_mechanism': 'GABAsyn',
        # 'default_GLU_mechanism': 'GLUsyn',
        'calculate_lfp': False,
        'lfp_sigma_extracellular': 0.3,
        'lfp_electrode_x': 100.0,
        'lfp_electrode_y': 100.0,
        'lfp_electrode_z': 100.0,
        'tau_m_scale': 1.0,
        'membrane_noise_std': 0.0,
        'max_num_gpe_syn': 19,
        'max_num_ctx_syn': 30,
        'max_num_stn_syn': 10,
        'default_GABA_mechanism': np.array('GABAsyn2'),
        'default_GLU_mechanism': np.array('GLUsyn'),
    }
    # extra_parameters = {}
    default_initial_values = {'v': -65.0}
    # recordable = ['spikes', 'v', 'lfp']

    # Combined with self.model.regions by MorphCellType constructor
    receptor_types = ['AMPA', 'NMDA', 'AMPA+NMDA',
                      'GABAA', 'GABAB', 'GABAA+GABAB']

    def get_schema(self):
        """
        Get mapping of parameter names to allowed parameter types.
        """
        # Avoids specifying default values for each scale parameter and
        # thereby calling the property setter for each of them
        schema = super(StnCellType, self).get_schema()
        schema.update({v+'_scale': float for v in self.model.rangevar_names})
        return schema


    def can_record(self, variable):
        """
        Override or it uses pynn.neuron.record.recordable_pattern.match(variable)
        """
        return super(StnCellType, self).can_record(variable)


class StnReducedType(StnCellType):
    """
    PyNN cell type for reduced STN model.
    """
    # All parameters the same as superclass
    model = StnCellReduced
    
    default_parameters = dict(StnCellType.default_parameters)
    default_parameters['cell_pickle_file'] = np.array(
        '~/workspace/bgcellmodels/bgcellmodels/models/STN/GilliesWillshaw/reduced/stn-cell_Gillies2005_reduced-BushSejnowski.pkl')


def test_stn_cells_multiple(export_locals=True):
    """
    Test creation of multiple instances of the Gillies STN model.
    """
    stn_cells = []
    for i in range(5):
        cell_idx = h.make_stn_cell_global()
        cell_idx = int(cell_idx)
        stn_cells.append(h.SThcells[cell_idx])

    if export_locals:
        print("Adding to global namespace: {}".format(locals().keys()))
        globals().update(locals())


def test_stn_pynn_population(export_locals=True):
    """
    Test creation of PyNN population of STN cells.
    """
    from pyNN.utility import init_logging
    import pyNN.neuron as nrn

    init_logging(logfile=None, debug=True)
    nrn.setup()

    # STN cell population
    cell_type = StnCellType()
    p1 = nrn.Population(5, cell_type)

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
    # test_stn_cells_multiple()
    # test_stn_pynn_population()

    # Make single cell
    import bluepyopt.ephys as ephys

    # cell = StnCellModel()
    cell = StnCellReduced(
        cell_pickle_file='/home/luye/cloudstore_m/simdata/Gillies2005_reduced/bush_sejnowski_tapered/stn-cell_Gillies2005_reduced-BushSejnowski.pkl',
        owning_gid=1)

    icell = cell.icell
