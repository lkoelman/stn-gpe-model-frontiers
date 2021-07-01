"""
Morphological axon models for PyNN

@author     Lucas Koelman

@date       07/05/2019

"""

import logging

# Third party modules
import numpy as np
from neuron import h
from pyNN.neuron.cells import NativeCellType
from pyNN.parameters import ArrayParameter

# Our custom modules
from bgcellmodels.extensions.pynn import cell_base
from bgcellmodels import emfield # Hoc code
from bgcellmodels.morphology import morph_3d
from bgcellmodels.extensions.pynn import cell_base
from bgcellmodels.models.axon.foust2011 import AxonFoust2011

logger = logging.getLogger('ext.pynn.cell_base')


class AxonalRelay(object):
    """
    Axon that relays incoming spikes to its terminal section as a propagating
    action potential.
    """

    def __init__(self, **kwargs):
        """
        As opposed to the original CellModel class,
        this class instantiates the cell in its __init__ method.

        ARGUMENTS
        ---------

        @param      **kwargs : dict(str, object)
                    Parameter name, value pairs
        """
        # Save parameters as attributes
        for param_name, param_val in kwargs.items():
            if isinstance(param_val, ArrayParameter):
                param_val = param_val.value
            setattr(self, param_name, param_val)

        # Build axon
        axon_builder = self.axon_class(
            self.streamline_coordinates_mm,
            termination_method=self.termination_method,
            interp_method='arclength',
            tolerance_mm=1e-4,
            without_extracellular=not (self.with_extracellular_rec or self.with_extracellular_stim))


        # Axon collaterals are specified as Nx3 matrix of branch points
        self.with_collaterals = (hasattr(self, 'collateral_branch_points_um') and
            not all([n[0] == 0 for n in self.collateral_lvl_num_branches]) and
            not all(np.isnan(self.collateral_branch_points_um).flatten())) # NaN when None is specified

        if self.with_collaterals:
            assert self.collateral_branch_points_um.ndim == 2
            assert self.collateral_branch_points_um.shape[1] == 3
            for i, branch_point_um in enumerate(self.collateral_branch_points_um):
                target_point_um = self.collateral_target_points_um[i]
                level_lengths_um = self.collateral_lvl_lengths_um[i]
                levels_num_branches = self.collateral_lvl_num_branches[i]
                step_length_um = 10.0
                step_angles_deg = np.zeros_like(level_lengths_um) + 30.0
                axon_builder.add_collateral_definition(
                                branch_point_um, target_point_um, level_lengths_um,
                                step_length_um, levels_num_branches, step_angles_deg)
        
        self.icell = axon_builder.build_axon()
                        
        # Initial section for receiving NetCon events (spikes)
        initial_sec = self.icell.main_branch[0]

        # Set terminal section for spike detection
        terminal_sec = axon_builder.get_terminal_section(self.netcon_source_spec)


        # Change source for NetCons (see pyNN.neuron.simulator code)
        self.source_section = terminal_sec
        self.source = terminal_sec(0.5)._ref_v


        # Attributes required for recording (see pyNN.neuron.recording.Recorder._record())
        self.rec = h.NetCon(self.source, None,
                            self.get_threshold(), 0.0, 0.0,
                            sec=self.source_section)
        self.spike_times = h.Vector(0)
        self.traces = {}
        self.recording_time = False

        # Create targets for NetCons
        # (property names must correspond to entries in cell_type.receptor_types)
        self.excitatory = h.Exp2Syn(initial_sec(0.5))
        # Synapse properties together with weight of 1 cause good following up to 200 Hz
        self.excitatory.tau1 = 0.1
        self.excitatory.tau2 = 0.2
        self.excitatory.e = 0.0

        # Initialize extracellular stim & rec
        self._init_emfield()


    def memb_init(self):
        """
        Set initial values for all variables in this cell.

        @override   memb_init() required by PyNN interface for cell models.
        """
        for sec in self.icell.all:
            for seg in sec:
                seg.v = self.v_init # set using pop.init(v=v_init) or default_initial_values


    def get_threshold(self):
        """
        Get spike threshold for self.source variable (usually points to membrane
        potential). This threshold is used when creating NetCon connections.

        @override   get_threshold() required by pyNN interface for cell models.

        @return     threshold : float
        """
        return -10.0


    # Bind interface method
    resolve_section = cell_base.irec_resolve_section
    _init_extracellular_stim_rec = cell_base.init_extracellular_stim_rec


    def _init_emfield(self):
        """
        Set up extracelullar stimulation and recording.

        @pre    mechanism 'extracullular' must be inserted and parameters set
                in all compartments that should contribute to the LFP and are
                targets for stimulation
        """
        # Check if we need extracellular layers for recording/stimulation
        if not (self.with_extracellular_rec or self.with_extracellular_stim):
            return
        elif self.with_extracellular_stim and self.with_extracellular_rec:
            stim_reclist = self.icell.all
            rec_seclist = self.icell.all
        elif self.with_extracellular_stim:
            stim_reclist = self.icell.all
            rec_seclist = h.SectionList()
        elif self.with_extracellular_rec:
            stim_reclist = h.SectionList()
            rec_seclist = self.icell.all

        # NOTE: 'extracellular' is alread inserted by AxonBuilder
        self._init_extracellular_stim_rec(self.icell.all, stim_reclist, rec_seclist,
                                          self.icell.main_branch[0](0.5))


    def get_all_sections(self):
        """
        Get all neuron.Section objects that make up this cell.

        @return     neuron.SectionList containing all sections
        """
        return self.icell.all


class AxonRelayType(NativeCellType):
    """
    PyNN cell type for use with AxonalRelay cell model.
    """
    model = AxonalRelay

    # Queried by Population.find_units() for recordings
    units = cell_base.UnitFetcherPlaceHolder()

    # Properties of cell model containing NetCon targets
    receptor_types = ['excitatory']

    default_parameters = {
        'transform': ArrayParameter([]),
    }
    default_parameters.update(cell_base.MorphCellType._axon_parameters)
    default_parameters.update(cell_base.MorphCellType._emf_parameters)

    # NOTE: extra_parameters supports non-numpy types. 
    extra_parameters = {
        'axon_class': AxonFoust2011,
    }

    default_initial_values = {'v': -68.0}

    def __init__(self, **parameters):
        """
        Trick for allowing extra parameters as kwargs.
        """
        self.extra_parameters = {
            k: parameters.pop(k, v) for k,v in AxonRelayType.extra_parameters.items()
        }
        NativeCellType.__init__(self, **parameters)


################################################################################
# Testing
################################################################################


def test_simulate_population(export_locals=False):
    """
    Test PyNN model creation, running, and recording.

    @see    Based on test in:
            https://github.com/NeuralEnsemble/PyNN/blob/master/test/system/test_neuron.py
    """
    from pyNN.utility import init_logging
    import pyNN.neuron as nrn
    import numpy as np

    init_logging(logfile=None, debug=True)
    nrn.setup()

    # Artificial axon trajectory
    axon_coords = np.concatenate([np.arange(20).reshape((-1,1))]*3, axis=1)
    elec_coords = np.array([1e6, 1e6, 1e6])

    # GPe cell population
    num_cell = 5
    axon_params = {
        'transform': ArrayParameter(np.eye(4)), # [np.eye(4)] * num_cell,
        'streamline_coordinates_mm': ArrayParameter(axon_coords), # [axon_coords] * num_cell,
        'axon_class': AxonFoust2011,
        'with_extracellular_stim': True,
        'with_extracellular_rec': True,
        'electrode_coordinates_um': ArrayParameter(elec_coords), # [elec_coords] * num_cell,
    }
    cell_type = AxonRelayType(**axon_params)
    pop_dest = nrn.Population(5, cell_type)
    pop_dest.initialize(v=-68.0)

    # Spike source population
    pop_src = nrn.Population(5, nrn.SpikeSourcePoisson(rate=100.0))
    
    # Connect populations
    connector = nrn.OneToOneConnector()
    syn = nrn.StaticSynapse(weight=1.0, delay=2.0)
    prj_alpha = nrn.Projection(pop_src, pop_dest,
                               connector, syn, 
                               receptor_type='excitatory')

    # Recording
    # p1.record(['apical(1.0).v', 'soma(0.5).ina'])
    nrn.run(250.0)

    if export_locals:
        globals().update(locals())


if __name__ == '__main__':
    test_simulate_population(export_locals=True)