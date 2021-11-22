"""
PyNN compatible version of Golomb et al. (2007) interneuron model.

@author     Lucas Koelman

@date       20/08/2018
"""

import neuron
h = neuron.h

# Load NEURON libraries, mechanisms
import os, os.path
script_dir = os.path.dirname(__file__)
neuron.load_mechanisms(os.path.join(script_dir, 'mechanisms'))

# Load Hoc functions for cell model
prev_cwd = os.getcwd()
os.chdir(script_dir)
h.xopen("fsi_createcell.hoc") # instantiates all functions & data structures on Hoc object
os.chdir(prev_cwd)

from bgcellmodels.extensions.pynn.ephys_models import (
    MorphModelBase, MorphCellType)

# Set error tolerances for adaptive integrator
h.golomb_set_state_tolerances()


class GolombFsiModel(MorphModelBase):
    """
    Model class for Mahon/Corbit MSN cell.


    EXAMPLE
    -------

    >>> cell = GolombFsiModel()
    >>> nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)
    >>> icell = cell.instantiate(sim=nrnsim)
    
    """

    regions = ['proximal']

    # Must define 'default_parameters' in associated cell type
    # parameter_names = []

    # Workaround: set directly as property on the class because
    # PyNN only allows numerical parameters
    default_GABA_mechanism = 'GABAsyn'
    default_GLU_mechanism = 'GLUsyn'
    allow_synapse_reuse = False

    spike_threshold = {'soma': 0.0}


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator

        @override       ephys.models.CellModel.instantiate()
        """
        self.icell = h.GolombFSI()
        self.icell.setparams_corbit2016()


    def get_synapses(self, region, receptors, num_contacts, **kwargs):
        """
        Get synapse in subcellular region for given receptors.
        Called by Connector object to get synapse for new connection.

        @override   MorphModelBase.get_synapse()
        """
        syns = [self.make_new_synapse(receptors, self.icell.soma[0](0.5), **kwargs) for i in range(num_contacts)]
        synmap_key = tuple(sorted(receptors))
        self._synapses['proximal'].setdefault(synmap_key, []).extend(syns)
        return syns


class FsiCellType(MorphCellType):
    """
    Encapsulates an MSN model described as a BluePyOpt Ephys model 
    for interoperability with PyNN.
    """

    # The encapsualted model available as class attribute 'model'
    model = GolombFsiModel

    # NOTE: default_parameters is used to make 'schema' for checking & converting datatypes
    default_parameters = {}
    # extra_parameters = {}
    default_initial_values = {'v': -70.038}
    # recordable = ['spikes', 'v']

    # Combined with self.model.regions by MorphCellType constructor
    receptor_types = ['AMPA', 'NMDA', 'AMPA+NMDA',
                      'GABAA', 'GABAB', 'GABAA+GABAB']


    def can_record(self, variable):
        """
        Override or it uses pynn.neuron.record.recordable_pattern.match(variable)
        """
        return super(FsiCellType, self).can_record(variable)


def test_msn_population(export_locals=True):
    """
    Test creation of PyNN population of MSN cells.
    """
    from pyNN.utility import init_logging
    import pyNN.neuron as nrn

    init_logging(logfile=None, debug=True)
    nrn.setup()

    # STN cell population
    cell_type = FsiCellType()
    p1 = nrn.Population(5, cell_type)

    p1.initialize(v=-70.038)

    # Stimulation electrode
    current_source = nrn.StepCurrentSource(times=[50.0, 110.0, 150.0, 210.0],
                                           amplitudes=[0.4, 0.6, -0.2, 0.2])
    p1.inject(current_source)

    # Recording
    # p1.record(['apical(1.0).v', 'soma(0.5).ina'])
    nrn.run(250.0)

    if export_locals:
        print("Adding to global namespace: {}".format(locals().keys()))
        globals().update(locals())


if __name__ == '__main__':
    test_msn_population()