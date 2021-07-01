"""
PyNN compatible version of Corbit/Mahon Striatal MSN model.

@author     Lucas Koelman

@date       14/08/2018
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
h.load_file("mahon_createcell.hoc")
os.chdir(prev_cwd)

from bgcellmodels.extensions.pynn.ephys_models import (
    MorphModelBase, MorphCellType)

class MsnCellModel(MorphModelBase):
    """
    Model class for Mahon/Corbit MSN cell.


    EXAMPLE
    -------

    >>> from mahon_pynn_model import MsnCellModel
    >>> MsnCellModel.default_GABA_mechanism = 'MyMechanism'
    >>> cell = MsnCellModel()
    >>> nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)
    >>> icell = cell.instantiate(sim=nrnsim)
    
    """

    # Combined with celltype.receptors in MorphCellType constructor
    # to make celltype.receptor_types in format 'region.receptor'
    regions = ['proximal']

    # Must define 'default_parameters' in associated cell type
    # parameter_names = []

    # Workaround: set directly as property on the class because
    # PyNN only allows numerical parameters
    default_GABA_mechanism = 'GABAsyn'
    default_GLU_mechanism = 'GLUsyn'
    allow_synapse_reuse = False


    def instantiate(self):
        """
        Instantiate cell in simulator

        @override       ephys.models.CellModel.instantiate()
                        
                        Since the wrapped model is a pure Hoc model completely
                        defined by its Hoc template, i.e. without ephys
                        morphology, parameters, or mechanisms definitions,
                        we have to override instantiate().
        """
        self.icell = h.MahonMSN()
        for seg in self.icell.soma[0]:
            seg.el_Leakm = -90 # Corbit (2016) changes -75 to -90


    def get_threshold(self):
        """
        Get spike threshold for soma membrane potential (used for NetCon)
        """
        return 0.0


    def get_synapses(self, region, receptors, num_contacts, **kwargs):
        """
        Get synapse in subcellular region for given receptors.
        Called by Connector object to get synapse for new connection.

        @override   MorphModelBase.get_synapse()
        """
        syns = [self.make_new_synapse(receptors, self.icell.soma[0](0.5), **kwargs) for i in xrange(num_contacts)]
        synmap_key = tuple(sorted(receptors))
        self._synapses['proximal'].setdefault(synmap_key, []).extend(syns)
        return syns



class MsnCellType(MorphCellType):
    """
    Encapsulates an MSN model described as a BluePyOpt Ephys model 
    for interoperability with PyNN.
    """

    # The encapsualted model available as class attribute 'model'
    model = MsnCellModel

    # NOTE: default_parameters is used to make 'schema' for checking & converting datatypes
    default_parameters = {}
    # extra_parameters = {}
    default_initial_values = {'v': -77.4}
    # recordable = ['spikes', 'v']

    # Combined with self.model.regions by MorphCellType constructor
    receptor_types = ['AMPA', 'NMDA', 'AMPA+NMDA',
                      'GABAA', 'GABAB', 'GABAA+GABAB']


    def can_record(self, variable):
        """
        Override or it uses pynn.neuron.record.recordable_pattern.match(variable)
        """
        return super(MsnCellType, self).can_record(variable)


def test_msn_population(export_locals=True):
    """
    Test creation of PyNN population of MSN cells.
    """
    from pyNN.utility import init_logging
    import pyNN.neuron as nrn

    init_logging(logfile=None, debug=True)
    nrn.setup()

    # STN cell population
    cell_type = MsnCellType()
    p1 = nrn.Population(5, cell_type)

    p1.initialize(v=-77.4)

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