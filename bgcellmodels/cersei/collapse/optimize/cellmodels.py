"""
Extensions to BluePyOpt classes for optimizing reduced morphology models.

@author Lucas Koelman

@date   13/09/2017


ARCHITECTURE

Different approaches to creating a cell model

- create own ephys.Model or ephys.models.CellModel

    + see CellModel class
        * https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/ephys/models.py
    
    + see external use cases
        * https://github.com/BlueBrain/BluePyOpt/wiki
        * particularly the example at https://github.com/apdavison/BluePyOpt/blob/pynn-models/bluepyopt/ephys_pyNN/models.py

    + see example of complex model at https://github.com/BlueBrain/BluePyOpt/tree/master/examples/l5pc

    + see dummy cell model at https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/tests/test_ephys/testmodels/dummycells.py

"""

import collections
import bluepyopt.ephys as ephys


class Cell(object):
    """
    Instantiated cell class, substitute/mock for Hoc template.
    """

    def __init__(self):
        """Constructor"""
        # list(Section) for each region
        self.soma = None
        self.dendrite = None

        # SectionList for each region
        self.somatic = None
        self.dendritic = None

        # list(SectionRef) for each region
        self._soma_refs = None
        self._dend_refs = None
        self._all_refs = None


class CollapsableCell(ephys.models.Model):
    '''
    Wraps the Gillies & Willshaw model so it can be used by BluePyOpt.

    This is version of ephys.models.CellModel without a morphology file,
    that instantiates the STN cell using the Gillies STN module.

    Based on:
        - ephys.models.CellModel
        - dummy cell template at https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/tests/test_ephys/testmodels/dummycells.py

    '''

    def __init__(self, 
            name=None,
            mechs=None,
            params=None,
            gid=0):
        """
        Constructor

        Args:
            mechs (list of Mechanisms):
                Mechanisms associated with the cell
            
            params (list of Parameters):
                Parameters of the cell model
        """

        super(CollapsableCell, self).__init__(name)

        # BluePyOpt variables
        self.mechanisms = [] if (mechs is None) else mechs
        self.params = collections.OrderedDict()
        
        if params is not None:
            for param in params:
                self.params[param.name] = param

        self.param_values = None
        self.gid = gid

        self.seclist_names = ['all', 'somatic', 'dendritic'] # SectionList variables defining regions
        self.secarray_names = ['soma', 'dendrite']

        self.icell = None

        # SelfContainedProtocol interaction
        self.proto_setup_funcs = []
        self.proto_setup_kwargs = {}


    def set_params(self, params):
        """
        Convenience function to be able to set parameters after making the model
        """
        self.params = collections.OrderedDict()
        if params is not None:
            for param in params:
                self.params[param.name] = param


    def set_mechs(self, mechs):
        """
        Convenience function to be able to set mechanisms after making the model
        """
        self.mechanisms = mechs


    def add_params(self, params, warn_override=True, raise_override=False):
        """
        Add parameters but warn if they exist
        """
        for param in params:
            if param.name in self.params.keys():
                if raise_override:
                    raise ValueError('Parameter {} already exist.'.format(param.name))
                if warn_override:
                    print('WARNING: Overriding existing parameter {}'.format(param.name))

            self.params[param.name] = param


    def add_mechs(self, mechs):
        """
        Append a list of mechanisms
        """
        self.mechanisms.extend(mechs)


    def params_by_names(self, param_names):
        """
        Get parameter objects by name
        """

        return [self.params[param_name] for param_name in param_names]


    def freeze(self, param_dict):
        """
        Set params
        """

        for param_name, param_value in param_dict.items():
            self.params[param_name].freeze(param_value)


    def unfreeze(self, param_names):
        """
        Unset params
        """

        for param_name in param_names:
            self.params[param_name].unfreeze()


    def check_nonfrozen_params(self, param_names):
        """
        Check if all nonfrozen params are set
        """

        for param_name, param in self.params.items():
            if not param.frozen:
                raise Exception(
                    'CellModel: Nonfrozen param %s needs to be '
                    'set before simulation' %
                    param_name)


    def exec_proto_setup_funcs(self, icell=None):
        """
        Execute protocol setup functions passed by SelfContainedProtocol.
        """
        setup_kwargs = { 'icell': icell }
        setup_kwargs.update(self.proto_setup_kwargs)

        for setup_func in self.proto_setup_funcs:
            setup_func(**setup_kwargs)


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator

        @note   this method should only be called after the subclass instantiate function

        @post   all mechanisms and parameters are instantiated

        @post   stored all all Sections referred to by icell's SectionRef objects
                in following attributes:

                    soma        list(Section)
                    dendrite    list(Section)

                    somatic     Hoc.SectionList()
                    dendritic   Hoc.SectionList()
                    all         Hoc.SectionList()
        """

        # Location of instantiate() in call graph during optimization:
        #       - see ephys.protocols.SweepProtocol._run_func
        #       - operations: model.freeze(params) > model.instantiate() > proto.instantiate() > sim.run () > model.unfreeze()
        #           - model.freeze()
        #               - calls param.freeze() on each param
        #               - this calls param.value.setter
        #               - how this value is used: see param.instantiate()
        #           - model.instantiate()
        #               - should call mechanism.instantiate(), param.instantiate()

        # Copy sections to instantiated cell object
        self.icell.soma = [ref.sec for ref in self.icell._soma_refs]
        self.icell.dendrite = [ref.sec for ref in self.icell._dend_refs]

        # Make SectionLists (for identification of regions)
        self.icell.all = sim.neuron.h.SectionList()
        
        ## Somatic region
        self.icell.somatic = sim.neuron.h.SectionList()
        for sec in self.icell.soma:
            self.icell.somatic.append(sec=sec)
            self.icell.all.append(sec=sec)

        ## Dendritic region
        self.icell.dendritic = sim.neuron.h.SectionList()
        for sec in self.icell.dendrite:
            self.icell.dendritic.append(sec=sec)
            self.icell.all.append(sec=sec)

        # Instantiate mechanisms and parameters
        for mechanism in self.mechanisms:
            mechanism.instantiate(sim=sim, icell=self.icell)
        
        for param in self.params.values():
            param.instantiate(sim=sim, icell=self.icell)

        return self.icell


    def destroy(self, sim=None):
        """
        Destroy cell from simulator
        """

        # FIXME: uncommenting below causes crash/hang
        # self.proto_setup_kwargs = None
        # self.proto_setup_funcs = None

        self.icell = None

        for mechanism in self.mechanisms:
            mechanism.destroy(sim=sim)
        
        for param in self.params.values():
            param.destroy(sim=sim)