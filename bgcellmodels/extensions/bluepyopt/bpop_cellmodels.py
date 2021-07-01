"""
Custom cell models for BluePyOpt.

@author Lucas Koelman
"""

import collections
import bluepyopt.ephys as ephys

from bgcellmodels.common import nrnutil
from bgcellmodels.morphology import morph_io
import cPickle as pickle


class PickledCellModel(ephys.models.Model):
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
            cell_pkl_file=None,
            gid=0):
        """
        Constructor

        Args:
            mechs (list of Mechanisms):
                Mechanisms associated with the cell
            
            params (list of Parameters):
                Parameters of the cell model
        """

        super(PickledCellModel, self).__init__(name)

        # BluePyOpt variables
        self.mechanisms = [] if (mechs is None) else mechs
        self.params = collections.OrderedDict()
        
        if params is not None:
            for param in params:
                self.params[param.name] = param

        self.param_values = None
        self.gid = gid
        self.cell_pkl_file = cell_pkl_file

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

        # Load cell
        with open(self.cell_pkl_file, 'rb') as file:
            cell_data = pickle.load(file)
        saved_seclists = morph_io.cell_from_dict(cell_data)
        self.icell = nrnutil.ICell(**saved_seclists)

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