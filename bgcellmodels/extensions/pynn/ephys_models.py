"""
Module for working with BluePyOpt cell models in PyNN.

@author     Lucas Koelman

@date       14/02/2018


USEFUL EXAMPLES
---------------

https://github.com/NeuralEnsemble/PyNN/blob/master/test/system/test_neuron.py
https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/neuron/cells.py
https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/neuron/standardmodels/cells.py
https://github.com/apdavison/BluePyOpt/blob/pynn-models/bluepyopt/ephys_pyNN/models.py

"""

# Standrad library
import re
import math
from copy import deepcopy
import logging

# Third party libraries
import bluepyopt.ephys as ephys
import pyNN.neuron
from pyNN.parameters import ArrayParameter

# Our custom modules
from bgcellmodels.extensions.pynn import cell_base

# Global variables
logger = logging.getLogger('ephys_models')
ephys_nrn_sim = None
h = pyNN.neuron.h

rng_structural_variability = h.Random(pyNN.neuron.state.mpi_rank + 
                                      pyNN.neuron.state.native_rng_baseseed)


def make_valid_attr_name(name):
    """
    Make name into a valid attribute name.

    @param      name : str

    @return     modified_name : str
                String that is a valid attribute name.
    """
    return name.replace(".", "_")



class CellModelMeta(type):
    """
    Create new type that is subclass of ephys.models.CellModel and
    automatically registers Ephys parameters as python properties.
    """

    def __new__(this_meta_class, new_class_name, new_class_bases, new_class_namespace):
        """
        Create new type for a class definition that has this metaclass
        as its metaclass.

        @effect     Converts each Ephys.Parameter object stored in subclass
                    attribute '_ephys_parameters' into a class property,
                    and stores its name in the 'parameter_names' attribute.

        @note       This method is called once every time a class is defined with 
                    MechanismType as its metaclass.
        """
        # Process mechanism variables declared in class definition
        # modified_bases = tuple(list(new_class_bases) + [ephys.models.CellModel])
        
        # Parameter names defined in class namespace (not Ephys parameters)
        parameter_names = new_class_namespace.get("parameter_names", [])
        
        for e_param in new_class_namespace.get("_ephys_parameters", []):
            param_name = e_param.name
            
            # NOTE: self.params is set in __init__()
            def __get_ephys_param(self):
                return self.params[param_name].value
            
            def __set_ephys_param(self, value):
                self.params[param_name].value = value
                self.params[param_name].instantiate(self.sim, self.icell)

            # Change name of functions for clarity
            param_name_nodots = make_valid_attr_name(param_name)
            __get_ephys_param.__name__ = "__get_" + param_name_nodots
            __set_ephys_param.__name__ = "__set_" + param_name_nodots
            
            parameter_names.append(param_name_nodots)

            # Insert into namespace as properties
            new_class_namespace[param_name_nodots] = property(
                                                    fget=__get_ephys_param,
                                                    fset=__set_ephys_param)

        # Make Ephys locations into class properties
        # location_names = []
        # for ephys_loc in new_class_namespace.get("_ephys_locations", []):

        #     loc_name = ephys_loc.name
        #     loc_name_nodots = make_valid_attr_name(loc_name)
            
        #     # NOTE: self.params is set in __init__()
        #     def __get_location(self):
        #         return self.locations[loc_name].instantiate(self.sim, self.icell)

        #     # Change name of getter function
        #     __get_location.__name__ = "__get_" + loc_name_nodots
        #     location_names.append(loc_name_nodots)

        #     # Insert into namespace as properties
        #     new_class_namespace[loc_name_nodots] = property(fget=__get_location)

        # Parameter names for pyNN
        new_class_namespace['parameter_names'] = parameter_names
        # new_class_namespace['location_names'] = location_names

        return type.__new__(this_meta_class, new_class_name, 
                            new_class_bases, new_class_namespace)


#     def __init__(new_class, new_class_name, new_class_bases, new_class_namespace):
#         """
#         Transform ephys parameters defined in subclass to properties for pyNN
# 
#         @param      new_class : type
#                     the class object returned by __new__
# 
#         @note       This method is called once every time a class is defined with 
#                     this class as its metaclass.
#         """
#         setattr(new_class, 'myprop', property(fget=__get_func, fset=__set_func))


class EphysModelWrapper(ephys.models.CellModel, cell_base.MorphModelBase):
    """
    Subclass of Ephys CellModel that conforms to the interface required
    by the 'model' attribute of a PyNN CellType class.

    This is a modified version of an Ephys CellModel that allows multiple
    instances of a cell to be created. As opposed to the original CellModel class,
    this class instantiates the cell in its __init__ method.


    @note   Called by ID._build_cell() defined in
            https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/neuron/simulator.py
            as follows:
            
                > ID._cell = Population.cell_type.model(**cell_parameters)


    @see    Based on definition of SimpleNeuronType and standardized cell types in:
                https://github.com/NeuralEnsemble/PyNN/blob/master/test/system/test_neuron.py
                https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/neuron/standardmodels/cells.py

            And on documentation at:
                http://neuralensemble.org/docs/PyNN/backends/NEURON.html#using-native-cell-models
    
    ATTRIBUTES
    ----------

    @attr   <location_name> : nrn.Section or nrn.Segment

            Each location defined on the cell has a corresponding attribute.
            This is a dynamic attribute: the Segment or Section is sampled
            from the region.


    @attr   <param_name> : float

            Each parameter defined for thsi cell has a corresponding attribute
    

    @attr   locations : dict[str, Ephys.location]

            Dict containing all locations defined on the cell.
            Keys are the same as the location attributes of this cell.


    @attr   synapses : dict[str, list(nrn.POINT_PROCESS)]

            Synapse object that synapse onto this cell.
    """

    __metaclass__ = CellModelMeta

    def __init__(self, *args, **kwargs):
        """
        As opposed to the original CellModel class,
        this class instantiates the cell in its __init__ method.

        @param      **kwargs : dict(str, object)
                    Parameter name, value pairs

        @post       self.icell contains the instantiated Hoc cell model

        """
        # Get parameter definitions from class attributes of subclass.
        model_name = self.__class__.__name__

        # Ensure self.params has valid names, and does not refer to class params
        ephys_param_defs = deepcopy(getattr(self, '_ephys_parameters', None))
        if ephys_param_defs is not None:
            for e_param in ephys_param_defs:
                e_param.name = make_valid_attr_name(e_param.name)

        # Call constructor of first base class
        # super(EphysModelWrapper, self).__init__(
        ephys.models.CellModel.__init__(self,
            model_name,
            morph=getattr(self, '_ephys_morphology', None),
            mechs=getattr(self, '_ephys_mechanisms', None),
            params=ephys_param_defs)
        ephys_param_names = self.params.keys() # ephys param dict

        # Don't pass ephys parameters to base class
        ephys_params_args = {
            k : kwargs.pop(k) for k in kwargs.keys() if k in ephys_param_names
        }

        # Call constructor of second base class
        cell_base.MorphModelBase.__init__(self, *args, **kwargs)

        # Handle post-instantiation parameters
        for param_name, param_value in ephys_params_args.iteritems():
            if isinstance(param_value, ArrayParameter):
                param_value = param_value.value
            if self.params[param_name].value != param_value:
                setattr(self, param_name, param_value)

        # # Parameters that are not dynamic must be set before instantiate
        # instantiated_params = set()
        # for param_name, param_value in kwargs.iteritems():
        #     if (param_name not in ephys_param_names) and (not param_name.endswith('_scale')):
        #         setattr(self, param_name, param_value)
        #         instantiated_params.add(param_name)

        # # Cell must be instantiated _before_ applying parameters
        # self.sim = cell_base.ephys_sim_from_pynn()
        # self.instantiate(sim=self.sim)

        # # Parse parameters passed by cell type
        # for param_name, param_value in kwargs.iteritems():
        #     if param_name in instantiated_params:
        #         continue
        #     elif (param_name in ephys_param_names) and \
        #          (self.params[param_name].value == param_value):
        #         # Ephys parameters are already set in self.instantiate()
        #         continue
        #     elif (param_name in self.parameter_names) or param_name.endswith('_scale'):
        #         # Dynamic parameters are dispatched to approprate method
        #         setattr(self, param_name, param_value)
        #     else:
        #         logger.warning("Unrecognized parameter {}. Ignoring.".format(param_name))

        # self._init_synapses()
        # self._post_instantiate()

        # # Attributes required by PyNN
        # self.source_section = self.icell.soma[0]
        # self.source = self.icell.soma[0](0.5)._ref_v
        
        # self.rec = h.NetCon(self.source, None,
        #                     self.get_threshold(), 0.0, 0.0,
        #                     sec=self.source_section)
        # self.spike_times = h.Vector(0) # see pyNN.neuron.recording.Recorder._record()
        # self.traces = {}
        # self.recording_time = False

    @staticmethod
    def create_empty_cell(
            template_name,
            sim,
            seclist_names=None,
            secarray_names=None):
        """
        Create an empty cell in Neuron

        @override   CellModel.create_empty_cell

                    The original function tries to recreate the Hoc template every
                    time so isn't suitable for instantiating multiple copies.
        """
        if not hasattr(sim.neuron.h, template_name):
            template_string = ephys.models.CellModel.create_empty_template(
                                                        template_name,
                                                        seclist_names,
                                                        secarray_names)
            sim.neuron.h(template_string)

        template_function = getattr(sim.neuron.h, template_name)

        return template_function()


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator

        The default behaviour implemented here works for cells that have ephys
        morphology, mechanism, and parameter definitions. If you want to use
        a Hoc cell without these definitions, you should subclass to override
        this behaviour.

        @override       ephys.models.CellModel.instantiate()
        """

        # Instantiate NEURON cell
        icell = self.create_empty_cell(
                    self.name,
                    sim=sim,
                    seclist_names=self.seclist_names,
                    secarray_names=self.secarray_names)

        self.icell = icell
        # icell.gid = self.gid # gid = int(ID) where ID._cell == self

        self.morphology.instantiate(sim=sim, icell=icell)

        for mechanism in self.mechanisms:
            mechanism.instantiate(sim=sim, icell=icell)

        for param in self.params.values():
            param.instantiate(sim=sim, icell=icell)


    def _post_instantiate(self):
        """
        Hook for subclass to put post-instantiation code.
        This is code to be executed after morphology and all sections
        have been created.

        @pre    cell has been instantiated in NEURON: all Sections have been
                created and mechanisms inserted.

        @pre    all parameters passed to the cell have been applied
        """
        pass


#    def get_synapse(self, region, receptors, mark_used, **kwargs):
#        """
#        Get synapse in subcellular region for given receptors.
#        Called by Connector object to get synapse for new connection.
#
#        @override   MorphModelBase.get_synapse()
#        """
#        return self.get_existing_synapse(region, receptors, mark_used)


    def _post_build(self, population, pop_index):
        """
        Hook called after Population._create_cells() -> ID._build_cell()
        is executed.

        @param  index in Population
                (same as result of Population.id_to_index(ID))

        @note   called by our custom Population class that overrides
                Population._create_cells()
        """
        cell_base.MorphModelBase._post_build(self, population, pop_index)

        # Pass position of cell in Population and use it
        position = population.positions[:, pop_index]
        self._update_position(position)
        
        self._init_lfp()


    def _init_lfp(self):
        """
        Initialize LFP sources for this cell.
        Override in subclass to implement this.

        @return     lfp_tracker : nrn.POINT_PROCESS
                    Object with recordable variable that represents the cell's
                    summed LFP contributions
        """
        pass


    def _update_position(self, xyz):
        """
        Called when the cell's position is changed, e.g. when changing 
        the space/structure of the parent Population.

        @effect     Adds xyz to all coordinates of the root sections and then
                    calls h.define_shape() so that whole tree is translated.
        """

        raise NotImplementedError("Subclass should update cell position.")


    def _set_tau_m_scale(self, value):
        """
        Setter for parameter 'tau_m'. Sets membrane time constant.

        @pre    Subclass must define attributes 'gleak_name' and
                'tau_m_scaled_regions'.
        """
        if not hasattr(self, '_tau_m_scale'):
            self._tau_m_scale = 1.0
        if value == self._tau_m_scale:
            return
        
        # If not yet scaled, save the base values for g_leak and cm
        # if self._tau_m_scale == 1.0:
        #     for region_name in scaled_regions:
        #         if region_name == 'all':
        #             continue
        #         sl = list(getattr(self.icell, region_name))
        #         self._base_gl[region_name] = sum((getattr(sec, self.gleak_name) for sec in sl)) / len(sl)
        #         self._base_cm[region_name] = sum((sec.cm for sec in sl)) / len(sl)

        # Extract the mechanism name
        matches = re.search(r'.+_([a-zA-Z0-9]+)$', self.gleak_name)
        gleak_mech = matches.groups()[0]

        # Scale tau_m in all compartments
        for region_name in self.tau_m_scaled_regions:
            for sec in getattr(self.icell, region_name):
                # If section has passive conductance: distribute scale factor 
                # over Rm (gleak) and Cm
                distribute_scale = h.ismembrane(gleak_mech, sec=sec)
                if distribute_scale:
                    cm_factor = math.sqrt(value)
                    gl_factor = 1.0 / math.sqrt(value)
                else:
                    cm_factor = value
                for seg in sec:
                    if distribute_scale:
                        setattr(seg, self.gleak_name, 
                                gl_factor * getattr(seg, self.gleak_name))
                    seg.cm = cm_factor * seg.cm
        
        self._tau_m_scale = value


    def _get_tau_m_scale(self):
        """
        Getter for parameter 'tau_m'. Get average membrane time constant.
        """
        if not hasattr(self, '_tau_m_scale'):
            self._tau_m_scale = 1.0
        return self._tau_m_scale


    tau_m_scale = property(fget=_get_tau_m_scale, fset=_set_tau_m_scale)
            

