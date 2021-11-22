"""
Create synaptic connections using custom NEURON synapse models
and get their parameters from centralized parameter databvase.

@author     Lucas Koelman

@date       21/03/2018


Example
-------

>>> pop_A = Population(num_cell, cell_type, label, initial_vals)
>>> pop_B = Population(num_cell, cell_type, label, initial_vals)
>>>
>>> connector = ConnectorType(**connector_params) # connection pattern
>>> synapse = SynapseType(**synapse_params) # has Connection type in class definition
>>>
>>> proj = Projection(pop_A, pop_B, connector=connector
>>>                   synapse_type=synapse, receptor_type='GABAA')


DEVNOTES
--------

There are three central objects involved in a PyNN connection:

- Projection    : encapsulates projection between two populations
- Connector     : the connection pattern
- SynapseType   : the synapse type and its parameters for a Projection
- Connection    : facilitates creation of simulator-specific connection
                  to instantiated SynapseType

The call graph for creating a connection in PyNN is:

-> Projection.__init__(..., Connector)
    `-> Connector.connect(..., Projection)
    `-> MapConnector._standard_connect(...)
        `-> Connector._parameters_from_synapse_type(...)
            - parameters are fetched from Projection.synapse_type
            - parameters are optionally calculated from distance map
    `-> Projection._convergent_connect(..., **connection_parameters)
    `-> Projection.synapse_type.connection_type(...)
        `-> Connection.__init__(...)
"""

import re

from pyNN.neuron.simulator import Connection, state


class NativeSynToRegion(Connection):
    """
    Connect a synaptic mechanism (MOD file name) to a cell region
    pecified as Ephys location.

    This Connection class can represent a multi-synaptic connection and hence
    encapsulate multiple NEURON synapse and NetCon objects.

    This Connection class can obly be used with NativeSynapse defined in
    this module.

    USAGE
    -----

        >>> syn_params = {'netcon:weight[0]': 1.0, 'netcon:weight[1]: 0.5',
        >>>               'netcon:delay': 3.0, 'syn:e_rev': 80.0,
        >>>               'syn:tau_fall': 15.0, 'syn:tau_rise': 5.0}
        >>>
        >>> syn = NativeSynapse(mechanism='GLUsyn',
        >>>                     mechanism_parameters=syn_params)
        >>>
        >>> proj = sim.Projection(pop_pre, pop_post, connector, syn,
                                  receptor_type='distal_region.AMPA+NMDA')

    """

    def __init__(self, projection, pre, post, **parameters):
        """
        Create a new connection.

        @pre    The region specified in the first part of the
                projection.receptor_type must be present as an attribute
                on the post-synaptic cell.
        """
        #logger.debug("Creating connection from %d to %d, weight %g" % (pre, post, parameters['weight']))
        self.presynaptic_index = pre
        self.postsynaptic_index = post
        self.presynaptic_cell = projection.pre[pre]
        self.postsynaptic_cell = projection.post[post]
        post_cell = self.postsynaptic_cell._cell # CellType.model instance

        # Get the target region on the cell and the receptor type
        # syn_mech_name = parameters.pop('mechanism', None)
        # syn_mech_name = projection.synapse_type.mechanism
        region, receptor = projection.receptor_type.split(".")
        receptors = receptor.split("+")


        # TODO: handle distribution of synapses better
        #   - cell is responsible for maintaining realistic distribution of synapses
        #       - maintain spacing, etc.
        #   - function to get multiple segments or synapses
        #       - args: pre_cell, receptor, region, number of synaptic contacts
        #   - can give existing synapses if parameters are the same
        #       - no overwrite of params in this case
        #       - check by comparing same pre-cel, receptor, region,


        # Ask cell for segment in target region
        synapse, used = post_cell.get_synapse(region, receptors, True)
        self.synapse = synapse
        if used > 0:
            raise Exception("No unused synapses on target cell {}".format())


        # Create NEURON synapse and NetCon
        pre_gid = int(self.presynaptic_cell)
        self.nc = state.parallel_context.gid_connect(pre_gid, self.synapse)

        # Interpret connection parameters
        param_targets = {'syn': self.synapse, 'netcon': self.nc}
        parameters.update(projection.synapse_type.mechanism_parameters)

        # Don't allow weight and delay in **parameters
        pynn_delay = parameters.pop('delay')
        parameters.setdefault('netcon:delay', pynn_delay)

        pynn_weight = parameters.pop('weight')
        parameters.setdefault('netcon:weight[0]', pynn_weight)

        for param_spec, value_spec in parameters.items():
            # Convert value specification to actual parameter value
            if callable(value_spec):
                # Generate distance map lazily and cache it
                distance_map = getattr(projection, '_distance_map', None)
                if distance_map is None:
                    distance_map = projection._connector._generate_distance_map(
                                                            projection)
                    projection._distance_map = distance_map
                # distance_map has projection.shape
                distance = distance_map[pre, post]
                param_value = value_spec(distance)
            elif isinstance(value_spec, (float, int)):
                param_value = value_spec
            else:
                raise ValueError("Cannot interpret parameter specification "
                    "<{}> for parameter {}".format(value_spec, param_spec))

            # interpret parameter specification in format "target:attribute[index]"
            matches = re.search(
                r'^(?P<target>\w+):(?P<parname>\w+)(\[(?P<index>\d+)\])?',
                param_spec)
            target_name = matches.group('target')
            param_name = matches.group('parname')
            param_index = matches.group('index')

            # Set attribute
            target = param_targets[target_name]
            if param_index is None:
                setattr(target, param_name, param_value)
            else:
                getattr(target, param_name)[int(param_index)] = param_value


class ConnectionNrnWrapped(Connection):
    """
    Connection to a wrapped NEURON synapse that exposes certain
    parameters to PyNN

    The advantage is that you can specify parameters using PyNN's built-in
    mechanisms. E.g.
    """

    def __init__(self, projection, pre, post, **parameters):
        """
        Create a new connection.
        """
        #logger.debug("Creating connection from %d to %d, weight %g" % (pre, post, parameters['weight']))
        self.synapse_type = projection.synapse_type
        self.presynaptic_index = pre
        self.postsynaptic_index = post
        self.presynaptic_cell = projection.pre[pre]     # ID instance
        self.postsynaptic_cell = projection.post[post]  # ID instance
        post_cell = self.postsynaptic_cell._cell # CellType.model instance

        # Get the target region on the cell and the receptor type
        # syn_mech_name = projection.synapse_type.mechanism
        regions_spec, receptor = projection.receptor_type.split(".")
        receptors = receptor.split("+")
        regions = regions_spec.split("2")
        if len(regions) == 2:
            pre_region, post_region = regions
        elif len(regions) == 1:
            pre_region = None
            post_region = regions[0]
        elif len(regions) != 1:
            raise ValueError('Region specification must be "<post>" or "<pre>2<post>"')

        # Plastic synapses use SynapseType.model to store weight adjuster,
        # we use it for the NEURON synapse mechanism
        mech_name = projection.synapse_type.model

        # Ask cell for segment in target region
        num_contacts = getattr(projection.synapse_type, 'num_contacts', 1)
        self.synapses = post_cell.get_synapses(post_region, receptors, num_contacts,
                                               mechanism=mech_name)
        # if used > 0 and not post_cell.allow_synapse_reuse:
        #     raise Exception("No unused synapses on target cell {}".format(type(post_cell)))

        # Get pre-synaptic source
        pre_gid = int(self.presynaptic_cell)
        if (pre_region is not None) and state.cell_has_multiple_sources(pre_gid):
            # TODO: export gid->gid pairs besides connection matrix
            pre_gid = state.query_spkgid(pre_gid, pre_region)

        # Save the actual GIDs (different from ID in case of multicompartment cell)
        self.presynaptic_gid = pre_gid
        self.postsynaptic_gid = int(self.postsynaptic_cell)
            

        # Create NEURON NetCon
        self.ncs = []
        weight = parameters.pop('weight')
        delay = parameters.pop('delay')
        for syn in self.synapses:
            nc = state.parallel_context.gid_connect(pre_gid, syn)
            nc.weight[0] = weight
            nc.delay = delay
            self.ncs.append(nc)

        # PyNN expects attritutes 'synapse' and 'nc'
        # Store first element, should all have same parameter values for Projection.get(...)
        self.synapse = self.synapses[0]
        self.nc = self.ncs[0]

        # if we have a mechanism (e.g. from 9ML) that includes multiple
        # synaptic channels, need to set nc.weight[1] here
        if self.nc.wcnt() > 1 and hasattr(self.postsynaptic_cell._cell, "type"):
            self.nc.weight[1] = self.postsynaptic_cell._cell.type.receptor_types.index(projection.receptor_type)

        # Plastic synapses use SynapseType.model to store weight adjuster,
        # we use it for the synapse model
        if mech_name is not None:
            self._setup_nrn_synapse(projection.synapse_type, parameters)
            # self._setup_plasticity(projection.synapse_type, parameters)


    def _setup_nrn_synapse(self, synapse_type, parameters):
        """
        Set parameters on the NEURON synapse object.

        @param      parameters : dict(str, float)
                    Parameters are already evaluated by PyNN
        """
        # parameters = synapse_type.translate(parameters<ParameterSpace>, copy=copy)
        # parameters.evaluate(simplify=True)
        # for name, value in parameters.items():
        #   setattr(target, name, value)
        synapse_param_names = synapse_type.get_parameter_names()
        for name, value in parameters.items():
            if name in synapse_param_names:
                pinfo = synapse_type.translations[name]
                pname = pinfo['translated_name']
                # setattr(self.synapse, pname, value)
                for synapse in self.synapses:
                    setattr(self.synapse, pname, value)


    def __setattr__(self, name, value):
        """
        Support setting properties of connection's synapse
        through Projection.set(param_name=param_val).

        NOTES
        -----

        - this hides properties of the base class (see Python data model)

        - in original Connector class properties are added for each attribute
          of the associated synapse type  and weight adjuster. See
          See https://github.com/NeuralEnsemble/PyNN/blob/master/pyNN/neuron/simulator.py#L340

        - We catch the property assignments here so we don't have to create
          explicit properties.
        """
        if hasattr(self, 'synapse_type') and (name in self.synapse_type.get_parameter_names()):
            pinfo = self.synapse_type.translations[name]
            pname = pinfo['translated_name']
            # setattr(self.synapse, pname, value)
            for synapse in self.synapses:
                setattr(self.synapse, pname, value)
        elif name in ('weight', 'delay'):
            for nc in self.ncs:
                setattr(nc, name, value)
        else:
            super(ConnectionNrnWrapped, self).__setattr__(name, value)

    def __getattr__(self, name):
        """
        Called as last resort if property not found.
        """
        try:
            stype = super(ConnectionNrnWrapped,self).__getattribute__("synapse_type")
            pinfo = stype.translations[name]
            pname = pinfo['translated_name']
            return getattr(self.synapse, pname)
        except AttributeError:
            return super(ConnectionNrnWrapped,self).__getattribute__(name)

class MultiMechanismConnection(Connection):
    """
    Connection that contacts the post-synaptic cell through multiple synapse
    mechanisms that can be distinct from eachother.
    """

    def __init__(self, projection, pre, post, **parameters):
        """
        Create a new connection.
        """
        #logger.debug("Creating connection from %d to %d, weight %g" % (pre, post, parameters['weight']))
        self.synapse_type = projection.synapse_type
        self.presynaptic_index = pre
        self.postsynaptic_index = post
        self.presynaptic_cell = projection.pre[pre]
        self.postsynaptic_cell = projection.post[post]
        post_cell = self.postsynaptic_cell._cell # CellType.model

        # Get mechanism for each receptor
        self.synapses = {} # dict[str, list[nrn.Synapse]]
        all_synapses = []

        for mech_name, subcell_conn in projection.synapse_type.subcellular_conn.items():
            region, mech_receptors = subcell_conn['receptors'].split(".")
            receptor_list = mech_receptors.split("+")
            synapses = post_cell.get_synapses(region, receptor_list,
                                              subcell_conn['num_contacts'],
                                              mechanism=mech_name,
                                              synapse_reuse=subcell_conn['num_converge'])
            self.synapses[mech_name] = synapses
            all_synapses.extend(synapses)

        # Create NEURON NetCon
        pre_gid = int(self.presynaptic_cell)
        self.ncs = []
        weight = parameters.pop('weight')
        delay = parameters.pop('delay')
        for syn in all_synapses:
            nc = state.parallel_context.gid_connect(pre_gid, syn)
            nc.weight[0] = weight
            nc.delay = delay
            self.ncs.append(nc)

        # Apply synapse parameters
        wrapper_params = self.synapse_type.get_parameter_names()
        for name, value in parameters.items():
            if name in wrapper_params:
                self._set_synapse_parameter(name, value)

        # PyNN expects attributes 'synapse' and 'nc'
        # Store first element, should all have same parameter values for Projection.get(...)
        self.synapse = all_synapses[0]
        self.nc = self.ncs[0]

        # if we have a mechanism (e.g. from 9ML) that includes multiple
        # synaptic channels, need to set nc.weight[1] here
        if self.nc.wcnt() > 1 and hasattr(self.postsynaptic_cell._cell, "type"):
            self.nc.weight[1] = self.postsynaptic_cell._cell.type.receptor_types.index(
                                    projection.receptor_type)


    def _set_synapse_parameter(self, pname, value):
        """
        Apply synapse parameter to all synapses that define that parameter.
        """
        if pname in ('weight', 'delay'):
            for nc in self.ncs:
                setattr(nc, pname, value)
        else:
            param_meta = self.synapse_type.translations[pname]
            match = re.match(r'^([a-zA-Z0-9]+)_(\w+)', param_meta['translated_name'])
            mech_name, nrn_pname = match.groups()
            for synapse in self.synapses[mech_name]:
                setattr(synapse, nrn_pname, value)


    def __setattr__(self, name, value):
        """
        Support setting properties of connection's synapse
        through Projection.set(param_name=param_val).
        """
        if name in ('weight', 'delay'):
            for nc in self.ncs:
                setattr(nc, name, value)
        elif hasattr(self, 'synapse_type') and (name in self.synapse_type.get_parameter_names()):
            self._set_synapse_parameter(name, value)
        else:
            super(MultiMechanismConnection, self).__setattr__(name, value)


    def __getattr__(self, name):
        """
        Called as last resort if property not found.
        """
        try:
            stype = super(MultiMechanismConnection,self).__getattribute__("synapse_type")
            pinfo = stype.translations[name]
            pname = pinfo['translated_name']
            return getattr(self.synapse, pname)
        except AttributeError:
            return super(MultiMechanismConnection,self).__getattribute__(name)
