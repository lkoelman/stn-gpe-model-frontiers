"""
Create synaptic connections using custom NEURON synapse models
and get their parameters from centralized parameter databvase.

@author     Lucas Koelman

@date       21/03/2018


Example
-------

See module 'connection' for usage and developer notes.

DEVNOTES
--------

- The synapse base classes are, from abstract to concrete:

    pyNN/models.py
        - BaseModelType
        - BaseSynapseType

    pyNN/standardmodels/__init__.py
        - StandardModelType
        - StandardSynapseType

    pyNN/standardmodels/synapses.py
        - StaticSynapse

    pyNN/neuron/standardmodels/synapses.py
        - BaseSynapse
        - StaticSynapse

- parameter names must be eval()-able so no special chars
"""

import pyNN.neuron
from pyNN.neuron.simulator import state
from pyNN.standardmodels import synapses, build_translations
from . import connection as custom_conn
from bgcellmodels.mechanisms import synapses as synmech


class NativeSynapse(pyNN.neuron.StaticSynapse):
    """
    A NEURON native synapse model with an explicit MOD mechanism
    and mechanism parameters.

    You can specify any mechanism name as the 'mechanism' argument,
    and its parameters as 'mechanism_parameters'

    @see        pyNN.standardmodels.synapses.StaticSynapse
    """
    connection_type = custom_conn.NativeSynToRegion

    # default_parameters = {
    #     'weight': 0.0,
    #     'delay': None,
    #     'mechanism': 'Exp2Syn',
    #     'mechanism_parameters': None,
    # }

    def __init__(self, **parameters):
        """
        Make new NEURON native synapse


        @param      mechanism : str
                    NEURON synapse mechanism name


        @param      **parameters : (keyword arguments)
                    
                    Only 'weight' 'delay' and 'mechanism_parameters'
                    are recognized keywords.


        @param      mechanism_parameters : dict(str, float)
                    
                    Parameters of neuron mechanism and NetCon, in format
                    "syn:attribute[index]" and "netcon:attribute[index]".
        """
        # Don't pass native mechanism parameters
        self.mechanism = parameters.pop('mechanism')
        self.mechanism_parameters = parameters.pop('mechanism_parameters', {})
        self.multi_synapse_rule = parameters.pop('multi_synapse_rule', 1)
        # self.physiological_parameters = parameters.pop('physiological_parameters', {})

        # Convert distance-based string expressions to callable functions
        converted_expressions = {}
        for param_spec, value_spec in self.mechanism_parameters.items():
            if not isinstance(value_spec, str):
                continue
            d_expression = value_spec
            try:
                # Check for singularities at large and small distances
                d = 0; assert 0 <= eval(d_expression), eval(d_expression)
                d = 1e12; assert 0 <= eval(d_expression), eval(d_expression)
            except ZeroDivisionError as err:
                raise ZeroDivisionError("Error in the distance expression %s. %s" % (d_expression, err))
            converted_expressions[param_spec]= eval("lambda d: {}".format(d_expression))
        self.mechanism_parameters.update(converted_expressions)
        
        # Insert dummy variables so pyNN doesn't complain
        parameters.setdefault('delay', 1.0)
        parameters.setdefault('weight', 0.0)
        super(NativeSynapse, self).__init__(**parameters)


class NativeMultiSynapse(pyNN.neuron.BaseSynapse, synapses.StaticSynapse):
    """
    Represents one axon that synapses onto a post-synaptic cell via
    multiple synapses that may have distinct synaptic mechanisms.

    The synaptic mechanisms are specified as NEURON mechanism names.
    """

    connection_type = custom_conn.MultiMechanismConnection

    # Map mechanism name to dict describing subcellular connectivity:
    # - 'receptors' : receptor/region specifier
    # - 'num_contacts': number of contacts per afferent with this mechanism
    # - 'num_converge': how many afferents should converge onto a single synapse of this mechanism
    subcellular_conn = {}

    # PyNN internal name to NEURON name
    translations = None # set in __init__

    default_parameters = {
        'weight':     1.0,
        'delay':      0.5,
    }

    def __init__(self, **kwargs):
        """

        Arguments
        ---------

        @param  mechanisms_receptors : dict[str,str]
                Maps NEURON mechanism name onto the receptor/region specification
                for the post-synaptic cell.
        """
        # TODO: provide dict "contact_types" that has entries mechanism, num_contacts, ...
        mechs_receptors = kwargs.pop('mechanisms_receptors')
        num_contacts = kwargs.pop('num_contacts', 1)
        translation_pairs = [('weight', 'weight'), ('delay', 'delay')]

        # Build default parameters and translations for mechanism
        for mechname, receptor in mechs_receptors.items():
            self.subcellular_conn[mechname] = {
                'receptors': receptor,
                'num_contacts': num_contacts,
                'num_converge': 1,
            }
            for pname, pinfo in synmech.SYN_MECH_DB[mechname]['parameters'].items():
                prefixed_pname = mechname + '_' + pname
                translation_pairs.append((prefixed_pname, prefixed_pname))
                self.default_parameters[prefixed_pname] = pinfo['value']

        self.translations = build_translations(*translation_pairs)
        super(NativeMultiSynapse, self).__init__(**kwargs)


    def _get_minimum_delay(self):
        return state.min_delay


class GluSynapse(pyNN.neuron.BaseSynapse, synapses.StaticSynapse):
    """
    Wrapper for NEURON GLUsyn mechanism defined in .mod file
    """

    connection_type = custom_conn.ConnectionNrnWrapped
    model = 'GLUsyn'

    # PyNN internal name to NEURON name
    translations = build_translations(
        # NetCon parameters
        ('weight', 'weight'),
        ('delay', 'delay'),
        # Conductance time course
        ('gmax_AMPA', 'gmax_AMPA'),     # Weight conversion factor (from nS to uS)
        ('gmax_NMDA', 'gmax_NMDA'),     # Weight conversion factor (from nS to uS)
        ('tau_r_AMPA', 'tau_r_AMPA'),   # Dual-exponential conductance profile
        ('tau_d_AMPA', 'tau_d_AMPA'),   # IMPORTANT: tau_r < tau_d
        ('tau_r_NMDA', 'tau_r_NMDA'),   # Dual-exponential conductance profile
        ('tau_d_NMDA', 'tau_d_NMDA'),    # IMPORTANT: tau_r < tau_d
        # Short-term Depression/Facilitation
        ('tau_rec', 'tau_rec'),         # time constant of recovery from depression
        ('tau_facil', 'tau_facil'),     # time constant of facilitation
        ('U1', 'U1'),                   # baseline release probability
        # Magnesium block for NMDA
        ('e', 'e'),                     # AMPA and NMDA reversal potential
        ('mg', 'mg'),                   # Initial concentration of mg2+
    )

    default_parameters = {
        'weight':       1.0,
        'delay':        0.5,
        'tau_r_AMPA':   0.2,
        'tau_d_AMPA':   1.7,
        'tau_r_NMDA':   0.29,
        'tau_d_NMDA':   43.0,
        'e':            0.0,
        'mg':           1.0,
        'gmax_AMPA':    0.001,
        'gmax_NMDA':    0.001,
        'tau_rec':      200.0,
        'tau_facil':    200.0,
        'U1':           0.5,
    }

    def _get_minimum_delay(self):
        return state.min_delay


class GabaSynapse(pyNN.neuron.BaseSynapse, synapses.StaticSynapse):
    """
    Wrapper for NEURON GLUsyn mechanism defined in .mod file
    """

    connection_type = custom_conn.ConnectionNrnWrapped
    model = 'GABAsyn' # defined in GABAsyn.mod

    # PyNN internal name to NEURON name
    translations = build_translations(
        # NetCon parameters
        ('weight', 'weight'),
        ('delay', 'delay'),
        # Conductance time course
        ('gmax_GABAA', 'gmax_GABAA'),     # Weight conversion factor (from nS to uS)
        ('gmax_GABAB', 'gmax_GABAB'),     # Weight conversion factor (from nS to uS)
        ('tau_r_GABAA', 'tau_r_GABAA'),   # Dual-exponential conductance profile
        ('tau_d_GABAA', 'tau_d_GABAA'),   # IMPORTANT: tau_r < tau_d
        ('tau_r_GABAB', 'tau_r_GABAB'),   # Dual-exponential conductance profile
        ('tau_d_GABAB', 'tau_d_GABAB'),    # IMPORTANT: tau_r < tau_d
        # Short-term Depression/Facilitation
        ('tau_rec', 'tau_rec'),         # time constant of recovery from depression
        ('tau_facil', 'tau_facil'),     # time constant of facilitation
        ('U1', 'U1'),                   # baseline release probability
        # Reversal potentials
        ('Erev_GABAA', 'Erev_GABAA'),                     # GABAA and GABAB reversal potential
        ('Erev_GABAB', 'Erev_GABAB'),                   # Initial concentration of mg2+
        # TODO: include GABA-B signaling cascade if necessary
    )

    default_parameters = {
        'weight':     1.0,
        'delay':      0.5,
        'tau_r_GABAA':   0.2,
        'tau_d_GABAA':   1.7,
        'tau_r_GABAB':   0.2,
        'tau_d_GABAB':   1.7,
        'Erev_GABAA':   -80.0,
        'Erev_GABAB':   -95.0,
        'gmax_GABAA':    0.001,
        'gmax_GABAB':    0.001,
        'tau_rec':      200.0,
        'tau_facil':    200.0,
        'U1':           0.5,
    }

    def _get_minimum_delay(self):
        return state.min_delay


class GabaSynTmHill(pyNN.neuron.BaseSynapse, synapses.StaticSynapse):
    """
    Tsodyks-Markram synapse for GABA-A and GABA-B receptors with GABA-B
    conductance expressed as hill function applied to Tsodyks-Markram
    conductance variable.
    
    @see    mechanism GABAsyn2.mod
    """

    connection_type = custom_conn.ConnectionNrnWrapped
    model = 'GABAsyn2' # defined in GABAsynTmGprot2.mod

    # PyNN internal name to NEURON name
    translations = build_translations(
        # NetCon parameters
        ('weight', 'weight'),
        ('delay', 'delay'),
        # Conductance time course
        ('gmax_GABAA', 'gmax_GABAA'),     # Weight conversion factor (from nS to uS)
        ('gmax_GABAB', 'gmax_GABAB'),     # Weight conversion factor (from nS to uS)
        ('tau_r_GABAA', 'tau_r_GABAA'),   # Dual-exponential conductance profile
        ('tau_d_GABAA', 'tau_d_GABAA'),   # IMPORTANT: tau_r < tau_d
        ('tau_r_GABAB', 'tau_r_GABAB'),   # Dual-exponential conductance profile
        ('tau_d_GABAB', 'tau_d_GABAB'),    # IMPORTANT: tau_r < tau_d
        # Short-term Depression/Facilitation
        ('tau_rec', 'tau_rec'),         # time constant of recovery from depression
        ('tau_facil', 'tau_facil'),     # time constant of facilitation
        ('U1', 'U1'),                   # baseline release probability
        # Reversal potentials
        ('Erev_GABAA', 'Erev_GABAA'),                     # GABAA and GABAB reversal potential
        ('Erev_GABAB', 'Erev_GABAB'),                   # Initial concentration of mg2+
        # GABA-B signaling cascade
        ('K3', 'K3'),
        ('K4', 'K4'),
        ('KD', 'KD'),
        ('n', 'n'),
    )

    default_parameters = {
        'weight':       1.0,
        'delay':        0.5,
        'tau_r_GABAA':  0.2,
        'tau_d_GABAA':  1.7,
        'tau_r_GABAB':  0.2,
        'tau_d_GABAB':  1.7,
        'Erev_GABAA':   -80.0,
        'Erev_GABAB':   -95.0,
        'gmax_GABAA':   0.001,
        'gmax_GABAB':   0.001,
        'tau_rec':      200.0,
        'tau_facil':    200.0,
        'U1':           0.5,
        'K3':           0.098,
        'K4':           0.033,
        'KD':           100.0,
        'n':            4,
    }

    def _get_minimum_delay(self):
        return state.min_delay


class GABAAsynTM(pyNN.neuron.BaseSynapse, synapses.StaticSynapse):
    """
    Wrapper for NEURON GLUsyn mechanism defined in .mod file
    """

    connection_type = custom_conn.ConnectionNrnWrapped
    model = 'GABAAsynTM' # defined in GABAsyn.mod

    # PyNN internal name to NEURON name
    translations = build_translations(
        # NetCon parameters
        ('weight', 'weight'),
        ('delay', 'delay'),
        # Conductance time course
        ('gmax_GABAA', 'gmax_GABAA'),     # Weight conversion factor (from nS to uS)
        ('tau_r_GABAA', 'tau_r_GABAA'),   # Dual-exponential conductance profile
        ('tau_d_GABAA', 'tau_d_GABAA'),   # IMPORTANT: tau_r < tau_d
        # Short-term Depression/Facilitation
        ('tau_rec', 'tau_rec'),         # time constant of recovery from depression
        ('tau_facil', 'tau_facil'),     # time constant of facilitation
        ('U1', 'U1'),                   # baseline release probability
        # Reversal potentials
        ('Erev_GABAA', 'Erev_GABAA'),
    )

    default_parameters = {
        'weight':     1.0,
        'delay':      0.5,
        'tau_r_GABAA':   0.2,
        'tau_d_GABAA':   1.7,
        'Erev_GABAA':   -80.0,
        'gmax_GABAA':    0.001,
        'tau_rec':      200.0,
        'tau_facil':    200.0,
        'U1':           0.5,
    }

    def _get_minimum_delay(self):
        return state.min_delay


class GabaSynTm2(pyNN.neuron.BaseSynapse, synapses.StaticSynapse):
    """
    Wrapper for NEURON GABAsynTM2 mechanism defined in .mod file
    """

    connection_type = custom_conn.ConnectionNrnWrapped
    model = 'GABAsynTM2'

    # PyNN internal name to NEURON name
    translations = build_translations(
        # NetCon parameters
        ('weight', 'weight'),
        ('delay', 'delay'),
        # Conductance time course
        ('gmax_GABAA', 'gmax_GABAA'),
        ('tau_r_GABAA', 'tau_r_GABAA'),
        ('tau_d_GABAA', 'tau_d_GABAA'),
        ('gmax_GABAB', 'gmax_GABAB'),
        ('tau_r_GABAB', 'tau_r_GABAB'),
        ('tau_d_GABAB', 'tau_d_GABAB'),
        # Short-term Depression/Facilitation
        ('tau_rec_A', 'tau_rec_A'),         # time constant of recovery from depression
        ('tau_facil_A', 'tau_facil_A'),     # time constant of facilitation
        ('U1_A', 'U1_A'),                   # baseline release probability
        ('tau_rec_B', 'tau_rec_B'),         # time constant of recovery from depression
        ('tau_facil_B', 'tau_facil_B'),     # time constant of facilitation
        ('U1_B', 'U1_B'),                   # baseline release probability
        # Reversal potentials
        ('Erev_GABAA', 'Erev_GABAA'),
        ('Erev_GABAB', 'Erev_GABAB'),
    )

    default_parameters = {
        'weight':     1.0,
        'delay':      0.5,
        'tau_r_GABAA':   0.2,
        'tau_d_GABAA':   1.7,
        'tau_rec_A':      200.0,
        'tau_facil_A':    200.0,
        'U1_A':           0.5,
        'Erev_GABAA':   -80.0,
        'gmax_GABAA':    0.001,
        'tau_r_GABAB':   5.0,
        'tau_d_GABAB':   25.0,
        'tau_rec_B':      5.0,
        'tau_facil_B':    100.0,
        'U1_B':           0.05,
        'Erev_GABAB':   -95.0,
        'gmax_GABAB':    0.001,
    }

    def _get_minimum_delay(self):
        return state.min_delay


class GabaABMultiMech(pyNN.neuron.BaseSynapse, synapses.StaticSynapse):
    """
    One synapse class that create two distinct synapse objects for
    the GABA-A and GABA-B receptors.
    """

    connection_type = custom_conn.MultiMechanismConnection

    # Map receptors to synapse mechanisms
    mechanisms = {
        'GABAA': 'GABAAsynTM',
        'GABAB': 'GABAAsynTM',
    }

    # For each receptor: how many NetCons converge on one synape implementing
    # that receptor.
    num_converging = {
        'GABAA': 1,
        'GABAB': 3,
    }

    # PyNN internal name to NEURON name
    translations = build_translations(
        # NetCon parameters
        ('weight', 'weight'),
        ('delay', 'delay'),
        # GABA-A receptor
        ('Erev_GABAA',  'GABAA:Erev_GABAA'),
        ('gmax_GABAA',  'GABAA:gmax_GABAA'),
        ('tau_r_GABAA', 'GABAA:tau_r'),
        ('tau_d_GABAA', 'GABAA:tau_d'),
        ('U1_A',        'GABAA:U1'),
        ('tau_rec_A',   'GABAA:tau_rec'),
        ('tau_facil_A', 'GABAA:tau_facil'),
        # GABA-B receptor
        ('Erev_GABAB',  'GABAB:Erev_GABAA'),
        ('gmax_GABAB',  'GABAB:gmax_GABAA'),
        ('tau_r_GABAB', 'GABAB:tau_r'),
        ('tau_d_GABAB', 'GABAB:tau_d'),
        ('U1_B',        'GABAB:U1'),
        ('tau_rec_B',   'GABAB:tau_rec'),
        ('tau_facil_B', 'GABAB:tau_facil'),
    )

    default_parameters = {
        'weight':     1.0,
        'delay':      0.5,
        'tau_r_GABAA':   0.2,
        'tau_d_GABAA':   1.7,
        'tau_rec_A':      200.0,
        'tau_facil_A':    200.0,
        'U1_A':           0.5,
        'Erev_GABAA':   -80.0,
        'gmax_GABAA':    0.001,
        'tau_r_GABAB':   5.0,
        'tau_d_GABAB':   25.0,
        'tau_rec_B':      5.0,
        'tau_facil_B':    100.0,
        'U1_B':           0.05,
        'Erev_GABAB':   -95.0,
        'gmax_GABAB':    0.001,
    }

    def _get_minimum_delay(self):
        return state.min_delay
