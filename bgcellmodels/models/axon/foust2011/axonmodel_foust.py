"""
Axon model from 'Foust AJ, Yu Y, Popovic M, Zecevic D, McCormick DA (2011) Somatic membrane potential and Kv1 channels control spike repolarization in cortical axon collaterals and presynaptic boutons. J Neurosci 31:15490-8'

@author     Lucas Koelman
@date       24/04/2019
"""

from __future__ import division # float division for int literals, like Hoc
import os.path
import math

import neuron
h = neuron.h

from bgcellmodels.models.axon.axon_base import AxonBuilder

# Load NEURON mechanisms
script_dir = os.path.dirname(__file__)
neuron.load_mechanisms(script_dir)


class AxonFoust2011(AxonBuilder):
    """
    Axon parameters 
    """

    def __init__(self, *args, **kwargs):
        # Set required properties before super call
        self._set_axon_params()
        self._set_compartment_types()
        
        super(AxonFoust2011, self).__init__(*args, **kwargs)


    def _set_axon_params(self):
        """
        Set axon paramters
        """

        self.axon_parameters = {
            # Axon initial segment
            'ais_L': 20,
            'ais_diam': 1.2,
            'ais_nseg': 5,
            'ais_Ra': 150,
            'ais_cm': 0.8,
            # Myelinated compartment
            'myelin_L': 200,
            'myelin_L_0': 80,
            'myelin_diam': 1.4,
            'myelin_Ra': 150,
            'myelin_cm': 0.04,
            'myelin_nseg': 11,
            'myelin_nseg_0': 5,
            # Node of ranvier
            'node_L': 2,
            'node_diam': 1.2,
            'node_nseg': 1,
            'node_Ra': 150,
            'node_cm': 0.8,
            # Collateral section
            'collateral_L': 500,
            'collateral_diam': 0.5,
            'collateral_nseg': 11,
            'collateral_Ra': 150,
            'collateral_cm': 0.8,
            # Other parameters
            'num_axon_compartments': 1,
        }

        # Calculation of extracellular space parameters from McIntyre (2002)
        node_diam = self.axon_parameters['node_diam']
        space_p1 = 0.002

        myelin_diam = self.axon_parameters['myelin_diam']
        space_i = 0.004

        self.axon_parameters['node_Rax_extra'] = (0.7e-8) / (
            math.pi*((((node_diam/2)+space_p1)**2) - ((node_diam/2)**2)))

        self.axon_parameters['myelin_Rax_extra'] = (0.7e-8) / (
            math.pi*((((myelin_diam/2)+space_i)**2) - ((myelin_diam/2)**2)))

        self.axon_parameters['node_gm_extra'] = 1e10
        self.axon_parameters['node_cm_extra'] = 0.0

        myelin_sheet_cm = 0.1       # uF/cm2/lamella membrane
        myelin_sheet_gm = 0.001     # uF/cm2/lamella membrane

        self.axon_parameters['myelin_gm_extra'] = myelin_sheet_gm / (30*2)
        self.axon_parameters['myelin_cm_extra'] = myelin_sheet_cm / (30*2)


    def _set_compartment_types(self):
        """
        Define compartments types that constitute the axon model.

        @pre    electrical parameters have been set

        @pre    morphology parameters have been set
        """
        # Sequence of compartment types that define repeating structure of axon
        self.initial_comp_sequence = ['aisnode', 'aismyelin']
        self.terminal_comp_sequence = ['node', 'collateral']
        self.repeating_comp_sequence = ['node', 'myelin']
        self.collateral_comp_sequence = ['node', 'collateral']

        # Name of compartment type representing node of Ranvier
        self.nodal_compartment_type = 'node'

        self.compartment_defs = {
            # Axon initial segment is like a longer node with modified conductances
            # (higher NaF, lower leak, higher KDR and Kd)
            'aisnode': {
                'mechanisms' : {
                    'pas_Foust': {
                        'g': 3.3e-5,
                    },
                    'NaF_Foust': {
                        'g': 0.4,
                    },
                    'Kv_Foust': {
                        'g': 0.002,
                    },
                    'Kd_Foust': {
                        'g': 0.015,
                    },
                    'extracellular': {
                        'xraxial': self.axon_parameters['node_Rax_extra'],
                        'xg': self.axon_parameters['node_gm_extra'],
                        'xc': self.axon_parameters['node_cm_extra'],
                    }
                }, 
                'passive' : {
                    'Ra': self.axon_parameters['ais_Ra'],
                    'cm': self.axon_parameters['ais_cm'],
                },
                'morphology' : {
                    'diam': self.axon_parameters['ais_diam'],
                    'L': self.axon_parameters['ais_L'],
                    'nseg': self.axon_parameters['ais_nseg']
                },
            },
            # First myelinated section has different length and nseg
            'aismyelin': {
                'mechanisms' : {
                    'pas_Foust': {
                        'g': 0.0,
                    },
                    'NaF_Foust': {
                        'g': 0.001,
                    },
                    'extracellular': {
                        'xraxial': self.axon_parameters['myelin_Rax_extra'],
                        'xg': self.axon_parameters['myelin_gm_extra'],
                        'xc': self.axon_parameters['myelin_cm_extra'],
                    }
                }, 
                'passive' : {
                    'Ra': self.axon_parameters['myelin_Ra'],
                    'cm': self.axon_parameters['myelin_cm'],
                },
                'morphology' : {
                    'diam': self.axon_parameters['myelin_diam'],
                    'L': self.axon_parameters['myelin_L_0'],
                    'nseg': self.axon_parameters['myelin_nseg_0']
                },
            },
            # Nodes of Ranvier
            'node': {
                'mechanisms' : {
                    'pas_Foust': {
                        'g': 0.000033,  # original: 0.02,
                    },
                    'NaF_Foust': {
                        'g': 0.4,       # original: 0.28,
                    },
                    'Kv_Foust': {
                        'g': 0.002,     # original: 0.0003
                    },
                    'Kd_Foust': {
                        'g': 0.015,     # original: 0.0072
                    },
                    'extracellular': {
                        'xraxial': self.axon_parameters['node_Rax_extra'],
                        'xg': self.axon_parameters['node_gm_extra'],
                        'xc': self.axon_parameters['node_cm_extra'],
                    }
                }, 
                'passive' : {
                    'Ra': self.axon_parameters['node_Ra'],
                    'cm': self.axon_parameters['node_cm'],
                },
                'morphology' : {
                    'diam': self.axon_parameters['node_diam'],
                    'L': self.axon_parameters['node_L'],
                    'nseg': self.axon_parameters['node_nseg']
                },
            },
            # Myelinated sections
            'myelin': {
                'mechanisms' : {
                    'pas_Foust': {
                        'g': 0.0,
                    },
                    'NaF_Foust': {
                        'g': 0.001,
                    },
                    'extracellular': {
                        'xraxial': self.axon_parameters['myelin_Rax_extra'],
                        'xg': self.axon_parameters['myelin_gm_extra'],
                        'xc': self.axon_parameters['myelin_cm_extra'],
                    }
                }, 
                'passive' : {
                    'Ra': self.axon_parameters['myelin_Ra'],
                    'cm': self.axon_parameters['myelin_cm'],
                },
                'morphology' : {
                    'diam': self.axon_parameters['myelin_diam'],
                    'L': self.axon_parameters['myelin_L'],
                    'nseg': self.axon_parameters['myelin_nseg']
                },
            },
            # Collateral sections
            'collateral': {
                'mechanisms' : {
                    'pas_Foust': {
                        'g': 0.0000333333,
                    },
                    'NaF_Foust': {
                        'g': 0.133333333,
                    },
                    'Kv_Foust': {
                        'g': 0.0010,
                    },
                    'Kd_Foust': {
                        'g': 0.0060,
                    },
                    'extracellular': {
                        'xraxial': self.axon_parameters['node_Rax_extra'],
                        'xg': self.axon_parameters['node_gm_extra'],
                        'xc': self.axon_parameters['node_cm_extra'],
                    }
                }, 
                'passive' : {
                    'Ra': self.axon_parameters['collateral_Ra'],
                    'cm': self.axon_parameters['collateral_cm'],
                },
                'morphology' : {
                    'diam': self.axon_parameters['collateral_diam'],
                    'L': self.axon_parameters['collateral_L'],
                    'nseg': self.axon_parameters['collateral_nseg']
                },
            }
        }