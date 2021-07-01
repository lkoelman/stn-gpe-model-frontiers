"""
Axon model from 'McIntyre CC, Richardson AG, and Grill WM. Modeling the 
excitability of mammalian nerve fibers: influence of afterpotentials on the
recovery cycle. J Neurophysiol 87: 995-1006, 2002'

@author     Lucas Koelman
@date       27/11/2018
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

# Load Hoc functions for cell model
# prev_cwd = os.getcwd()
# os.chdir(script_dir)
# h.xopen("axonmodel.hoc")
# os.chdir(prev_cwd)

PI = math.pi


class AxonMcintyre2002(AxonBuilder):
    """
    See diagram in article Fig. 1 for explanations of compartment types:
    
    - 'node' = node of Ranvier
    
    - 'MYSA' = paranodal myelin attachment segment
             = initial segment connecting unmyelinated node of Ranvier to 
               myelinated internodal part of the axon
    
    - 'FLUT' = paranodal main segment
             = transitional segment between connecting segment 'MYSA' and the
               myelinated internodal part of the axon
    
    - 'STIN' = internodal segment
    """

    def __init__(self, **kwargs):
        # Only pass keyword arguments
        super(AxonMcintyre2002, self).__init__(**kwargs)

        self._set_elec_params()
        self._set_morpho_params()
        self._set_compartment_types()


    def _set_elec_params(self):
        """
        Set electrical paramters
        """
        self.Rho_axial = 0.7e-6          # kilo-Ohm-um
        self.myelin_sheet_cm = 0.1       # uF/cm2/lamella membrane
        self.myelin_sheet_gm = 0.001     # uF/cm2/lamella membrane
        self.Ra = self.Rho_axial * 1e-4


    def _set_morpho_params(self):
        """
        Set Morphological parameters
        """

        self.fiberD = 2                    # fiber diameter
        self.paralength1=3               # MYSA length
        self.nodelength=1
        self.space_p1=0.002  
        self.space_p2=0.004
        self.space_i=0.004

        if (self.fiberD == 1):
            self.axonD=0.8
            self.nodeD=0.7 
            self.paraD1=0.7
            self.paraD2=0.8
            self.deltax=200
            self.paralength2=10
            self.nl=20
        elif (self.fiberD == 2):
            self.axonD=1.6
            self.nodeD=1.4
            self.paraD1=1.4
            self.paraD2=1.6
            self.deltax=200
            self.paralength2=10
            self.nl=30
        else:
            raise ValueError()

        self.Rpn0 = (self.Rho_axial*.01) / (
            PI*((((self.nodeD/2)+self.space_p1)**2) - ((self.nodeD/2)**2)))
        self.Rpn1 = (self.Rho_axial*.01) / (
            PI*((((self.paraD1/2)+self.space_p1)**2) - ((self.paraD1/2)**2)))
        self.Rpn2 = (self.Rho_axial*.01) / (
            PI*((((self.paraD2/2)+self.space_p2)**2) - ((self.paraD2/2)**2)))
        self.Rpx  = (self.Rho_axial*.01) / (
            PI*((((self.axonD/2)+self.space_i)**2) - ((self.axonD/2)**2)))


    def _set_compartment_types(self):
        """
        Define compartments types that constitute the axon model.

        @pre    electrical parameters have been set

        @pre    morphology parameters have been set
        """
        # Sequence of compartment types that define repeating structure of axon
        self.initial_comp_sequence = []
        self.repeating_comp_sequence = ['node', 'MYSA', 'FLUT'] + \
                                       ['STIN', 'STIN', 'STIN'] + \
                                       ['FLUT', 'MYSA']

        # Name of compartment type representing node of Ranvier
        self.nodal_compartment_type = 'node'

        self.compartment_defs = {
            # Nodes of Ranvier
            'node': {
                'mechanisms' : {
                    'axnode75': {
                        'gnabar': 2.0,
                        'gnapbar': 0.05,
                        'gkbar': 0.07,
                        'gl': 0.005,
                        'ek': -85,
                        'ena': 55,
                        'el': -60,
                        'vshift': 15,
                        'vtraub': -80,
                    },
                    'extracellular': {
                        'xraxial': self.Rpn0,
                        'xg': 1e10,
                        'xc': 0.0,
                    }
                }, 
                'passive' : {
                    'Ra': self.Ra,
                    'cm': 2.0,
                },
                'morphology' : {
                    'diam': 1.4,
                    'L': 1.0,
                },
            },
            'MYSA': {
                'mechanisms' : {
                    'pas' : {},
                    'extracellular': {
                        'xraxial': self.Rpn1,
                        'xg': self.myelin_sheet_gm / (self.nl*2),
                        'xc': self.myelin_sheet_cm / (self.nl*2),
                    }
                }, 
                'passive' : {
                    'g_pas': 0.0001,
                    'e_pas': -65,
                    'Ra': self.Ra,
                    'cm': 2.0,
                },
                'morphology' : {
                    'diam': 1.4,
                    'L': 3.0,
                },
            },
            'FLUT': {
                'mechanisms' : {
                    'parak75': {
                        'gkbar': 0.02,
                        'ek': -85,
                        'vshift': 15,
                    },
                    'pas' : {},
                    'extracellular': {
                        'xraxial': self.Rpn2,
                        'xg': self.myelin_sheet_gm / (self.nl*2),
                        'xc': self.myelin_sheet_cm / (self.nl*2),
                    }
                }, 
                'passive' : {
                    'g_pas': 0.0001,
                    'e_pas': -60,
                    'Ra': self.Ra,
                    'cm': 2.0,
                },
                'morphology' : {
                    'diam': 1.6,
                    'L' : 10.0,
                },
            },
            'STIN': {
                'mechanisms' : {
                    'pas' : {},
                    'extracellular': {
                        'xraxial': self.Rpx,
                        'xg': self.myelin_sheet_gm / (self.nl*2),
                        'xc': self.myelin_sheet_cm / (self.nl*2),
                    }
                }, 
                'passive' : {
                    'g_pas': 0.0001,
                    'e_pas': -65,
                    'Ra': self.Ra,
                    'cm': 2.0,
                },
                'morphology' : {
                    'diam': 1.6,
                    'L': 29.0,
                },
            },
        }