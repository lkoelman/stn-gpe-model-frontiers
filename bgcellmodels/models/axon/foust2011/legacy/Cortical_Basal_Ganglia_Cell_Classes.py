# -*- coding: utf-8 -*-
""" 
------------------------------------------------------------------------------------
Cortical Axon and Basal Ganglia Neurons: file containing classes for defining network neurons
------------------------------------------------------------------------------------

                            Model References 
------------------------------------------------------------------------------------

Cortical Pyramidal Cell Axon: 
Foust, A.J., Yu, Y., Popovic, M., Zecevic, D. and McCormick, D.A., 
2011. "Somatic membrane potential and Kv1 channels control spike 
repolarization in cortical axon collaterals and presynaptic boutons." 
Journal of Neuroscience, 31(43), pp.15490-15498.


STN Neurons:
Otsuka, T., Abe, T., Tsukagawa, T. and Song, W.J., 2004. 
"Conductance-based model of the voltage-dependent generation 
of a plateau potential in subthalamic neurons."
Journal of neurophysiology, 92(1), pp.255-264.

GP Neurons:




Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18

Edits: 16-01-19: Created classes for cell models so they can be 
                 utilized in PyNN.

Created on Tues Jan 15 12:51:26 2019

"""

from math import pi
from neuron import h
from nrnutils import Mechanism, Section
from pyNN.neuron import NativeCellType
from pyNN.parameters import Sequence
import numpy as np

try:
    reduce
except NameError:
    from functools import reduce


def _new_property(obj_hierarchy, attr_name):
    """
    Returns a new property, mapping attr_name to obj_hierarchy.attr_name.

    For example, suppose that an object of class A has an attribute b which
    itself has an attribute c which itself has an attribute d. Then placing
            e = _new_property('b.c', 'd')
    in the class definition of A makes A.e an alias for A.b.c.d
    """

    def set(self, value):
        obj = reduce(getattr, [self] + obj_hierarchy.split('.'))
        setattr(obj, attr_name, value)

    def get(self):
        obj = reduce(getattr, [self] + obj_hierarchy.split('.'))
        return getattr(obj, attr_name)
    return property(fset=set, fget=get)


class Cortical_Neuron_Axon(object):

    def __init__(self, **parameters):

        # Multiplication by 1e-4 is to convert units from Eleanor's model
        # (pS/um2 -> S/cm2)
        self.ais = Section(
            L=parameters['ais_L'],
            diam=parameters['ais_diam'],
            nseg=parameters['ais_nseg'],
            Ra=parameters['ais_Ra'],
            cm=parameters['ais_cm'],
            mechanisms=(
                Mechanism('cortical_axon_i_leak', g_l=3.3e-5),
                Mechanism('cortical_axon_i_na', g_Na=4000e-4),
                Mechanism('cortical_axon_i_kv', g_Kv=20e-4),
                Mechanism('cortical_axon_i_kd', g_Kd=0.015)))

        # Use loop to create myelin and node sections of axon
        self.myelin = []
        self.node = []
        for i in np.arange(parameters['num_axon_compartments']):
            if i == 0:
                self.myelin.append(
                    Section(L=parameters['myelin_L_0'],
                        diam=parameters['myelin_diam'],
                        nseg=5,
                        Ra=parameters['myelin_Ra'],
                        cm=parameters['myelin_cm'],
                        mechanisms=(Mechanism('cortical_axon_i_leak', g_l=0),
                                    Mechanism('cortical_axon_i_na', g_Na=0.001)),
                        parent=self.ais))
            else:
                self.myelin.append(
                    Section(
                        L=parameters['myelin_L'],
                        diam=parameters['myelin_diam'],
                        nseg=11,
                        Ra=parameters['myelin_Ra'],
                        cm=parameters['myelin_cm'],
                        mechanisms=(Mechanism('cortical_axon_i_leak', g_l=0),
                                    Mechanism('cortical_axon_i_na', g_Na=10e-4)),
                        parent=self.node[i - 1]))

            self.node.append(
                Section(
                    L=parameters['node_L'],
                    diam=parameters['node_diam'],
                    nseg=parameters['node_nseg'],
                    Ra=parameters['node_Ra'],
                    cm=parameters['node_cm'],
                    mechanisms=(Mechanism('cortical_axon_i_leak', g_l=0.02),
                                Mechanism('cortical_axon_i_na', g_Na=2800e-4),
                                Mechanism('cortical_axon_i_kv', g_Kv=5e-4),
                                Mechanism('cortical_axon_i_kd', g_Kd=0.0072)),
                    parent=self.myelin[i]))

        self.collateral = Section(
            L=parameters['collateral_L'],
            diam=parameters['collateral_diam'],
            nseg=parameters['collateral_nseg'],
            Ra=parameters['collateral_Ra'],
            cm=parameters['collateral_cm'],
            mechanisms=(Mechanism('cortical_axon_i_leak'),
                        Mechanism('cortical_axon_i_na', g_Na=1333.33333e-4),
                        Mechanism('cortical_axon_i_kv', g_Kv=10e-4),
                        Mechanism('cortical_axon_i_kd', g_Kd=0.006)),
            parent=self.node[-1])

        self.end_node = self.node[0]
        self.end_myelin = self.myelin[0]

        # Add extracellular (and xtra) mechanism to collateral
        self.collateral.insert('extracellular')
        self.collateral.insert('xtra')

        # Assign rx values to the segments rx_xtra
        for seg in self.collateral:
            seg.xtra.rx = seg.x * 3e-1

        # Maybe above was wrong when setting pointers
        for seg in self.collateral:
            h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)
            h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)

        # insert synaptic noise
        self.noise = h.SynNoise(0.5, sec=self.ais)
        self.noise.f0 = 0
        self.noise.f1 = 0.3

        # Add AMPA and GABAa synapses to the cell, i.e. add to the ais section
        self.AMPA = h.AMPA_S(0.5, sec=self.ais)
        self.GABAa = h.GABAa_S(0.5, sec=self.ais)

        # needed for PyNN
        self.source_section = self.collateral
        self.source = self.collateral(0.5)._ref_v
        # Needed to clear the simulator
        self.rec = h.NetCon(self.source, None, sec=self.collateral)
        self.spike_times = h.Vector(0)
        self.traces = {}
        self.recording_time = False
        self.parameter_names = ()

    def memb_init(self):
        for seg in self.ais:
            seg.v = self.v_init

    def _set_collateral_rx(self, sequence_values):
        rx_values = sequence_values.value
        for ii, seg in enumerate(self.collateral):
            seg.xtra.rx = rx_values[ii]

    def _get_collateral_rx(self):
        print("Getter Working!")
        rx_values = np.zeros((1, self.collateral.nseg))
        for i, seg in enumerate(self.collateral):
            rx_values[0, i] = seg.xtra.rx
        print(Sequence(rx_values.flatten()))
    collateral_rx = property(fget=_get_collateral_rx, fset=_set_collateral_rx)


class Cortical_Neuron_Axon_Type(NativeCellType):
    default_parameters = {
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
        'num_axon_compartments': 1
    }

    # Define initial vector of transfer resistances for the collateral segments
    initial_collateral_rx = np.zeros(
        (1, default_parameters['collateral_nseg'])).flatten()
    initial_collateral_rx_Sequence = Sequence(initial_collateral_rx)
    default_parameters['collateral_rx'] = initial_collateral_rx_Sequence

    default_initial_values = {'v': -68.0}
    recordable = ['collateral(0.5).v', 'collateral(0.5).i_membrane_',
                  'ais(0.5).v', 'end_node(0.5).v', 'end_myelin(0.5).v']
    units = {'collateral(0.5).v': 'mV', 'collateral(0.5).i_membrane_': 'nA',
             'ais(0.5).v': 'mV', 'end_node(0.5).v': 'mV', 'end_myelin(0.5).v': 'mV'}
    receptor_types = ['AMPA', 'GABAa']
    model = Cortical_Neuron_Axon


class STN_Neuron(object):

    def __init__(self, **parameters):

        # Create single compartment Otsuka STN cell section, i.e. soma section
        # Note: 100um2 area so point process currents (nA) are equivalent to
        # density currents (mA/cm2)
        self.soma = Section(L=parameters['L'], diam=parameters['diam'], nseg=parameters['nseg'], Ra=parameters['Ra'], cm=parameters['cm'],
                            mechanisms=(Mechanism('stn_i_leak'), Mechanism('stn_cai_dynamics'), Mechanism('stn_i_a'), Mechanism('stn_i_ca_k'),
                                        Mechanism('stn_i_k'), Mechanism('stn_i_l'), Mechanism('stn_i_na'), Mechanism('stn_i_t')))

        # Add bias current to neuron model - current density is in terms of original model paper, uA/cm2
        # insert current source
        self.stim = h.IClamp(0.5, sec=self.soma)
        self.stim.delay = 0
        self.stim.dur = 1e12
        # bias current density = uA/cm2
        self.stim.amp = parameters[
            'bias_current_density'] * (self.area() * (1e-8)) * (1e3)
        # area = pi*10000um2 -> 1um2 = 1e-8cm2 -> area = pi*10000*(1e-8)cm2 (conversion factor)
        # (1e3) is conversion factor so uA -> nA

        # insert synaptic noise
        self.noise = h.SynNoise(0.5, sec=self.soma)
        self.noise.f0 = 0
        self.noise.f1 = 0.3

        # Add AMPA and GABAa synapses to the cell, i.e. add to the soma section
        self.AMPA = h.AMPA_S(0.5, sec=self.soma)
        self.GABAa = h.GABAa_S(0.5, sec=self.soma)

        # needed for PyNN
        self.source_section = self.soma
        self.source = self.soma(0.5)._ref_v
        # Needed to clear the simulator
        self.rec = h.NetCon(self.source, None, sec=self.soma)
        self.spike_times = h.Vector(0)
        self.parameter_names = ('L', 'diam', 'nseg',
                                'Ra', 'cm', 'bias_current_density')
        self.traces = {}
        self.recording_time = False

    L = _new_property('soma', 'L')
    diam = _new_property('soma', 'diam')
    nseg = _new_property('soma', 'nseg')
    Ra = _new_property('soma', 'Ra')
    cm = _new_property('soma', 'cm')
    bias_current_amp = _new_property('stim', 'amp')

    def area(self):
        """Membrane area in µm²"""
        return pi * self.soma.L * self.soma.diam                                            # pi*L*diam -> um*um -> 1um^2 = 1e-8cm^2

    def memb_init(self):
        for seg in self.soma:
            seg.v = self.v_init


class STN_Neuron_Type(NativeCellType):
    default_parameters = {'L': 60, 'diam': 60, 'nseg': 1,
                          'Ra': 150, 'cm': 1, 'bias_current_density': 0.25}
    default_initial_values = {'v': -68.0}
    recordable = ['soma(0.5).v', 'AMPA.i', 'GABAa.i']
    units = {'soma(0.5).v': 'mV', 'AMPA.i': 'nA', 'GABAa.i': 'nA'}
    receptor_types = ['AMPA', 'GABAa']
    model = STN_Neuron


class GP_Neuron(object):

    def __init__(self, **parameters):

        # Create single compartment Rubin and Terman GP cell section, i.e. soma
        # section
        self.soma = Section(L=parameters['L'], diam=parameters['diam'], nseg=parameters['nseg'], Ra=parameters['Ra'], cm=parameters['cm'],
                            mechanisms=(Mechanism('gp_i_leak'), Mechanism('gp_i_naf'), Mechanism('gp_i_nap'), Mechanism('gp_i_kv2'),
                                        Mechanism('gp_i_kv3'), Mechanism('gp_i_kv4f'), Mechanism(
                                            'gp_i_kv4s'), Mechanism('gp_i_kcnq'),
                                        Mechanism('gp_i_cah'), Mechanism('gp_i_hcn'), Mechanism('gp_i_sk'), Mechanism('gp_cai_dynamics')))

        # insert current source
        self.stim = h.IClamp(0.5, sec=self.soma)
        self.stim.delay = 0
        self.stim.dur = 1e12
        # (0.001 or 1e-3) is conversion factor so pA -> nA
        self.stim.amp = parameters[
            'bias_current_density'] * (self.area()) * (0.001)

        # Add AMPA and GABAa synapses to the cell, i.e. add to the soma section
        self.AMPA = h.AMPA_S(0.5, sec=self.soma)
        self.GABAa = h.GABAa_S(0.5, sec=self.soma)

        # needed for PyNN
        self.source_section = self.soma
        self.source = self.soma(0.5)._ref_v
        # Needed to clear the simulator
        self.rec = h.NetCon(self.source, None, sec=self.soma)
        self.spike_times = h.Vector(0)
        self.parameter_names = ('L', 'diam', 'nseg',
                                'Ra', 'cm', 'bias_current_density')
        self.traces = {}
        self.recording_time = False

    L = _new_property('soma', 'L')
    diam = _new_property('soma', 'diam')
    nseg = _new_property('soma', 'nseg')
    Ra = _new_property('soma', 'Ra')
    cm = _new_property('soma', 'cm')
    bias_current_density = _new_property('stim', 'amp')

    def area(self):
        """Membrane area in µm²"""
        return pi * self.soma.L * self.soma.diam

    def memb_init(self):
        for seg in self.soma:
            seg.v = self.v_init


class GP_Neuron_Type(NativeCellType):
    default_parameters = {'L': 10, 'diam': 10, 'nseg': 1,
                          'Ra': 150, 'cm': 2.4, 'bias_current_density': 0}
    default_initial_values = {'v': -68.0}
    recordable = ['soma(0.5).v']
    units = {'soma(0.5).v': 'mV'}
    receptor_types = ['AMPA', 'GABAa']
    model = GP_Neuron
