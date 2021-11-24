"""
Fujita, Kitano et al. (2011) GPe cell model for use in NEURON.

@author     Lucas Koelman

@date       14/09/2018

Usage
-----

Importing this module will make the cell model available in NEURON by loading
the cell template ('class' in usual OO terms).
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
h.xopen("fujita_createcell.hoc") # instantiates all functions & data structures on Hoc object
os.chdir(prev_cwd)

# Channel mechanisms (key = suffix of mod mechanism) : max conductance parameters
gbar_dict = {
    # Nonspecific channels
    'HCN':      ['gmax'],
    'leak':     ['gmax'],
    # Na channels
    'NaF':      ['gmax'],
    'NaP':      ['gmax'],
    # K-channels
    'Kv2':      ['gmax'],
    'Kv3':      ['gmax'],
    'Kv4f':     ['gmax'],
    'Kv4s':     ['gmax'],
    'KCNQ':     ['gmax'],
    'SK':       ['gmax'],
    # Calcium channels / buffering
    'CaH':      ['gmax'],
}
gleak_name = 'gmax_leak'

# Mechanism and optional parameters that are modified in original model code
mechs_params_dict = {
    # Nonspecific channels
    'HCN':      ['gmax', 'e'], # HCN channel
    'pas':      ['g', 'e'],
    # Na channels
    'NaF':      ['gmax'],
    'NaP':      ['gmax'],
    # K-channels
    'Kv2':      ['gmax'],
    'Kv3':      ['gmax'],
    'Kv4f':     ['gmax'],
    'Kv4s':     ['gmax'],
    'KCNQ':     ['gmax'],
    'SK':       ['gmax'],
    # Calcium channels / buffering
    'Calcium':  [],
    'CaH':      ['gmax', 'e'], # high-voltage-activated calcium channel

}

# All mechanism parameters that are not conductances
mechs_params_nogbar = dict(mechs_params_dict)
for mech, params in mechs_params_nogbar.items():
    for gbar_param in gbar_dict.get(mech, []):
        try:
            params.remove(gbar_param)
        except ValueError:
            pass

# GLOBAL mechanism parameters (assigned using h.param = val)
global_params_list = [
    'ena', 'ek'
]

# List of mechanisms, max conductance params, active conductances
mechs_list = list(mechs_params_dict.keys()) # all mechanisms
gbar_list = [gname+'_'+mech for mech,chans in gbar_dict.items() for gname in chans]
active_gbar_names = [gname for gname in gbar_list if gname != gleak_name]