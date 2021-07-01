"""
Synapse mechanisms for NEURON.

Mechanisms are loaded into NEURON when module is imported.
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import submodule as submodule_alias
# from .submodule import myclass

import neuron
import os, os.path, json
from bgcellmodels.common import nrnmodutil

here = os.path.abspath(os.path.dirname(__file__))

# Load mechanisms in module directory
neuron.load_mechanisms(here)

# Load Hoc scripts in module directory
# neuron.h.load_file(os.path.join(here, 'script.hoc'))

SYN_DB_FILE = os.path.join(here, 'synapse_mechanism_DB.json')
SYN_MECH_DB = None # will contain json file as dict

def rebuild_synapse_database():
    """
    Create JSON file containing information about all .mod files in
    this directory.
    """
    # find all mod files in this directory
    mod_files = [os.path.join(here, f) for f in os.listdir(here) if f.endswith('.mod')]
    syn_db = {}
    for fname in mod_files:
        with open(fname, 'r') as modfile:
            modtext = modfile.read()
        mech_name = nrnmodutil.get_mod_mechname(modtext, filename=False)
        mod_params = nrnmodutil.get_mod_parameters(modtext, filename=False)
        syn_db[mech_name] = {
            'filename': fname,
            'parameters': {}
        }
        for pname, pval, units, comment in mod_params:
            syn_db[mech_name]['parameters'][pname] = {
                'value': pval, 'units': units, 'comment': comment
            }

    # write to JSON file
    with open(SYN_DB_FILE, 'w') as f:
            json.dump(syn_db, f, indent=2)


def load_synapse_database():
    """
    Load database of synaptic mechanisms.

    The database contains the mechanism names, parameters, and mod filenames.

    @return     db : dict
                Synaptic mechanism database
    """
    global SYN_MECH_DB, SYN_DB_FILE
    if not os.path.isfile(SYN_DB_FILE):
        rebuild_synapse_database()
    with open(SYN_DB_FILE, 'r') as f:
        SYN_MECH_DB = json.load(f)


load_synapse_database()