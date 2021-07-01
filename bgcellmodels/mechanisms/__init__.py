"""
NEURON mechanisms for use in Basal Ganglia models.
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import submodule as submodule_alias
# from .submodule import myclass

import neuron
import os.path

here = os.path.abspath(os.path.dirname(__file__))

# neuron.load_mechanisms(here)
# neuron.h.load_file(os.path.join(here, 'script.hoc'))

def load_synapse_mechanisms():
    neuron.load_mechanisms(os.path.join(here, 'synapses'))