"""
Noise mechanisms for NEURON.

Mechanisms are loaded into NEURON when module is imported.
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import submodule as submodule_alias
# from .submodule import myclass

import neuron
import os

here = os.path.abspath(os.path.dirname(__file__))

# Load mechanisms in module directory
neuron.load_mechanisms(here)

# List Hoc files that should be loaded
hoc_files = [
    f for f in os.listdir(here) if f.endswith('.hoc')
]

# Load Hoc files
for hoc_file in hoc_files:
    neuron.h.load_file(os.path.join(here, hoc_file))

# from . import xtra_utils as xtra