"""
GABA-B synapse mechanisms with GPCR and exctracellular GABA diffusion.

Mechanisms are loaded into NEURON when module is imported.
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import submodule as submodule_alias
# from .submodule import myclass

import neuron
import os.path

here = os.path.abspath(os.path.dirname(__file__))

# Load mechanisms and scripts
neuron.load_mechanisms(here)
# neuron.h.load_file(os.path.join(here, 'script.hoc'))