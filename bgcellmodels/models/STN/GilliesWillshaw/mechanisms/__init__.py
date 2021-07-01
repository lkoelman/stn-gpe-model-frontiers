"""
NEURON NMODL mechanisms for Gillies & Willshaw (2005) STN cell model.

Mechanisms are loaded into NEURON when module is imported.
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import submodule as submodule_alias
# from .submodule import myclass

import os.path
import neuron

here = os.path.abspath(os.path.dirname(__file__))
neuron.load_mechanisms(here)