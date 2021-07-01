"""
Single compartmental model of a Globus pallidus neuron

Tomohiro Fujita, Tomoki Fukai and Katsunori Kitano (2011)
"Influences of membrane properties on phase response curve and
synchronization stability in a model globus pallidus neuron"
Journal of Computational Neuroscience
DOI: 10.1007/s10827-011-0368-2

The original model is proposed by Gunay et al.(2008) as a
multi-compartmental model.  We modified values of ionic conductances
so that the single compartmental model shows the similar firing
activity to that of the original model.
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import thispkg_submodule
# from .mymodule import myclass
# from . import gillies_pynn_model as pyNN_types