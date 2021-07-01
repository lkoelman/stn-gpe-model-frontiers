# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from mymodule import myclass

### LOAD HOC FILES IN PACKAGE ###

# import neuron
# import os

# here = os.path.abspath(os.path.dirname(__file__))

# # List Hoc files that should be loaded
# hoc_files = [
#     f for f in os.listdir(here) if f.endswith('.hoc')
# ]

# # Load Hoc files
# for hoc_file in hoc_files:
#     neuron.h.load_file(os.path.join(here, hoc_file))