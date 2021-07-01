"""
Subthalamic nucleus (STN) cell models.

@author     miscellaneous, see individual models
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import submodule as submodule_alias
# from .submodule import myclass

import os
from neuron import h

pkg_dir = os.path.abspath(os.path.dirname(__file__))

templates = {
    "STN_morph_arcdist": "stn_proto_arcdist.hoc",
    "STN_morph_cartdist": "stn_proto_cartdist.hoc",
    "STN_morph_type1RD": "stn_proto_type1RD.hoc",
}

_loaded_templates = []

def load_template(template_name):
    """
    Load Hoc code with template definition.
    """
    global loaded_templates

    if template_name not in templates:
        raise ValueError(
            "Unknown template: {}. Available templates are: {}".format(
                template_name, ", ".join(templates.keys())))

    # Load Hoc template
    if template_name not in _loaded_templates:
        prev_wd = os.getcwd()
        os.chdir(pkg_dir)
        h.load_file(templates[template_name])
        os.chdir(prev_wd)
        _loaded_templates.append(template_name)


swc2gillies = {
    0: 0,
    1: 1,
    2: 3,
    3: 5,
    4: 6,
    5: 4,
    6: 2,
    7: 7,
    8: 8,
    9: 9,
    10: 10,
    11: 0,
    12: 1,
    13: 3,
    14: 5,
    15: 7,
    16: 8,
    17: 6,
    18: 4,
    19: 9,
    20: 10,
    21: 11,
    22: 12,
    23: 2,
    24: 13,
    25: 15,
    26: 17,
    27: 18,
    28: 16,
    29: 14,
    30: 19,
    31: 20,
    32: 21,
    33: 22,
}


def swc_to_gillies_index(index, secarray_name='soma'):
    """
    Mapping from SWC section indices to compartmental indices in the
    original Gillies & Willshaw model.

    This is useful for loading the saved conductance values for each tree.

    @return     (tree_index, array_index) : tuple
    """
    if secarray_name == 'soma':
        # All somatic compartments map to single soma compartment
        return -1, 0
    elif index <= 10:
        # Indices up to 10 map to 'dend1'
        return 1, swc2gillies[index]
    else:
        # Indices over 10 map to 'dend0' and have nonlinear mapping
        return 0, swc2gillies[index]