"""
Model of the rat subthalamic projection neuron published and coded by 
Gillies, A. and Willshaw, D. (2006).

This package contains the original model code as well as adapted versions
of the cell model that are more suitable for network simulations.

@see        Gillies, A. and Willshaw, D. (2006)
            "Membrane channel interactions underlying rat subthalamic projection
            neuron rhythmic and bursting activity"
            J. Neurophys. 95, 2352-2365


@author     Lucas Koelman

@date       20/03/2017
"""
# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from . import thispkg_submodule
# from .mymodule import myclass
# from . import gillies_pynn_model as pyNN_types

"""
Subthalamic nucleus (STN) cell models.

@author     miscellaneous, see individual models
"""

import os
from neuron import h

pkg_dir = os.path.abspath(os.path.dirname(__file__))

templates = {
    "STN_gillies_singleton": "gillies_create_singleton.hoc",
    "STN_gillies_network": "gillies_create_factory.hoc",
}


def load_hoc(hoc_filename):
    """
    Load hoc code in package directory.
    """
    prev_wd = os.getcwd()
    os.chdir(pkg_dir)
    h.xopen(hoc_filename)
    os.chdir(prev_wd)


def load_mechanisms():
    """
    Load NMODL mechanisms.
    """
    from . import mechanisms


def load_template(template_name):
    """
    Load Hoc code with template definition.
    """
    if template_name not in templates:
        raise ValueError(
            "Unknown template: {}. Available templates are: {}".format(
                template_name, ", ".join(templates.keys())))

    load_hoc(templates[template_name])


def load_gbar_dist(gbar_name):
    """
    Return conductance values saved in 'sth-data' files for specific
    ion channel.

    @return gbar : np.array[N, 4]
            Array where each row is (tree_index, array_index, segment_x, gbar)
            with tree_index -1 for soma, 0 for 'dend0', 1 for 'dend1'. 
    """
    import numpy as np
    gbar_filepath = os.path.join(pkg_dir, 'sth-data', 'cell-{}'.format(gbar_name))
    gbar_mat = np.loadtxt(gbar_filepath)
    return gbar_mat