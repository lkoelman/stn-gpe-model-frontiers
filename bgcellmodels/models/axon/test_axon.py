"""
Tests for axon builder and axon models.
"""

# Load NEURON
from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# External modules
import numpy as np

# Our modules
import bgcellmodels.morphology.morph_ni as morph_ni
from bgcellmodels.models.axon.mcintyre2002 import AxonMcintyre2002

def test_gpe_axon():
    """
    Test case: add axon to GPe neuron.
    """
    raise NotImplementedError("See file 'calibrate_morphologies_GPe.ipnb'")

def test_without_collaterals():
    """
    Test case without axon collaterals.
    """

    streamline_path = '/home/luye/Documents/mri_data/Waxholm_rat_brain_atlas/WHS_DTI_v1_ALS/S56280_track_filter-ROI-STN.tck'

    
    axon_builder = AxonMcintyre2002()

    # Build the axon
    tracks_coords = morph_ni.load_streamlines(streamline_path, max_num=1, min_length=2.0)
    axon_coords = tracks_coords[0]
    axon_conn_sec = h.Section(name='parent')
    axon = axon_builder.build_along_streamline(axon_coords,
                terminate='nodal_cutoff', interp_method='cartesian',
                parent_cell=None, parent_sec=None)


    # Write original and reconstructed axon
    build_results = {
        'original': [np.array(axon_builder.streamline_pts)],
        'reconstructed': [np.array(axon_builder.interp_pts)],
    }
    import pickle
    out_fpath = 'compare_axon_pre-post-build.pkl'
    with open(out_fpath, 'wb') as file:
        pickle.dump(build_results, file)
    print("Wrote streamlines to file: " + out_fpath)

    # Export variables
    globals().update(locals())


if __name__ == '__main__':
    test_without_collaterals()