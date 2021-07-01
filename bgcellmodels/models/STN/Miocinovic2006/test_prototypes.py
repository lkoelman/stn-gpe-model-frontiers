"""
Run Miocinovic STN neuron with extracellular stimulation

@author Lucas Koelman
@date   06-04-2017
@note   must be run from script directory or .hoc files not found

"""

# Make sure other modules are on Python path
import os
script_dir = os.path.dirname(__file__)
main_cwd = os.getcwd()

# Load NEURON
import neuron
from neuron import h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# Load external libs
import bluepyopt.ephys as ephys
import bgcellmodels.common.electrotonic as electrotonic
import bgcellmodels.morphology.morph_ni as morph_ni

# Load NMODL and Hoc code
import bgcellmodels.models.STN.GilliesWillshaw.mechanisms
os.chdir(script_dir)
h.xopen("stn_proto_cartdist.hoc")
h.xopen("stn_proto_arcdist.hoc")
os.chdir(main_cwd)

# Global variables
nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)


def test_stn_prototype(template_name='STN_morph_arcdist',
                       axon_builder=None,
                       streamline_path=None,
                       export_locals=False):
    """
    Test generic STN prototype based on cartesian distance calculation.
    """

    # Instantiate cell template
    template_function = getattr(h, template_name)
    icell = template_function()
    icell.with_extracellular = 0

    # Load morphology into template
    morph_path = os.path.join(script_dir, 
                              'morphologies/type1RD_axonless-with-AIS.swc')
    morphology = ephys.morphologies.NrnFileMorphology(morph_path,
                                                      do_replace_axon=False)
    morphology.instantiate(sim=nrnsim, icell=icell)
    
    # Setup biophysical properties
    icell.del_unused_sections()
    icell.insert_biophys()
    nseg_extra = electrotonic.set_min_nseg_hines(icell.all, f_lambda=100.0)
    print("Created {} extra segments to satisfy Hines' rule".format(nseg_extra))
    icell.set_biophys_spatial()

    # Build the axon
    tracks_coords = morph_ni.load_streamlines(streamline_path, max_num=1, min_length=2.0)
    axon_coords = tracks_coords[0]
    axon_conn_sec = icell.axon_terminal_ref().sec
    axon = axon_builder.build_along_streamline(axon_coords,
                terminate='nodal_cutoff', interp_method='cartesian',
                parent_cell=icell, parent_sec=axon_conn_sec)

    
    # append_axon(self.icell)
    if export_locals:
        globals().update(locals())


if __name__ == '__main__':
    # Test all templates
    # templates = ['STN_morph_cartdist', 'STN_morph_arcdist']
    # for template_name in templates:
    #     test_stn_prototype(export_locals=False)

    streamline_path = '/home/luye/Documents/mri_data/Waxholm_rat_brain_atlas/WHS_DTI_v1_ALS/S56280_track_filter-ROI-STN.tck'

    from bgcellmodels.models.axon.mcintyre2002 import AxonMcintyre2002
    axon_builder = AxonMcintyre2002()

    # Test single template
    test_stn_prototype(template_name='STN_morph_arcdist',
        axon_builder=axon_builder,
        streamline_path=streamline_path,
        export_locals=True)