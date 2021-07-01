# Enable interactive plots (%matplotlib -l to list backends)
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from bgcellmodels.common import analysis, units, morphology, treeutils
import bgcellmodels.models.STN.GilliesWillshaw.gillies_pynn_model as gillies_model
import neuron; h = neuron.h
import bluepyopt.ephys as ephys

# Load custom synapse mechanisms
neuron.load_mechanisms('/home/luye/workspace/bgcellmodels/mechanisms/synapses')


################################################################################
# Create cell

cell = gillies_model.StnCellModel()
icell = cell.icell
nrnsim = cell.sim

named_seclists =  {listname: list(getattr(icell, listname)) for listname in cell.seclist_names}
for k, v in named_seclists.items():
    if len(v)==0:
        named_seclists.pop(k) # don't include empty SectionLists
    else:
        print("{} : {} sections".format(k, len(v)))

somatic = named_seclists['somatic']
dendritic = named_seclists['basal']

soma = somatic[0]
dend = dendritic[0]

nseg = sum((sec.nseg for sec in icell.all))
print("Total number of compartments: {}".format(nseg))

# Give our cell 3D coordinates and shape info
h.define_shape(sec=soma)

################################################################################
# LFP Electrode
import lfpsim

sigma = extracellular_conductivity = 0.3
coords = electrode_coordinates = h.Vector([50.0, 50.0, 50.0])

#-------------------------------------------------------------------------------
# Method 1 : Summator mechanism
# summator = h.insert_lfp_summator(soma)
# h.add_lfp_sources(summator, "PSA", sigma, coords, icell.somatic, icell.basal)

#-------------------------------------------------------------------------------
# Method 1 : lfp_src mechanism + custom fadvance()

# h.insert_lfp_sources(icell.somatic, icell.basal)
# h.initialize_lfp_factors("PSA", sigma, coords, icell.somatic, icell.basal)

################################################################################
# Record & Run

# Define traces
# rec_secs = {
#     'soma': soma,
#     'soma_spiker': h.NetCon(soma(0.5)._ref_v, None, -10.0, 0, 0),
#     'lfp_sum': summator,
# }

# trace_specs = {
#     't_global': {'var': 't'},
#     'V_soma': {'var':'v', 'sec':'soma', 'loc':0.5},
#     'AP_soma': {'netcon': 'soma_spiker'},
#     # LFP signal
#     'LFP': {'pointp':'lfp_sum', 'var':'summed'},
# }

# rec_dt = 0.05
# vec_dict, markers = analysis.recordTraces(rec_secs, trace_specs, rec_dt)


# Init and run simulation
h.cvode_active(0)
h.dt = 0.025
h.celsius = 35.0
h.set_aCSF(4) # Hoc function defined in Gillies code
h.v_init = -68.0
h.tstop = 2000.0
h.init()

import time
tstart = time.time()
nrnsim.run(h.tstop, h.dt)
tstop = time.time()
cputime = tstop - tstart
num_segments = sum((sec.nseg for sec in h.allsec()))
print("Simulated {} segments for {} ms in {} s CPU time".format(
        num_segments, h.tstop, cputime))