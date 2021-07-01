"""
Test Van Rossum (2000) activity-dependent weight scaler
"""

from neuron import h
import numpy as np

h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library


cell = h.Section()
cell.insert('hh')
cell.insert('pas')

izhi = h.IntFire2()
izhi.ib = 1.05

# Make weight scaler
print("Making scaler...")
hpwa = h.VanRossumHPWA(cell(0.5))
hpwa.scaling = 0
hpwa.set_target_rate(20.0)

# Send cell spikes to weight scaler
print("Sending spikes...")
if izhi is None:
    spike_rec = h.NetCon(cell(0.5)._ref_v, hpwa, -10, 1, 1, sec=cell)
    spike_rec.threshold = -10
else:
    spike_rec = h.NetCon(izhi, hpwa)
spike_rec.delay = 1
spike_rec.weight[0] = 1


spike_vec = h.Vector()
spike_rec.record(spike_vec)

# Start weight scaling after a few seconds
print("Make init handlers...")
# nc = h.NetCon(None, hpwa) # NetCon to turn on/off NetStim
# nc.weight[0] = 2.0 # signal turn on
def enable_scaling():
    hpwa.scaling = 1
    # nc.event(200)
fih = h.FInitializeHandler(enable_scaling)

# Stimulator for cell
stim = h.IClamp(cell(0.5))
stim.dur = 1e9
stim.delay = 5
stim.amp = 1e3

print("Initial weight value is {}".format(stim.amp))

# Add an excitatory weight for scaling
print("Adding weight ref...")
if izhi is None:
    h.setpointer(stim._ref_amp, 'temp_wref', hpwa)
else:
    h.setpointer(izhi._ref_ib, 'temp_wref', hpwa)
hpwa.add_wref(1)

# Record variables
print("Recording variables...")
w_rec = h.Vector()
w_rec.record(stim._ref_amp)
v_rec = h.Vector()
if izhi is None:
    v_rec.record(cell(0.5)._ref_v)
else:
    v_rec.record(izhi._ref_m)

print("Simulating...")
h.dt = 0.025
h.v_init = -68
h.celsius = 35
# h.finitialize()
h.tstop = 10e3
h.run()


print("Final weight value is {}".format(stim.amp))

import matplotlib.pyplot as plt
fig, axes = plt.subplots(2,1)

wvec = w_rec.as_numpy()
tvec = np.arange(wvec.size) * h.dt
vvec = v_rec.as_numpy()

ax = axes[0]
ax.plot(tvec, wvec)
ax.vlines(spike_vec.as_numpy(), wvec.min(), wvec.max(), color='r')

ax = axes[1]
ax.plot(tvec, vvec)

plt.show(block=False)