"""
Test functionality of STN neuron model (Otsuka (2004)) as it is used
in Kamaravelu2016
"""
import numpy as np
import matplotlib.pyplot as plt

import neuron
from neuron import h
nrn = neuron
hoc = h

from bgcellmodels.common import analysis
import collections

# Load NEURON mechanisms
from sys import platform
import os.path
if platform == "linux" or platform == "linux2":
    raise Exception('Linux paths not set')
    # use neuron.load_mechanisms(<path containining x86_64>)
elif platform == "darwin":
    raise Exception('Mac paths not set')
elif platform == "win32":
	# add this line to nrn/lib/python/neuron/__init__.py/load_mechanisms()
	# from sys import platform as osplatform
	# if osplatform == 'win32':
	# 	lib_path = os.path.join(path, 'nrnmech.dll')
	NRN_MECH_PATH = 'C:\\Users\\lkoelman\\workspace\\bgcellmodels\\Otsuka\\nrn_mechs'
	neuron.load_mechanisms(NRN_MECH_PATH)

# Load NEURON function libraries
hoc.load_file("stdlib.hoc") # Load the standard library
hoc.load_file("stdrun.hoc") # Load the standard run library

# global variables
init_handlers = []

def firingrate(vrec, timerange, dt, thresh=-10):
    """ calculate firing rate
    @param vrec			membrane voltage (np.array)
    @param timerange	e.g. [500, 800]
    """
    istart = int(timerange[0]/dt)
    istop = int(timerange[1]/dt)
    zcross = np.where(np.diff(np.sign(vrec[istart:istop]-thresh)))[0]
    apc = len(zcross)/2
    rate = apc/((timerange[1]-timerange[0])*1e-3)
    return rate

def stn_cell():
	cell = hoc.Section()
	# Section parameters
	cell.diam = 60
	cell.L = 60
	cell.nseg = 1
	cell.Ra = 200
	cell.cm = 1e3
	# Otsuka mechanism parameters
	cell.insert('ionsBG') # prevent NEURON from modifying k & na concentrations
	cell.insert('stn')
	return cell

def init_ions(cell):
	# set initial concentrations on per-cell basis
	cell.push()
	def stn_initions():
		# Initial [Ca]_{i/o} is same in Kumaravelu2016
		cell.cai = 5e-6
		cell.cao = 2
		cell.ki = 105
		cell.ko = 3
		cell.nai = 10
		cell.nao = 108
	return hoc.FInitializeHandler(stn_initions)

def simple_setup():
	""" Make a cell, record its voltage and time axis """
	# make set-up
	cell = stn_cell()

	# set initial concentrations on per-cell basis
	fih = init_ions(cell)
	init_handlers.append(fih)

	# start recording
	vrec = h.Vector()
	vrec.record(cell(0.5)._ref_v)
	trec = h.Vector()
	trec.record(h._ref_t)
	return cell, vrec, trec
	
def test_spontaneous():
	""" Test spontaneous activity of STN cells

	PAPER
	- fig 1.A,B: model neurons fire spontaneuously at ~ 10Hz
	- after AP: memmbrane is repolarized quickly, and then depolarized 
	  to spike threshold following a slow hyperpolarization

	TEST RESULTS

	CURRENTS
	- only INa and ICaT are significant during ISI, ICaT is 1.85 mA/cm2
	  and INa and about 4.2 mA/cm2 while the other currents are < 0.1 mA/cm2
	  	- both INa and ICaT seems to drive spontaneous depolarization

	"""
	
	# Set simulation parameters
	dur = 2000
	h.dt = 0.025
	h.celsius = 30 # Otsuka2004, methods section
	h.v_init = -58 # Otsuka2004, fig 1,2: RMP = -58/60
	
	# make set-up
	cell = stn_cell()
	fih = init_ions(cell)

	# Record
	secs = {'soma': cell}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)
	traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
	traceSpecs['I_Na'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ina'}
	traceSpecs['I_KDR'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikDR'}
	traceSpecs['I_KA'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikA'}
	traceSpecs['I_KCa'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikCA'}
	traceSpecs['I_CaL'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'icaL'}
	traceSpecs['I_CaT'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'icaT'}
	recordStep = 0.05
	recData = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Simulate
	h.tstop = dur
	h.init() # calls finitialize() and fcurrent()
	h.run()

	# Analyze
	analysis.plotTraces(recData, recordStep)
	recI = collections.OrderedDict()
	for key in reversed(recData): recI[key] = recData[key]
	recI.pop('V_soma')
	analysis.cumulPlotTraces(recI, recordStep, cumulate=False)
	# raise Exception('breakpoint')

def test_continuous_depol():
	""" Test STN cell response to continuous depolarizing current

	PAPER
	- fig 1.D,E: response to depolarizing current pulse of 1s/0.5s (10/50/100pA)
	- repetitive firing at < 100Hz _during_ pulse
	- firing frequency increased ~linearly with magnitude of pulse

	IDEAL
	- spike rate adaptation mediated by fast Kv3 or Na slow inactivation

	TEST RESULTS
	- OK: firing frequency increases approx. linearly with stimulation aplitude
	- NOT OK: there is no spike rate adaptation at all
		- e.g. ISI remains 60 ms at beginning and end of firing stimulated by 50 pA

	"""

	# Set simulation parameters
	dur = 2000
	h.dt = 0.025
	h.celsius = 30 # Otsuka2004, methods section
	h.v_init = -58 # Otsuka2004, fig 1,2: RMP = -58/60
	
	# make set-up
	cell = stn_cell()
	fih = init_ions(cell)

	# Insert electrode
	electrode = h.IClamp(cell(0.5))
	cur_dens = 50 # (10/50/100) desired current density (mA/cm2)
	cellarea = np.pi*cell.diam*cell.L # (micron^2)
	electrode.amp = cur_dens * cellarea * 1e-2 # (nA) (1e-2 = product converson factors)
	electrode.delay = 500
	electrode.dur = 1000

	# Record
	secs = {'soma': cell}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)
	traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
	traceSpecs['I_Na'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ina'}
	traceSpecs['I_KDR'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikDR'}
	traceSpecs['I_KA'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikA'}
	traceSpecs['I_KCa'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikCA'}
	traceSpecs['I_CaL'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'icaL'}
	traceSpecs['I_CaT'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'icaT'}
	recordStep = 0.05
	recData = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Simulate
	h.tstop = dur
	h.init() # calls finitialize() and fcurrent()
	h.run()

	# Analyze
	analysis.plotTraces(recData, recordStep)
	# recI = collections.OrderedDict()
	# for key in reversed(recData): recI[key] = recData[key]
	# recI.pop('V_soma')
	# analysis.cumulPlotTraces(recI, recordStep, cumulate=False)
	raise Exception('breakpoint')
	
def test_plateau():
	""" Test plateau potential evoked by applying depolarizing pulse 
		at hyperpolarized level of membrane potential

	PAPER
	- fig 2.A: at hyperpolarized potential, injection of depolarizing current 
	  induces polteau potential that outlasts current injection
	- hold Vm at -75 mV (using Iclamp or Vclamp), then apply short 
	  depolarizing pulse (50 pA or 5 mA/cm2 for 80 ms)
	- firing frequence should decrease towards end of plateau

	TEST RESULTS
	- OK: plateau outlasts current injection
	- OK: firing frequency decreases during plateau

	CURRENTS
	- KA peak decreases during burst (height of peaks during AP), as KA decreases the firing frequency also decreases
		- hence KA seems to responsible for rapid repolarization and maintenance of high-frequency firing
	- KCa peak increases over about 25 ms (height of peaks during AP), and decreases during last 100 ms of burst
	- CaT is the first depolarizing current that rises after realease from hyperpolarization and seems to be
	  responsible for initiation of the rebound burst
	  	- CaT bootstraps burst (bootstraps pos feedback of CaL entry)
	  	- it runs out of fuel during burst and thus may contribute to ending the burst
	  		- this is contradicted by burst at regular Vm: there drop in ICaL clearly ends burst
	- CaL reaces steady maximum peak after approx. 70 ms into the burst, after CaT is already past its peak
		- hypothesis that it seems responsible for prolonging th burst seems plausible
		- burst seems to go on as long as CaT+CaL remains approx. constant, and burst ends as long as CaT too low


	"""
	# Set simulation parameters
	dur = 2000
	h.dt = 0.025
	h.celsius = 30 # Otsuka2004, methods section
	h.v_init = -58 # Otsuka2004, fig 1,2: RMP = -58/60
	
	# make set-up
	cell = stn_cell()
	fih = init_ions(cell)

	# Insert Voltage clamp
	# vce = h.VClamp(cell(0.5))
	# vce.amp[0] = -75
	# vce.dur[0] = 500

	# Insert Current Clamp
	electrode = h.IClamp(cell(0.5))
	electrode.dur = dur
	# electrode.delay = 0
	# electrode.amp = 0
	cellarea = np.pi*cell.diam*cell.L # (micron^2)
	I_hyper = -5.7 * cellarea * 1e-2 # hyperpolarizes to ~-75 mV
	I_depol = 5.0 * cellarea * 1e-2 # 5 mA/cm2 equiv to  50pA
	i_ampvec = h.Vector([I_hyper, I_depol, I_hyper])
	i_tvec = h.Vector([0, 500, 580])
	i_ampvec.play(electrode, electrode._ref_amp, i_tvec)

	# Record
	secs = {'soma': cell}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)
	traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
	traceSpecs['I_Na'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ina'}
	traceSpecs['I_KDR'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikDR'}
	traceSpecs['I_KA'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikA'}
	traceSpecs['I_KCa'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikCA'}
	traceSpecs['I_CaL'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'icaL'}
	traceSpecs['I_CaT'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'icaT'}
	recordStep = 0.05
	recData = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Record electrode
	irec = h.Vector()
	irec.record(electrode._ref_amp)

	# Simulate
	h.tstop = dur
	h.init() # calls finitialize() and fcurrent()
	h.run()

	### Analyze ###

	# Electrode
	# plt.figure()
	# plt.plot(np.arange(irec.size())*h.dt, irec.as_numpy())
	# plt.ylabel('I_stim')

	# Cell voltage/currents
	analysis.plotTraces(recData, recordStep)
	recI = collections.OrderedDict()
	for key in reversed(recData): recI[key] = recData[key]
	recI.pop('V_soma')
	analysis.cumulPlotTraces(recI, recordStep, cumulate=False)

def test_rebound():
	""" Test rebound response to different duration hyperpolarizing pulses

	PAPER
	- fig 4.A: when the amplitude of the hyperpolarizing current pulse was kept constant,
	  the generation of a plateau potetnial was directly related to the duration of the pulse
	  	- stimulate with 50 pA after stim with -1.6e-4 mA/cm2 to prevent spont. activity
   		- [ ] for short duration current pulse (<= 100 ms), you get a small amplitude, short
   		  duration response
   		- [ ] for longer duration pulse (100-200 ms < T), you get long duration plateau potentials
   		- [ ] the duration of the plateau becaume longer with further increase of duration T >> 100-200 ms

   	CURRENTS
   	- CaT is the first depolarizing current that rises after realease from hyperpolarization and seems to be
	  responsible for initiation of the rebound burst
	  	- CaT bootstraps burst (bootstraps positive feedback of CaL entry)
	  	- it runs out of fuel during burst and thus may contribute to ending the burst
	  		- this is contradicted by burst at regular Vm: there drop in ICaL clearly ends burst
	- CaL reaces steady maximum peak after approx. 70 ms into the burst, after CaT is already past its peak
		- burst ends when CaL rapidly drops in peak amplitude
		- hypothesis that it seems responsible for prolonging th burst seems plausible
	"""

	# Set simulation parameters
	dur = 2000
	h.dt = 0.025
	h.celsius = 30 # Otsuka2004, methods section
	h.v_init = -58 # Otsuka2004, fig 1,2: RMP = -58/60
	
	# make set-up
	cell = stn_cell()
	fih = init_ions(cell)

	# Insert Voltage clamp
	# vce = h.VClamp(cell(0.5))
	# vce.amp[0] = -75
	# vce.dur[0] = 500

	# Insert Current Clamp
	electrode = h.IClamp(cell(0.5))
	electrode.dur = dur
	# electrode.delay = 0
	# electrode.amp = 0
	cellarea = np.pi*cell.diam*cell.L # (micron^2)
	I_hyper = -1.5 * cellarea * 1e-2 # prevents spontaneous firing
	I_depol = 5.0 * cellarea * 1e-2 # 50 pA equivalent to 5 mA/cm2
	T_pulse = 300 # test 100<200<300<400
	i_ampvec = h.Vector([I_hyper, I_depol, I_hyper])
	i_tvec = h.Vector([0, 500, 500+T_pulse])
	i_ampvec.play(electrode, electrode._ref_amp, i_tvec)

	# Record
	secs = {'soma': cell}
	traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)
	traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
	traceSpecs['I_Na'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ina'}
	traceSpecs['I_KDR'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikDR'}
	traceSpecs['I_KA'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikA'}
	traceSpecs['I_KCa'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'ikCA'}
	traceSpecs['I_CaL'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'icaL'}
	traceSpecs['I_CaT'] = {'sec':'soma','loc':0.5,'mech':'stn','var':'icaT'}
	recordStep = 0.05
	recData = analysis.recordTraces(secs, traceSpecs, recordStep)

	# Record electrode
	irec = h.Vector()
	irec.record(electrode._ref_amp)

	# Simulate
	h.tstop = dur
	h.init() # calls finitialize() and fcurrent()
	h.run()

	### Analyze ###

	# Electrode
	# plt.figure()
	# plt.plot(np.arange(irec.size())*h.dt, irec.as_numpy())
	# plt.ylabel('I_stim')

	# Cell voltage/currents
	analysis.plotTraces(recData, recordStep)
	recI = collections.OrderedDict()
	for key in reversed(recData): recI[key] = recData[key]
	recI.pop('V_soma')
	analysis.cumulPlotTraces(recI, recordStep, cumulate=False)


if __name__ == '__main__':
	test_plateau()
	# test_continuous_depol()
	# test_spontaneous()
	# test_deppulse_hypervm()