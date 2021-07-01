"""
Python code to run Khaliq & Raman (2003) Cerebellar Purkinje neuron model

@author Lucas Koelman
@date	25-10-2016

"""

NRN_MECHS_PATH = 'C:\\Users\\lkoelman\\cloudstore_m\\code\\KhaliqRamanBean'

import neuron
from neuron import h

import matplotlib.pyplot as plt
import numpy as np

neuron.load_mechanisms(NRN_MECHS_PATH) # Load custom mechanism binaries (see neuron __init__.py)

def runsim_hoc():
	""" Run original simulation provided with model files """
	h.load_file('resurgesim.hoc')

def runsim_python():
	""" Run same experiments and generate same graphs as in paper 
		but using Python-NEURON instead of Hoc.

		This is an implementation in Python of the scripts in 
		resurgesim.hoc and default.ses

		NOTE: initial concentrations are not set in these model
		files, that's why you see a transient regime in each plot
	"""
	# We don't start neuron via nrngui so load libraries
	h.load_file("stdlib.hoc") # Load the standard library
	h.load_file("stdrun.hoc") # Load the standard run library

	# make section
	soma = h.Section()
	soma.push() # equivalent to `access soma` in Hoc

	# insert membrane mechanisms/channels
	soma.insert('naRsg')
	soma.insert('kpkj')
	soma.insert('kpkj2')
	soma.insert('kpkjslow')
	soma.insert('bkpkj')
	soma.insert('cadiff')
	soma.insert('cap')
	soma.insert('lkpkj')
	soma.insert('hpkj')

	# Section parameters
	soma.L = 20
	soma.diam = 20
	soma.ena = 60
	soma.ek = -88

	### implementation of default.ses ###
	vce = h.VClamp(soma(0.5))
	vchold=-71
	vctestbase=-69

	# parameters for consecutive stimulation levels
	vce.dur[0]=5
	vce.amp[0]=vchold
	vce.dur[1]=100
	vce.amp[1]=vctestbase
	vce.dur[2]=5 
	vce.amp[2]=vchold

	# create stim sequence function
	vcincrement=10
	vcsteps=3

	def voltfamily():
		""" runs three experiments in a row with different holding potential """
		h.tstop = vce.dur[0]+vce.dur[1]+vce.dur[2]

		vce.amp[0] = vchold
		vce.amp[2] = vchold

		for j in range(vcsteps):
			x = vctestbase + j*vcincrement
			vce.amp[1] = x

			h.init()
			h.run()
			if h.stoppedrun():
				break

	def plotvoltfamily(vrec, label):
		""" Plot an experiment """
		plt.figure()
		plt.suptitle(label)
		plt.plot(np.arange(vrec.size())*h.dt, vrec.as_numpy())
		plt.show(block=False)

	## Figure 6B ##
	print("generating 6B (K fast)")

	vchold=-71
	vcincrement=10
	vcsteps=8
	vce.dur[1]=100
	vce.dur[2]=30

	h.steps_per_ms=1
	h.dt=1

	Kfast = h.Vector()
	Kfast.record(soma(0.5)._ref_ik_kpkj)
	voltfamily() # run experiment
	plotvoltfamily(Kfast, "Fig. 6B: K fast (highly TEA-sensitive K current)")

	## Figure 6D ##
	print("generating 6D (K mid)")

	vchold=-71
	vcincrement=10
	vcsteps=9
	vce.dur[1]=100
	vce.dur[2]=30

	h.steps_per_ms=1
	h.dt=1

	Kmid = h.Vector()
	Kmid.record(soma(0.5)._ref_ik_kpkj2)
	voltfamily() # run experiment
	plotvoltfamily(Kmid, "Fig. 6D: K mid (moderately TEA-sensitive K current")

	## Figure 6F ##
	print("generating 6F (K slow)")

	vchold=-71
	vcincrement=10
	vctestbase=-61
	vcsteps=8
	vce.dur[1]=100
	vce.dur[2]=30

	h.steps_per_ms=1
	h.dt=1

	Kslow = h.Vector()
	Kslow.record(soma(0.5).ik_kpkjslow)
	voltfamily() # run experiment
	plotvoltfamily(Kslow, "Fig. 6F: K slow (TEA-insensitive potassium current)")

	## Figure 6H ##
	print("generating 6H (BK)")

	vchold=-90
	vcincrement=10
	vcsteps=5

	vce.dur[0]=2
	vce.dur[1]=20
	vce.dur[2]=5

	h.steps_per_ms=40
	h.dt=.025

	BK = h.Vector()
	BK.record(soma(0.5).ik_bkpkj)
	voltfamily() # run experiment
	plotvoltfamily(BK, "Fig. 6H: BK (Iberiotoxin-sensitive KCa current)")

	## Figure 6J ##
	print("generating 6J (P-type Ca)")

	vchold=-90
	vcincrement=10
	vcsteps=11

	vchold=-90
	vcincrement=10
	vctestbase=-90
	vcsteps=11
	vce.dur[1]=10
	vce.dur[2]=5

	h.steps_per_ms=4
	h.dt=.25

	IP = h.Vector()
	IP.record(soma(0.5).ica_cap)
	voltfamily() # run experiment
	plotvoltfamily(IP, "Fig. 6J: P-type Ca current")

	## Figure 6L ##
	print("generating 6L (Ih)")

	vchold=-90
	vcincrement=10
	vcsteps=11

	vchold=-50
	vcincrement=-10
	vctestbase=-60
	vcsteps=7
	vce.dur[0]=100
	vce.dur[1]=1200
	vce.dur[2]=300

	h.steps_per_ms=1
	h.dt=1

	Ih = h.Vector()
	Ih.record(soma(0.5).i_hpkj)
	voltfamily() # run experiment
	plotvoltfamily(Ih, "Fig 6L: Ih (HCN: hyperpolarization-activated cationic current)")

	## Figure 6N ##
	print("generating 6N (leak)")

	vchold=-90
	vcincrement=10
	vcsteps=11

	vchold=-71
	vcincrement=10
	vctestbase=-91
	vcsteps=5
	vce.dur[0]=10
	vce.dur[1]=100
	vce.dur[2]=30

	h.steps_per_ms=1
	h.dt=1

	leak = h.Vector()
	leak.record(soma(0.5).i_lkpkj)
	voltfamily() # run experiment
	plotvoltfamily(leak, "Fig. 6N Leak currents")

if __name__ == '__main__':
	runsim_python()