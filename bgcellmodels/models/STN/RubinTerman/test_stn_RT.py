"""
Test STN cell model published in Rubin, Terman (2002).

Reproduce results in paper fig. 1 and described in first paragraph of
'Results' section.

PUBLISHED BEHAVIOUR:
    - [x] fires spontaneously at ~3 Hz
    - [x] f-I relation of graph 1.c
    - [0] AHP duration-I relation of graph 1.d
    - [x] burst duration in function of I=25 pA/um2 for 300/450/600 ms
    - [x] burst duration in function of I=20/30/40 pA/um2 for 300 ms

"""
import numpy as np
import matplotlib.pyplot as plt

import neuron
from neuron import h
nrn = neuron
hoc = h

# Load NEURON mechanisms
NRN_MECH_PATH = '/home/lkoelman/Workspace/bgmodel/nrn_mechs'
neuron.load_mechanisms(NRN_MECH_PATH)

# Load NEURON function libraries
hoc.load_file("stdlib.hoc") # Load the standard library
hoc.load_file("stdrun.hoc") # Load the standard run library

# Which cell parameters to use
SETPAR_RT2002 = True # use parameters published in 2002 paper
SETPAR_RT2004 = not SETPAR_RT2002 # use parameters from .ode file of 2004 paper

# global variables
init_handlers = []

def reset_sim(v_init=-58, temp=33):
    """ Re-initialize everything to run a new simulation """
    # Initialize global variables
    h.dt = 0.025
    h.v_init = v_init
    h.celsius = temp

    # reset simulator
    h.finitialize()
    h.fcurrent()

def calc_rate(vrec, dur, thresh=-10):
    """ calculate firing rate
    @param vrec membrane voltage vector
    @param dur  duration in ms
    """
    zcross = np.where(np.diff(np.sign(vrec-thresh)))[0]
    apc = len(zcross)/2
    rate = apc/(dur*1e-3)
    return rate

def make_setup():
    """ Make a cell, record its voltage and time axis """
    # make set-up
    cell = hoc.Section()

    # Section parameters
    cell.nseg = 1
    cell.diam = 60
    cell.L = 60
    cell.Ra = 200
    cell.cm = 1e3

    # STN mechanism parameters
    cell.insert('stnRT')
    # see parameters tagged 'RT2002' in .mod file
    if SETPAR_RT2002:
        cell.theta_b_stnRT = 0.4
        cell.tau_r0_stnRT = 40.0
        cell.phi_r_stnRT = 0.2
        cell.epsilon_stnRT = 3.75e-5
        cell.bias_amp_stnRT = 0
        # fix for Ca units
        # cell.kca_stnRT = 22.5e1
        # cell.k1_stnRT = 15.0e-1
    elif SETPAR_RT2004:
        cell.theta_b_stnRT = 0.25
        cell.tau_r0_stnRT = 7.1
        cell.phi_r_stnRT = 0.5
        cell.epsilon_stnRT = 5e-05
        cell.bias_amp_stnRT = 25
    else:
        ValueError('Must choose parameters of 2002 or 2004 paper')

    # Initial ion concentration
    cell.push()
    def stn_initions():
        cell.cai = 0.44 # from RT2005 .ode file
    init_handlers.append(hoc.FInitializeHandler(stn_initions))

    # start recording
    vrec = h.Vector()
    vrec.record(cell(0.5)._ref_v)
    trec = h.Vector()
    trec.record(h._ref_t)
    return cell, vrec, trec

def recordvars(cell, electrode):
    """ record al variables of interest """
    # Electrode current
    istim = h.Vector()
    istim.record(electrode._ref_i)

    # I_T related variables
    ainf = h.Vector()
    ainf.record(cell(0.5)._ref_a_inf_stnRT)
    binf = h.Vector()
    binf.record(cell(0.5)._ref_b_inf_stnRT)
    it = h.Vector()
    it.record(cell(0.5)._ref_icaT_stnRT)
    # I_AHP related variables
    cai = h.Vector()
    cai.record(cell(0.5)._ref_cai)
    iahp = h.Vector()
    iahp.record(cell(0.5)._ref_ikAHP_stnRT)
    # I_CaS related variables
    sinf = h.Vector()
    sinf.record(cell(0.5)._ref_s_inf_stnRT)
    icas = h.Vector()
    icas.record(cell(0.5)._ref_icaS_stnRT)

    return {'istim':istim, 'ainf':ainf, 'binf':binf, 'it':it, 
            'cai':cai, 'iahp':iahp, 'sinf':sinf, 'icas':icas}

def test_spont():
    """ Test GPe neuron spontaneous activity
    
    REPRODUCED BEHAVIOUR:
    - [x] fires spontaneously at ~3 Hz
    """
    dur = 5000
    cell, vrec, trec = make_setup()
    electrode = h.IClamp(cell(0.5))
    recvars = recordvars(cell, electrode)

    # stimulation
    idens = 0 # current density
    cellarea = np.pi*cell.diam*cell.L # (micron^2)
    electrode.amp = idens * cellarea * 1e-2 # (nA) (1e-2 = product converson factors)
    electrode.delay = 0
    electrode.dur = dur
    # bias current Iapp
    # cell.bias_amp_stnRT = 0.0

    # simulate for one second
    reset_sim()
    h.tstop = dur
    h.run()

    # Convert results
    tvec = trec.as_numpy()
    sigs = {k:v.as_numpy() for k,v in recvars.items()}

    # plot results
    interval = (0,dur)
    plt.figure()
    # Membrane voltage
    plt.subplot(4,1,1)
    plt.plot(tvec, vrec.as_numpy())
    plt.xlim(interval)
    plt.ylabel('$V_m$ (mV)')
    # I_T gating
    plt.subplot(4,1,2)
    plt.plot(tvec, sigs['ainf']**3*sigs['binf']**2)
    plt.xlim(interval)
    plt.ylabel('$I_T$ gating')
    # I_AHP gating
    plt.subplot(4,1,3)
    plt.plot(tvec, sigs['cai']/(sigs['cai']+cell.k1_stnRT))
    plt.xlim(interval)
    plt.ylabel('$I_{AHP}$ gating')
    # Ca level
    plt.subplot(4,1,4)
    plt.plot(tvec, sigs['cai'])
    plt.xlim(interval)
    plt.ylabel('$[Ca]_i$')
    plt.show()

    # report firing rate
    rate = calc_rate(vrec.as_numpy(), dur)
    print("GP cell firing rate is {0} Hz".format(rate))

def test_f_I_rel():
    """ Test current-frequency relationship of paper fig. 1c

    REPRODUCED BEHAVIOUR:
    - I=20 -> f=
    - I=60 -> f=
    - I=100 -> f=
    - I=160 -> f=

    """
    dur = 5000
    cell, vrec, trec = make_setup()
    electrode = h.IClamp(cell(0.5))

    # stimulation
    allcurrents = [20., 60., 100., 160.]
    idens = 20 # current density, between 0-180 in fig. 1c
    cellarea = np.pi*cell.diam*cell.L # (micron^2)
    electrode.delay = 0
    electrode.dur = dur
    # bias current Iapp
    # cell.bias_amp_stnRT = 0.0

    plt.figure()
    for i, Iapp in enumerate(allcurrents):
        # set stimulation
        electrode.amp = Iapp * cellarea * 1e-2 # (nA)

        # simulate
        reset_sim()
        neuron.run(dur)

        # Convert results
        tvec = trec.as_numpy()
        rate = calc_rate(vrec.as_numpy(), dur)

        # plot results
        interval = (1000,2000)
        plt.subplot(len(allcurrents),1,i+1)
        # Membrane voltage
        plt.plot(tvec, vrec.as_numpy(), label='I={0}/f={1}'.format(Iapp,rate))
        plt.xlim(interval)
        plt.ylabel('$V_m$ (mV)')
        plt.legend()

    plt.show()

def test_ahp_dur():
    """ Test AHP duration
    
    For fig. 1d: set allcurrents=[-20,-60,-80,-120]; electrode.dur=500
    For fig. 1e: set allcurrents=[-25.]; electrode.dur=[300,450,600]
    For fig. 1f: set alcurrents=[-20.,-30.,-40.]; electrode.dur=300
    

    REPRODUCED BEHAVIOUR:
    - I=20 -> f=
    - I=60 -> f=
    - I=100 -> f=
    - I=160 -> f=

    """
    dur = 2500
    cell, vrec, trec = make_setup()
    electrode = h.IClamp(cell(0.5))
    recvars = recordvars(cell, electrode)

    # stimulation
    # fig. 1d
    allcurrents=[-20,-60,-80,-120]
    alldur = [500]
    # fig. 1f
    # allcurrents = [-20., -30., -40.] # fig. 1f
    # alldur = [300.] # fig. 1f
    # fig. 1e
    # allcurrents = [-25.] # fig. 1e
    # alldur = [300., 450., 600.] # fig. 1e
    # electrode properties
    cellarea = np.pi*cell.diam*cell.L # (micron^2)
    electrode.delay = 500
    electrode.dur = 300
    # bias current Iapp
    # cell.bias_amp_stnRT = 0.0

    plt.figure(1) # for membrane voltage
    plt.figure(2) # for I_AHP
    for i, Iapp in enumerate(allcurrents): # fig. 1f
    # for i, edur in enumerate(alldur): # fig. 1e
        # set stimulation
        # Iapp = 25. # fig. 1e
        electrode.amp = Iapp * cellarea * 1e-2 # (nA)
        electrode.dur = alldur[0]

        # simulate
        reset_sim()
        neuron.run(dur)

        # Convert results
        tvec = trec.as_numpy()
        rate = calc_rate(vrec.as_numpy(), dur)

        # plot results
        interval = (0,dur)
        
        # Membrane voltage
        plt.figure(1)
        # plt.subplot(len(alldur),1,i+1) # fig. 1e
        plt.subplot(len(allcurrents),1,i+1) # fig. 1f
        plt.plot(tvec, vrec.as_numpy(), label='I={0}'.format(Iapp))
        plt.xlim(interval)
        plt.ylabel('$V_m$ (mV)')
        plt.legend()

        # I_AHP
        plt.figure(2)
        # plt.subplot(len(alldur),1,i+1) # fig. 1e
        plt.subplot(len(allcurrents),1,i+1) # fig. 1f
        plt.plot(tvec, recvars['iahp'].as_numpy(), label='I={0}'.format(Iapp))
        plt.xlim(interval)
        plt.ylabel('$I_{AHP}$')
        plt.legend()

    plt.show()

if __name__ == '__main__':
    test_ahp_dur()
