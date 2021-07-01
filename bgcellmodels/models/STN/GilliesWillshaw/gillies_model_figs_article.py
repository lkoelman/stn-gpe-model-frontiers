"""
Test suite for Gillies model that reproduces experiments and figures
from original article.
"""

import matplotlib.pyplot as plt
import gillies_model
from neuron import h

def runsim_paper(plotting='NEURON'):
    """
    Run original simulation provided with model files
    """
    if plotting=='NEURON':
        h('graphics = 1') # turn on plotting using NEURON
    else:
        h('graphics = 0') # turn off plotting using NEURON
    
    # Execute hoc file containing simulation
    h.xopen('sample.hoc')

    if plotting=='mpl' or plotting=='matplotlib':

        # Regular firing plots
        fig = plt.figure()
        plt.suptitle("Regular firing")

        plt.subplot(3,1,1)
        plt.plot(h.recapt.as_numpy(), h.recapv.as_numpy())
        plt.xlabel("Action Potential (30 degC)")

        plt.subplot(3,1,2)
        plt.plot(h.recsp1t.as_numpy(), h.recsp1v.as_numpy())
        plt.xlabel("Rest firing at 25 degC")

        plt.subplot(3,1,3)
        plt.plot(h.recsp2t.as_numpy(), h.recsp2v.as_numpy())
        plt.xlabel("Rest firing at 37 degC")

        # Burst firing plots
        plt.figure()
        plt.suptitle("Bursting")

        plt.subplot(4,1,1)
        plt.plot(h.recrbt.as_numpy(), h.recrbv.as_numpy())
        plt.xlabel("Rebound burst (at 35 degC)")

        plt.subplot(4,1,2)
        plt.plot(h.recsrt.as_numpy(), h.recsrv.as_numpy())
        plt.xlabel("Slow rhythmic bursting (Apamin, 37 degC)")

        plt.subplot(4,1,3)
        plt.plot(h.recfrt.as_numpy(), h.recfrv.as_numpy())
        plt.xlabel("Fast rhythmic bursting (Apamin, 37 degC, CaL-10%)")

        plt.subplot(4,1,4)
        plt.plot(h.recmrt.as_numpy(), h.recmrv.as_numpy())
        plt.xlabel("Mixed rhythmic bursting (Apamin, 37 degC, CaL+10%)")

        # if htmlplotting:
        #   plugins.clear(fig)
        #   plugins.connect(fig, plugins.Reset(), plugins.BoxZoom(), plugins.MousePosition())
        #   mpld3.show()
        plt.show()


def runtest_actionpotential():
    """
    Run AP experiment from original Hoc file
    """
    print("*** Action potential form\n")\

    # Set up recording
    recapt = h.Vector()
    recapv = h.Vector()
    recapt.record(h._ref_t)
    recapv.record(h.SThcell[0].soma(0.5)._ref_v)

    # Simulate
    h.celsius = 30
    h.set_aCSF(3) # this sets initial concentrations via global vars
    h.tstop = 500
    h.dt = 0.025
    h.init()
    h.run()

    fig = plt.figure()
    plt.plot(recapt.as_numpy(), recapv.as_numpy())
    plt.xlabel("Action Potential (30 degC)")
    plt.show(block=False)


def runtest_restfiring():
    """
    Run rest firing experiment from original Hoc file
    """
    print("*** Resting firing rate (at 25 & 37 degC) \n")\

    # Set up recording
    recsp1t = h.Vector()
    recsp1v = h.Vector()
    recsp1t.record(h._ref_t)
    recsp1v.record(h.SThcell[0].soma(0.5)._ref_v)

    # Simulate
    h.celsius = 25
    h.set_aCSF(4)
    h.tstop = 2100
    h.dt=0.025
    h.init()
    h.run()

    # Set up recording
    recsp2t = h.Vector()
    recsp2v = h.Vector()
    recsp2t.record(h._ref_t)
    recsp2v.record(h.SThcell[0].soma(0.5)._ref_v)

    # Simulate
    h.celsius = 37
    h.tstop = 2100
    h.dt=0.025
    h.init()
    h.run()

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(recsp1t.as_numpy(), recsp1v.as_numpy())
    plt.xlabel("Rest firing at 25 degC")
    plt.subplot(2,1,2)
    plt.plot(recsp2t.as_numpy(), recsp2v.as_numpy())
    plt.xlabel("Rest firing at 37 degC")
    plt.show(block=False)


def runtest_reboundburst():
    """
    Run rebound burst experiment from original Hoc file
    """
    print("*** Rebound burst (at 35 degC) \n")

    # Set up recording
    recrbt = h.Vector()
    recrbv = h.Vector()
    recrbt.record(h._ref_t)
    recrbv.record(h.SThcell[0].soma(0.5)._ref_v)

    recicat = h.Vector()
    recicat.record(h.SThcell[0].soma(0.5).CaT._ref_iCaT)
    reccai = h.Vector()
    reccai.record(h.SThcell[0].soma(0.5)._ref_cai)

    # Simulate
    h.celsius = 35

    h.stim1.delay = 0
    h.stim1.dur = 1000
    h.stim1.amp = 0.0

    h.stim2.delay = 1000
    h.stim2.dur = 500
    h.stim2.amp = -0.25

    h.stim3.delay = 1500
    h.stim3.dur = 1000
    h.stim3.amp = 0.0

    h.set_aCSF(4)
    h.tstop = 2500
    h.dt=0.025
    h.init()
    h.run()

    # Plot
    plt.figure()
    plt.plot(recrbt.as_numpy(), recrbv.as_numpy())
    plt.xlabel("Rebound (35 degC)")

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(recrbt.as_numpy(), recicat.as_numpy(), label='ICaT')
    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(recrbt.as_numpy(), reccai.as_numpy(), label="[Ca]_i")
    plt.legend()

    plt.show(block=False)


def runtest_slowbursting():
    """
    Run slow rhytmic bursting experiment from original Hoc file
    """
    print("*** Slow rhythmic bursting (at 37 degC) \n")

    # Set up recording
    recsrt = h.Vector()
    recsrv = h.Vector()
    recsrt.record(h._ref_t)
    recsrv.record(h.SThcell[0].soma(0.5)._ref_v)

    # Simulate
    h.celsius = 37

    h.stim1.delay = 0
    h.stim1.dur = 40000
    h.stim1.amp = -0.25

    h.stim2.delay = 0
    h.stim2.dur = 0
    h.stim2.amp = 0.0

    h.stim3.delay = 0
    h.stim3.dur = 0
    h.stim3.amp = 0.0 

    h.set_aCSF(4)
    h.tstop = 8000
    h.applyApamin()
    h.dt=0.025
    h.init()
    h.run()
    h.washApamin()

    plt.figure()
    plt.plot(recsrt.as_numpy(), recsrv.as_numpy())
    plt.xlabel("Slow rhythmic bursting (Apamin, 37 degC)")
    plt.show(block=False)


def runtest_fastbursting():
    """
    Run fast rhythmic bursting experiment from original Hoc file
    """
    print("*** Fast rhythmic bursting (at 37 degC) \n")

    # Set up recording
    recfrt = h.Vector()
    recfrv = h.Vector()
    recfrt.record(h._ref_t)
    recfrv.record(h.SThcell[0].soma(0.5)._ref_v)

    # Simulate
    h.celsius = 37

    h.stim1.delay = 0
    h.stim1.dur = 40000
    h.stim1.amp = -0.35

    h.set_aCSF(4)
    h.tstop = 4000
    h.cset(0,"gcaL_HVA","-dl0.9") # 10% decrease in dendritic linear CaL (see Figure 8A)
    h.applyApamin()
    h.dt=0.025
    h.init()
    h.run()
    h.washApamin()

    plt.figure()
    plt.plot(recfrt.as_numpy(), recfrv.as_numpy())
    plt.xlabel("Fast rhythmic bursting (Apamin, 37 degC, CaL-10%)")
    plt.show(block=False)


def runtest_mixedbursting():
    """
    Run mixed bursting experiment from original Hoc file
    """
    print("*** Mixed rhythmic bursting (at 37 degC) \n")

    # Set up recording
    recmrt = h.Vector()
    recmrv = h.Vector()
    recmrt.record(h._ref_t)
    recmrv.record(h.SThcell[0].soma(0.5)._ref_v)

    # Simulate
    h.celsius = 37

    h.stim1.delay = 0
    h.stim1.dur = 40000
    h.stim1.amp = -0.32

    h.set_aCSF(4)
    h.tstop = 8000
    h.cset(0,"gcaL_HVA","-dl1.1") # 10% increase in dendritic linear CaL (see Figure 8A,B)
    h.applyApamin()
    h.dt=0.025
    h.init()
    h.run()
    h.washApamin()

    plt.figure()
    plt.plot(recmrt.as_numpy(), recmrv.as_numpy())
    plt.xlabel("Mixed rhythmic bursting (Apamin, 37 degC, CaL+10%)")
    plt.show(block=False)


def runalltests():
    """
    Run all experiments from Hoc file
    """
    gillies_model.stn_cell_gillies()
    runtest_actionpotential()
    runtest_restfiring()
    runtest_reboundburst()
    runtest_slowbursting()
    runtest_fastbursting()
    runtest_mixedbursting()


if __name__ == '__main__':
    runalltests()