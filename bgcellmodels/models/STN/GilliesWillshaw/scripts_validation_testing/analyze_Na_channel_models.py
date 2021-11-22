"""
Test gating and dynamics of various Na channel models

@author Lucas Koelman
@date   27-02-2017
"""

from test_stn_Gillies import * # as if we're working in that file

def rec_currents_states(traceSpecs, sec_tag, sec_loc, mechs):
    """ Specify trace recordings for all ionic currents
        and activation variables 

    @param  traceSpecs  collections.OrderedDict of trace specifications

    @param  sec_tag     section tag, i.e. key in section dictionary passed
                        to analysis.recordTraces(). Note that the first two
                        characters must be unique

    @param  sec_loc     location in the section to record the variable
                        (location maps to a segment where the var is recorded)

    @effect             for each ionic current, insert a trace specification
                        for the current, open fractions, activation, 
                        and inactivation variables
    """

    # Derive suffix for traces from section tag
    suf = '_' + sec_tag[0:2]
    ts = traceSpecs

    for mech in mechs:
        # Na currents
        ts['I_'+mech+suf] =     {'sec':sec_tag,'loc':sec_loc,'mech':mech,'var':'ina'}

        # Na Channel open fractions
        ts['O_'+mech+suf] =     {'sec':sec_tag,'loc':sec_loc,'mech':mech,'var':'O'}

        # Na channel closes fractions
        ts['Cn_'+mech+suf] =  {'sec':sec_tag,'loc':sec_loc,'mech':mech,'var':'Ctot'}

        # Na channel inactivated fractions
        ts['In_'+mech+suf] =  {'sec':sec_tag,'loc':sec_loc,'mech':mech,'var':'Itot'}

        # Na channel blocked fractions
        ts['B_'+mech+suf] =     {'sec':sec_tag,'loc':sec_loc,'mech':mech,'var':'B'}
    
def plot_currents_states(recData, recordStep, channames, timeRange=None, sec_tag=None):
    """ Plot currents and (in)activation variable for each ionic
        current in the same axis. Ionic currents are grouped per ion
        in one figure, and the x-axes are synchronized for zooming
        and panning.

    @param recData      traces recorded using a trace specification provided
                        by rec_currents_activations()

    @return             tuple figs, cursors where figs is a list of figures
                        that were created and cursors a list of cursors
    """
    figs = []
    cursors = []

    # Plot activations-currents on same axis per current
    fig, axrows = plt.subplots(len(channames), 1, sharex=True) # one plot for each channel
    for i, mechname in enumerate(channames):
        if sec_tag is None:
            # no sec_tag: find all traces that are marked with channel (e.g. 'CaT')
            pat = re.compile(r'^[A-Za-z]{1,3}_' +mechname) # find 'char+_+abbr' at beginning of word
        else:
            # if sec_tag: find traces marked with channel AND sec_tag
            sec_suffix = sec_tag[0:2]
            pat = re.compile(r'^[A-Za-z]{1,3}_' + mechname + r'\w+' + sec_suffix + r'$')
        chanFilter = lambda x: re.search(pat, x) # variables plotted on left axis match this filter
        twinFilter = lambda x: x.startswith('I_') # vars plotted on right axis match this filter
        # Plot traces that match pattern
        analysis.cumulPlotTraces(recData, recordStep, showFig=False, 
                            fig=None, ax1=axrows[i], yRange=(-0.1,1.1), timeRange=timeRange,
                            includeFilter=chanFilter, twinFilter=twinFilter, solid_only=True)
        # Add figure interaction
        cursor = matplotlib.widgets.MultiCursor(fig.canvas, fig.axes, 
                    color='r', lw=1, horizOn=False, vertOn=True)
    figs.append(fig)
    cursors.append(cursor)

    return figs, cursors

def test_spontaneous():
    """
    Play back spontaneous action potentials through a voltage clamp and record
    how channels react.

    OBSERVATIONS
    - for voltage protocol to trigger slow inactivation, you see the following
      in the transient_resurgent models: during the 5 ms step to +30mV the occupation
      of the B state and Itot states saturates. Then when you drop to -40mV, occupation
      of the B state decreases exponentially while occupation of Itot rises to a new
      saturation level. This is the transition B->O->Ii, i.e. recovery from the blocked
      state through the open state (accompanied by a sizeable resurgent current) into
      the 'slow'/'tight' inactivated states.

    TODO
    - test natural firing under depolarizing stimulation, then bring to -40 to trigger
      transition into tight inactivated state, then stimulate again and see if firing
      is slower (less excitable)
    """
    # Make dummy cell - will be voltage clamped so channels unimportant for dynamics
    all_Ra = 150.224
    all_cm = 1.0
    soma_L = 18.8
    soma_diam = 18.3112

    # Set channels and conductances
    mechs_chans = dict(gillies_gdict)
    test_mechs = ['NaKB', 'NaTB', 'NaRB', 'NaKR']
    mechs_chans['NaKB'] = ['gbar']
    mechs_chans['NaTB'] = ['gbar']
    mechs_chans['NaRB'] = ['gbar']
    mechs_chans['NaKR'] = ['gbar']
    glist = [gname+'_'+mech for mech,chans in mechs_chans.items() for gname in chans]
    gbar_default = {
        'gna_Na':   1.483419823e-02, # global default var
        'gna_NaL':  1.108670852e-05, # global default var
        'gcaL_HVA': 0.0009500, # from file
        'gcaN_HVA': 0.0011539, # from file
        'gcaT_CaT': 0.0000000, # from file
        'gk_Ih':    0.0010059, # from file
        'gk_KDR':   0.0038429, # from file
        'gk_Kv31':  0.0134086, # from file
        'gk_sKCa':  0.0000684, # from file
    }

    # Create soma
    soma = h.Section()
    soma.nseg = 1
    soma.Ra = all_Ra
    soma.diam = soma_diam
    soma.L = soma_L
    soma.cm = all_cm
    for mech in test_mechs:
        soma.insert(mech)
    # for k,v in gbar_default.items():
    #     setattr(soma, k, v)
    # setionstyles_gillies(soma)

    # Insert voltage clamp (space clamp)
    clamp = h.SEClamp(soma(0.5)) # NOTE: use SEClamp instead of IClamp to hyperpolarize to same level independent of passive parameters
    clamp.dur1 = 2000
    ptr = clamp._ref_amp1
    ptr_parent = clamp

    # Play vector into voltage clamp
    replay = False
    if replay:
        # csv_path = "C:\\Users\\lkoelman\\cloudstore_m\\simdata\\fullmodel\\spont_fullmodel_Vm_dt25e-3_0ms_2000ms.csv"
        csv_path = "C:\\Users\\lkoelman\\cloudstore_m\\simdata\\fullmodel\\cstim_fullmodel_Vm_dt25e-3_0ms_2000ms.csv"
        tv = np.loadtxt(csv_path, delimiter=',', usecols=(0, 1), unpack=False)
        t = tv[:,0] * 1e3
        v = tv[:,1] * 1e3
    else:
        # See protocols in Do & Bean (2003) Figs. 4A & 7A
        v = [-90, +30, -40, -90]
        t = [0, 700, 705, 805]
    amp_vec = h.Vector(v)
    t_vec = h.Vector(t)
    amp_vec.play(ptr_parent, ptr, t_vec)

    # Run experiment
    # Set simulation parameters
    dur = 2000
    h.dt = 0.025
    h.celsius = 22 # 22 = room temperature (no Q10 adjustment)
    h.v_init = -60
    set_aCSF(4) # Set initial ion concentrations from Bevan & Wilson (1999)

    # Specify traces to record
    secs = {'soma': soma}
    traceSpecs = collections.OrderedDict() # for ordered plotting (Order from large to small)
    traceSpecs['V_soma'] = {'sec':'soma','loc':0.5,'var':'v'}
    traceSpecs['t_global'] = {'var':'t'}
    rec_currents_states(traceSpecs, 'soma', 0.5, test_mechs)

    # Connect pointers to recording vectors
    recordStep = 0.025
    recData, _ = analysis.recordTraces(secs, traceSpecs, recordStep)

    # Simulate
    h.tstop = dur
    h.init() # calls finitialize() and fcurrent()
    h.run()

    # Plot membrane voltages
    recV = collections.OrderedDict([(k,v) for k,v in recData.items() if k.startswith('V_')])
    figs_vm = analysis.plotTraces(recV, recordStep, yRange=(-80,40), traceSharex=True)
    vm_fig = figs_vm[0]
    vm_ax = figs_vm[0].axes[0]

    # Plot ionic currents, (in)activation variables
    figs, cursors = plot_currents_states(recData, recordStep, test_mechs, timeRange=[0.,1000.])
    plt.show(block=False)

    globals().update(locals())
    return recData, figs, cursors


if __name__ == '__main__':
    test_spontaneous()