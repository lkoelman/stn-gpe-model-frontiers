"""
Python code to run and plot Gillies & Willshaw (2006) model

@author Lucas Koelman
@date   21-10-2016

NOTE: this script relies on commenting out the 'graphics=1' line in sample.hoc
"""

import neuron
h = neuron.h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library


from bgcellmodels.common.nrnutil import ExtSecRef, getsecref

# Load NEURON mechanisms
import os.path
script_dir = os.path.dirname(__file__)
neuron.load_mechanisms(os.path.join(script_dir, 'mechanisms'))
# neuron.load_mechanisms(os.path.join(script_dir, 'mechanisms_extra'))

# Map of channel mechanisms to max conductance parameters (key = suffix of mod mechanism)
gillies_gdict = {
    'STh':  ['gpas'],                               # passive/leak channel
    'Na':   ['gna'], 'NaL':['gna'],                 # Na channels
    'KDR':  ['gk'], 'Kv31':['gk'], 'sKCa':['gk'],   # K channels
    'Ih':   ['gk'],                                 # nonspecific channels
    'CaT':  ['gcaT'], 'HVA':['gcaL', 'gcaN'],       # Ca channels
    'Cacum': [],                                    # No channels
}
gbar_dict = gillies_gdict
gleak_name = 'gpas_STh'

# Mechanism parameters that are changed from default values in original model code
mechs_params_dict = {
    'STh':  ['gpas'],
    'Na':   ['gna'],
    'NaL':  ['gna'],
    'KDR':  ['gk'],
    'Kv31': ['gk'],
    'sKCa': ['gk'],
    'Ih':   ['gk'],
    'CaT':  ['gcaT'],
    'HVA':  ['gcaL', 'gcaN'],
    'Cacum':[],
}

# All mechanism parameters that are not conductances
mechs_params_nogbar = dict(mechs_params_dict)
for mech, params in mechs_params_nogbar.iteritems():
    for gbar_param in gbar_dict.get(mech, []):
        try:
            params.remove(gbar_param)
        except ValueError:
            pass

# All mechanism names
gillies_mechs = list(gillies_gdict.keys()) # all mechanisms
mechs_list = gillies_mechs

# All gbar (max conductance) names
gillies_glist = [
    gname+'_'+mech for mech,chans in gillies_gdict.iteritems() 
                    for gname in chans
]
gbar_list = gillies_glist
active_gbar_names = [gname for gname in gillies_glist if gname != gleak_name]


def stn_cell_gillies():
    """
    Initialize Gillies & Willshaw cell model as singleton
    (only one copy of cell possible)
    """
    if not hasattr(h, 'SThcell'):
        neuron.h.load_file(os.path.join(script_dir, 'gillies_create_singleton.hoc'))
    else:
        print("Gillies STN cell already exists. Cannot create more than one instance.")
    
    soma = h.SThcell[0].soma
    dends = h.SThcell[0].dend0, h.SThcell[0].dend1
    stims = h.stim1, h.stim2, h.stim3

    return soma, dends, stims


def stn_cell_standardized():
    """
    Create STN cell using standardized BluePyOpt-compatible prototype
    that allows multiple copies of the cell to exist, e.g. for network
    simulations.
    """
    if not hasattr(h, 'SThcell'):
        neuron.h.load_file(os.path.join(script_dir, 'gillies_cell_factory.hoc'))
    cell_idx = h.make_stn_cell_global()
    cell_idx = int(cell_idx)
    return h.SThcells[cell_idx]


def get_stn_refs():
    """
    Make SectionRef for each section and assign identifiers
    """
    if not hasattr(h, 'SThcell'):
        stn_cell_gillies()
    
    somaref = ExtSecRef(sec=h.SThcell[0].soma)
    dendLrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend0] # 0 is left tree
    dendRrefs = [ExtSecRef(sec=sec) for sec in h.SThcell[0].dend1] # 1 is right tree
    allsecrefs = [somaref] + dendLrefs + dendRrefs
    
    for noderef in allsecrefs:
        
        # Assign indices in /sth-data/treeX-nom.dat
        if noderef in dendLrefs:
            noderef.tree_index = 0
            noderef.table_index = dendLrefs.index(noderef) + 1
        
        elif noderef in dendRrefs:
            noderef.tree_index = 1
            noderef.table_index = dendRrefs.index(noderef) + 1
        
        elif noderef is somaref:
            noderef.tree_index = -1
            noderef.table_index = 0

        # Assign a unique GID based on table and tree index
        noderef.gid = min(0,noderef.tree_index)*100 + noderef.table_index

    return somaref, dendLrefs, dendRrefs


def get_soma_refs(all_refs):
    """
    Return SectionRef to soma sections
    """
    return [ref for ref in all_refs if ref.sec.same(h.SThcell[0].soma)]


def get_each_dend_refs(all_refs):
    """
    Get one list of SectionRef for each dendrite.
    """
    dend0 = [getsecref(sec, all_refs) for sec in h.SThcell[0].dend0]
    dend1 = [getsecref(sec, all_refs) for sec in h.SThcell[0].dend0]
    return dend0, dend1


def get_all_dend_refs(all_refs):
    """
    Return list of SectionRef to unique dendritic sections.
    """
    dend0, dend1 = get_each_dend_refs(all_refs)
    return dend0 + dend1


def make_passive(sec_refs, save_gbar=True):
    """
    Make given sections passive by setting active conductances to zero.

    @param  sec_refs    list(SectionRef) of sections to make passive

    @param  save_gbar   if True, current conductance values will be savec on the
                        SectionRef object in a dict named 'initial_params'
    """
    from reducemodel import redutils
    for ref in sec_refs:
        # Store original conductances
        redutils.store_seg_props(ref, gillies_gdict, attr_name='initial_params')

        # Set active conductances to zero
        for seg in ref.sec:
            for gbar in active_gbar_names:
                setattr(seg, gbar, 0.0)

    # NOTE: reset as follows:
    # for ref in sec_refs:
    #     # Restore parameter dict stored on SectionRef
    #     redutils.set_range_props(ref, ref.initial_params)


def reset_channel_gbar():
    """
    Reset all channel conductances to initial state.

    NOTE: initialization copied from sample.hoc/gillies_cell_singleton.hoc
    """
    # Soma
    h.SThcells[0].soma.gna_Na = h.default_gNa_soma
    h.SThcells[0].soma.gna_NaL = h.default_gNaL_soma

    # Dendrites
    h.cset(0,"gk_KDR","")
    h.cset(0,"gk_Kv31","")
    h.cset(0,"gk_Ih","")
    h.cset(0,"gk_sKCa","")
    h.cset(0,"gcaT_CaT","")
    h.cset(0,"gcaN_HVA","")
    h.cset(0,"gcaL_HVA","")


def setionstyles_gillies(sec):
    """
    Set ion styles to work correctly with membrane mechanisms
    """
    sec.push()
    h.ion_style("na_ion",1,2,1,0,1)
    h.ion_style("k_ion",1,2,1,0,1)
    h.ion_style("ca_ion",3,2,1,1,1)
    h.pop_section()


def set_aCSF(req):
    """
    Set global initial ion concentrations (artificial CSF properties)

    This is a Python version of the Hoc function set_aCSF()

    @param req      int: identifier

                    0 = NEURON defaults

                    3 = Beurrier et al (1999)

                    4 = Bevan & Wilson (1999)

                    5 = equilibrium concentrations at 35 degrees celsius

    NOTE: only cai is actually changed during the simulation
    """

    if req == 3: # Beurrier et al (1999)
        h.nai0_na_ion = 15
        h.nao0_na_ion = 150

        h.ki0_k_ion = 140
        h.ko0_k_ion = 3.6

        h.cai0_ca_ion = 1e-04
        h.cao0_ca_ion = 2.4

        h('cli0_cl_ion = 4') # self-declared Hoc var
        h('clo0_cl_ion = 135') # self-declared Hoc var

    elif req == 4: # Bevan & Wilson (1999)
        h.nai0_na_ion = 15
        h.nao0_na_ion = 128.5

        h.ki0_k_ion = 140
        h.ko0_k_ion = 2.5

        h.cai0_ca_ion = 1e-04
        h.cao0_ca_ion = 2.0

        h('cli0_cl_ion = 4')
        h('clo0_cl_ion = 132.5')

    elif req == 0: # NEURON's defaults
        h.nai0_na_ion = 10
        h.nao0_na_ion = 140

        h.ki0_k_ion = 54
        h.ko0_k_ion = 2.5

        h.cai0_ca_ion = 5e-05
        h.cao0_ca_ion = 2

        h('cli0_cl_ion = 0')
        h('clo0_cl_ion = 0')

    elif req == 5: # equilibrium concentrations (average) at 35 degree celsius

        h.nai0_na_ion = 15
        h.nao0_na_ion = 128.5

        h.ki0_k_ion = 140
        h.ko0_k_ion = 2.5

        h.cai0_ca_ion = 0.04534908688919702 # min:0.019593383952621085 / max: 0.072908152365581
        h.cao0_ca_ion = 2.0

        h('cli0_cl_ion = 4')
        h('clo0_cl_ion = 132.5')


def applyApamin(soma, dends):
    """
    Apply apamin (reduce sKCa conductance)

    @param soma     soma Section

    @param dends    list of dendritic Sections objects
    
    NOTE: in paper they say reduce by 90 percent but in code
    they set everything to 0 except in soma where they divide
    by factor 10
    """
    soma(0.5).__setattr__('gk_sKCa', 0.0000068)
    for sec in dends:
        for iseg in range(1, sec.nseg+1):
            xnode = (2.*iseg-1.)/(2.*sec.nseg) # arclength of current node (segment midpoint)
            sec(xnode).__setattr__('gk_sKCa', 0.0)


def stn_init_physiology():
    """
    Initialize STN cell in biologically plausible physiological state.
    """
    h.celsius = 35
    h.v_init = -68.0
    h.set_aCSF(4)
    h.init()


################################################################################
# Simulations
################################################################################


def runtest_multithreaded(testfun, nthread):
    """
    Run a test procol in multithreaded mode
    """
    # make cell
    stn_cell_gillies()

    # enable multithreaded execution
    h.cvode_active(0)
    h.load_file('parcom.hoc')
    pct = h.ParallelComputeTool[0]
    pct.nthread(4)
    pct.multisplit(1)
    pct.busywait(1)

    # Do test
    t0 = h.startsw()
    testfun()
    t1 = h.startsw(); h.stopsw()
    print("Elapsed time: {} ms".format(t1-t0))


def runtest_singlethreaded(testfun, use_tables=True):
    """
    Run a test protocol in single-threaded mode
    """
    # make cell
    stn_cell_gillies()

    # Disable tables
    if not use_tables:
        for at in dir(h):
            if at.startswith('usetable_'):
                setattr(h, at, 0)
                print("Disabled TABLE {}".format(at))

    # Do test
    t0 = h.startsw()
    testfun()
    t1 = h.startsw(); h.stopsw()
    print("Elapsed time: {} ms".format(t1-t0))

if __name__ == '__main__':
    stn_cell_gillies()
    # runtest_reboundburst()

    # runtest_singlethreaded(runtest_reboundburst)
    # runtest_multithreaded(runtest_reboundburst, 6)