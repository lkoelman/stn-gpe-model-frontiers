"""
Setup code for Balbi et al (2015) motoneuron compartmental cell model for use
with NEURON + Python.

@author Lucas Koelman

@date   15-12-2017
"""

import neuron
h = neuron.h
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

# Load NEURON mechanisms
# import os.path
# scriptdir, scriptfile = os.path.split(__file__)
# NRN_MECH_PATH = os.path.normpath(os.path.join(scriptdir, 'channels'))
# neuron.load_mechanisms(NRN_MECH_PATH)

model_secarray_vars = [ # Sections created using hoc 'create secname[N]'
    'soma', 'dend', # somatic and dendritic sections
    'node', 'MYSA', 'FLUT', 'STIN' # axonal sections
]

model_sec_vars = [ # Sections created using hoc 'create secname'
    'AH', 'IS'
]

# Channel mechanisms (key = suffix of mod mechanism) : max conductance parameters
gbar_dict = {
    'gh':       ['ghbar'],      # H channel (Na + K ions)
    'kca2':     ['g'],          # Ca-dependent K channel (K ion)
    'kdrRL':    ['gMax'],       # Delayed Rectifier K Channel (K ion)
    'L_Ca':     ['gcabar'],     # L-type Ca channel (Ca ion, CaL virtual ion)
    'mAHP':     ['gkcamax', 'gcamax'], # Ca-dependent K channel + Ca channel (Ca + K ions)
    'na3rp':    ['gbar'],       # Fast Na current (Na ion)
    'naps':     ['gbar'],       # Persistent Na current (Na ion)
    'pas':      ['g'],      # Passive/leak channel
}
gleak_name = 'g_pas'

# Mechanism parameters that are changed from default values in original model code
mechs_params_dict = {
    'gh':       ['ghbar', 'htau', 'half', 'slp'], # H channel (Na + K ions)
    'kca2':     ['g', 'depth2', 'taur2'], # Ca-dependent K channel (K ion)
    'kdrRL':    ['gMax', 'mVh'], # Delayed Rectifier K Channel (K ion)
    'L_Ca':     ['gcabar', 'tau_m', 'theta_m'], # L-type Ca channel (Ca ion, CaL virtual ion)
    'mAHP':     ['gkcamax', 'gcamax', 'taur'], # Ca-dependent K + Ca channel (Ca + K ions)
    'na3rp':    ['gbar', 'sh', 'ar', 'Rd', 'qd', 'qg', 'thi1', 'thi2'], # Fast Na current (Na ion)
    'naps':     ['gbar', 'sh', 'ar'], # Persistent Na current (Na ion)
    'pas':      ['g', 'e'],  # Passive/leak channel
    'hh':       ['gnabar', 'gkbar', 'gl', 'el'], # Hodgkin-Huxley mechanisms
}

# All mechanism parameters that are not conductances
mechs_params_nogbar = dict(mechs_params_dict)
for mech, params in mechs_params_nogbar.items():
    for gbar_param in gbar_dict.get(mech, []):
        try:
            params.remove(gbar_param)
        except ValueError:
            pass

# GLOBAL mechanism parameters (assigned using h.param = val)
global_params_list = [
    'tmin_kdrRL', 'taumax_kdrRL',
    'qinf_na3rp', 'thinf_na3rp',
    'vslope_naps'
]

mechs_list = list(mechs_params_dict.keys()) # all mechanisms
gbar_list = [gname+'_'+mech for mech,chans in gbar_dict.items() for gname in chans]
active_gbar_names = [gname for gname in gbar_list if gname != gleak_name]


def make_cell_balbi(model_no=1):
    """
    Initialize Balbi et al. cell model

    @param  model_no    model number to load: integer in range 1-14

    @return             dict { region_name<str> : list(Section) }

    @effect             Following variables will be available on Hoc interpreter:

                        * somatic sections:
                            - soma[N] <Section> somatic sections

                        * dendritic sections:
                            - dend[M] <Section> dendritic sections
                        
                        * axon-soma interface
                            - AH[1] <Section> axon hillock
                            - IS[1] <Section> axon initial segment
                        
                        * axonal sections:
                            - node[axonnodes] <Section>
                            - MYSA[paranodes1] <Section>
                            - FLUT[paranodes2] <Section>
                            - STIN[axoninter] <Section>
    """
    varnames_used = model_secarray_vars + model_sec_vars
    if any((hasattr(h, attr) for attr in varnames_used)):
        raise Exception("Folowing global variables must be unallocated on Hoc interpreter object: {}".format(', '.join(varnames_used)))

    # load model
    h.xopen("createcell_balbi.hoc")
    h.choose_model(model_no)


def get_named_sec_lists():
    """
    Get one Python list per section group.

    @NOTE       node[0] is deleted in file '2_complete_cell.hoc': this crashes
                NEURON when trying to access

    @return     <dict(str:list(h.Section))> with section array names as keys,
                and their member sections stored in Python list as values
    """
    # NOTE: node[0] is deleted in file '2_complete_cell.hoc'
    sec_arrays = list(model_secarray_vars)
    sec_arrays.remove('node')
    named_sections = { arrayname: list(getattr(h, arrayname)) for arrayname in sec_arrays}
    named_sections.update({ secname: [getattr(h, secname)] for secname in model_sec_vars})

    # Add existing sections named 'node'
    node_sec_ids = range(1, int(h.axonnodes)) # node[0] does not exist
    named_sections['node'] = [h.node[i] for i in node_sec_ids]

    return named_sections


def motocell_steadystate(model_no):
    """
    Load steate-state condition for motoneuron cell model with given cell id

    @param  model_no    model number to load: integer in range 1-14

    @note   h.load_steadystate() uses SaveState.restore() which means 
            Between a save and a restore, it is important not to create or delete sections, NetCon objects, or point processes. Do not change the number of segments, insert or delete mechanisms, or change the location of point processes.
    """
    h.celsius = 37
    h.load_steadystate(model_no)


if __name__ == '__main__':
    named_sec = make_cell_balbi(model_no=1)
    # from neuron import gui