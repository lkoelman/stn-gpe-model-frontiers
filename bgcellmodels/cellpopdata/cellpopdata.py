"""
Parameters for STN connections

@author     Lucas Koelman
@date       12/07/2017


Based on structure of:

    bgmodel/models/kang/nuclei/cellpopdata.py
    bgmodel/models/kumaravelu/nuclei/cellpopdata.py
    bgmodel/models/netkang/kang_netparams.py

"""

from enum import Enum, unique
from textwrap import dedent
import collections

import numpy as np
import neuron
nrn = neuron.nrn # types nrn.Section and nrn.Segment
h = neuron.h

from physiotypes import PhysioState, Populations, NTReceptors, ParameterSource
Pop = Populations
Rec = NTReceptors
Src = ParameterSource

import synmechs

from bgcellmodels.common.configutil import interpretParamSpec
from bgcellmodels.common import logutils

logger = logutils.getBasicLogger(__name__, 'DEBUG')


@unique
class StnModel(Enum):
    """
    STN cell models
    """
    Gillies2005 = 0
    Miocinovic2006 = 1
    Gillies_GIF = 2

    # reduction using Gillies-specific code
    Gillies_BranchZip_Legacy = 3
    Gillies_FoldStratford_Legacy = 4
    Gillies_FoldMarasco_Legacy = 5

    # reduction using CERSEI toolbox
    Gillies_FoldMarasco_Tapered = 6
    Gillies_FoldBush_Tapered = 7


# Indicate which of the STN models are reduced models
ReducedModels = list(StnModel)
ReducedModels.remove(StnModel.Gillies2005)


class CellConnector(object):
    """
    Class for storing connection parameters and making connections.
    """

    def __init__(self, physio_state, rng, 
            preferred_mechanisms=None, preferred_sources=None):
        """

        @param      preferred_mechanisms : list(str)
                    
                    Ordered list of preferred synaptic mechanisms to use
                    when making connections. When connection is made with only
                    the receptor specified, mechanisms will be used in this
                    order.
        """
        if isinstance(physio_state, str):
            physio_state = PhysioState.from_descr(physio_state)
        
        self._physio_state = physio_state
        self._rng = rng

        # make some functions available as instance methods (lazy refactoring)
        self.getSynMechReceptors = synmechs.get_mechanism_receptors
        self.getSynMechParamNames = synmechs.get_mech_param_names
        self.getNrnConParamMap = synmechs.get_physio2mech_map

        if preferred_mechanisms is None:
            self.preferred_mechanisms = ['GLUsyn', 'GABAsyn', 'Exp2Syn']
        else:
            self.preferred_mechanisms = preferred_mechanisms

        if preferred_sources is None:
            self.preferred_sources = [Src.Custom, Src.Default]
        else:
            self.preferred_sources = preferred_sources


    def getFireParams(self, pre_pop, phys_state, use_sources, custom_params=None):
        """
        Get parameters describing firing statistics for given population.
        
        @param use_sources      ordered list of ParameterSource members
                                to indicate which literatur sources should
                                be used.
        """

        firing_params = dict(((pop, {}) for pop in list(Populations)))
        fp = firing_params

        # ---------------------------------------------------------------------
        # parameters from Bergman (2015)
        fp[Pop.GPE][Src.Bergman2015RetiCh3] = {}

        fp[Pop.GPE][Src.Bergman2015RetiCh3][PhysioState.NORMAL | PhysioState.AWAKE] = {
            'rate_mean': 60.0, # 50-70 Hz
            'rate_deviation': 10.0,
            'rate_units': 'Hz',
            'pause_dur_mean': 0.6, # average pause duration 0.5-0.7 s in monkeys
            'pause_dur_units': 's',
            'pause_rate_mean': 15.0/60.0, # 10-20 pauses per minute
            'pause_rate_units': 'Hz',
            'pause_rate_dist': 'poisson', # Can control rate NetStim with pause NetStim (stim.number resets, see mod file)
            # NOTE: discarge_dur = 1/pause_rate - pause_dur ~= 4 - 0.6
            'discharge_dur_mean': 1.0, # must be < 1/pause_rate_mean
            'discharge_dur_units': 's',
            'species': 'primates',
        }


        fp[Pop.CTX][Src.Bergman2015RetiCh3] = {}

        fp[Pop.CTX][Src.Bergman2015RetiCh3][PhysioState.NORMAL | PhysioState.AWAKE] = {
            'rate_mean': 2.5,
            'rate_deviation': 2.5,
            'rate_dist': 'poisson',
            'rate_units': 'Hz',
            'species': 'primates',
        }

        # ---------------------------------------------------------------------
        # parameters from Mallet (2016)
        fp[Pop.GPE][Src.Mallet2016] = {}

        fp[Pop.GPE][Src.Mallet2016][PhysioState.NORMAL | PhysioState.AWAKE] = {
            'rate_mean': 47.3,
            'rate_deviation': 6.1,
            'rate_units': 'Hz',
            'species': 'rats',
        }

        # ---------------------------------------------------------------------
        # parameters from Nambu (2014)
        fp[Pop.GPE][Src.Nambu2014] = {}

        fp[Pop.GPE][Src.Nambu2014][PhysioState.NORMAL | PhysioState.AWAKE] = {
            'rate_mean': 65.2,
            'rate_deviation': 25.8,
            'rate_units': 'Hz',
            'species': 'monkey',
            'subspecies': 'macaque',
        }

        # Get final parameters
        fp_final = {} # {Population -> params}
        use_sources = list(use_sources)


        # Use preferred sources to update final parameters
        for citation in reversed(use_sources):

            if citation in firing_params[pre_pop]:

                # Update params if requested Physiological state matches the described state
                for state_described, params in firing_params[pre_pop][citation].iteritems():
                    if phys_state.is_subset(state_described):
                        fp_final.update(params)

            elif (citation == Src.Custom) and (custom_params is not None):
                # If custom params provided, and custom is in list of preferred params: use them
                fp_final.update(custom_params)

        return fp_final

    def get_physiological_parameters(
            self,
            pre_pop,
            post_pop,
            use_sources=None,
            custom_params=None,
            adjust_gsyn_msr=None):
        """
        Get parameters for afferent connections onto given population,
        in given physiological state.

        @param  custom_params : dict<NTR, dict<str, object>>

                Custom parameters in the form of a dict
                {
                    NTR_0: {params...},
                    NTR_1: {params...},
                    ...
                }


        @param  adjust_gsyn_msr : int

                Number of synapses per axonal contact (int), 
                and whether the synaptic condcutance for each synapse
                should be scaled to take this into account.
                If an int > 0 is given, each synaptic condcutance
                is divided by this number.


        @return params : dict<NTR, dict<str, object>>
                
                Dictionary containing synapse parameters for each neurotransmitter
                involved in the given connection from pre to post-synaptic population.
                The parameter dict for each key (NTReceptor) may contain following
                entries:

                'Ipeak':        float       peak synaptic current
                'Epeak':        float       peak PSP amplitude
                'PSP_dV':       float       PSP deviation from baseline membrane voltage
                'gbar':         float       synaptic conductance
                'tau_rise_g':   float       exponential time constant for rising phase
                'tau_decay_g':  float       exponential rime constant for decay phase
                'Erev':         float       reversal potential for synapse
                'tau_rec_STD':  float       time constant of recovery from depression
                'tau_rec_STP':  float       time constant of facilitation
                'P_release_base':   float   baseline release probability


        """

        physio_state = self._physio_state
        rng = self._rng

        if rng is None:
            rng = np.random

        # Initialize parameters dict
        cp = {}

        # TODO: refactor this in table (excel-like): one table per projection,
        #       one row per observation, one column per parameter (and NT, comment).
        #       Goal should be that you can clearly see which information is missing,
        #       which is N/A, and what is known about the connection
        
        # TODO: add units and check them while setting, like Ephys parameters.
        #       Also plot a final connectivity matrix that shows all connections
        #       strengths in same units (easy to check visually, in color)

        # TODO: none of the time constants are temperature-corrected

        # TODO: remove IPSC/EPSC calculations: only save data that is
        #       actually present/shown in articles (PSP magnitudes etc.)
        if post_pop == Pop.GPE:

            # Inputs from Striatum
            cp[Pop.STR] =  {
                Rec.GABAA: {},
                Rec.GABAB: {},
            }

            # Inputs from Subthalamic Nucleus
            cp[Pop.STN] =  {
                Rec.AMPA: {},
                Rec.NMDA: {},
            }

            # ---------------------------------------------------------------------
            # Default parameters XYZ -> GPE.
            # TODO: Change these. Current values are copied from AMPA/NMDA CTX->STN
            
            # Default GLUergic (STN -> GPE)
            Erev = 0.0
            Ermp = -80.
            cp[Pop.STN][Rec.AMPA][Src.Default] = {
                'Ipeak': -275.,     # peak synaptic current (pA)
                'gbar': -275.*1e-3 / (Ermp - Erev), # gbar calculation (~ 3e-3 uS)
                'tau_rise_g': 1.0,
                'tau_decay_g': 4.0,
                'Erev': Erev,
                'tau_rec_STD': 200.,
                'tau_rec_STP': 1., # no facilitation
                'P_release_base': 0.7,
            }

            # NMDA from Chu (2015)
            Ermp = 30. # NMDA EPSC measured at 30 mV rmp to remove Mg2+ block
            I_peak = 210. # opposite sign: see graph
            cp[Pop.STN][Rec.NMDA][Src.Default] = {
                'Ipeak': I_peak,    # peak synaptic current (pA)
                'gbar': I_peak * 1e-3 / (Ermp - Erev), # gbar calculation
                'tau_rise_g': 3.7,
                'tau_decay_g': 80.0,
                'Erev': Erev,
                'tau_rec_STD': 200.,
                'tau_rec_STP': 1., # no facilitation
                'P_release_base': 0.7,
            }

            # Default GABAergic (STR -> GPE)
            Erev_GABAA = -85.
            Erev_GABAB = -93.
            Ermp = -70.

            cp[Pop.STR][Rec.GABAA][Src.Default] = {
                'gbar': 350.*1e-3 / (Ermp - Erev_GABAB), # Chu2015 in healthy state
                'Erev': Erev_GABAA,
                'delay': 1.0,
                'Vpre_threshold': 0.0,
                'f_med_PSP_single': 0.19, # median frequency of PSP triggered by burst
                'f_med_PSP_burst': 0.05, # median frequency of PSP triggered by single spike
            }

            cp[Pop.STR][Rec.GABAB][Src.Default] = {
                'Erev': Erev_GABAB,
                'gbar': 350.*1e-3 / (Ermp - Erev_GABAB), # Chu2015 in healthy state
                'delay': 1.0,
                'Vpre_threshold': 0.0,
                'f_med_PSP_single': 4.22,
                'f_med_PSP_burst': 8.72,
            }

            # ------------------------------------------------------------------
            # Levine, Hull, Buchenwald (1974) GP PSPs
            # Data is from cats

            # Figure 1: measurement of IPSP by stimulating Caudate Nucleus
            cp[Pop.STR][Rec.GABAB][Src.Levine1974] = {
                'PSP_dV': -5.0, # ~ 5 mV dip from baseline
                'tau_rise_g': 66.0, # ~ 66 ms to 63% of dip amplitude
                'tau_decay_g': 100.0, # ~ 100 ms to 63% recovery to baseline
            }

            # Figure 5: measurement of EPSP-IPSP by stimulating Caudate
            # This is the magnitude of the initial EPSP before the IPSP
            cp[Pop.STR][Rec.AMPA][Src.Levine1974] = {
                'PSP_dV': 2.65, # ~ +2.6 mV increase from baseline
                'tau_rise_g': 23.5, # ~ 23.5 ms to 63% of peak
                'tau_decay_g': 29.5, # ~ 29.5 ms to 63% recovery to baseline
            }

            # ------------------------------------------------------------------
            # Kita & Kitai (1991) GP intracellular recordings
            # Data is from rats

            # Figure 7: measurement of EPSP by stimulating STh (Subthalamic Nucleus)
            cp[Pop.STN][Rec.AMPA][Src.KitaKitai1991] = {
                'PSP_dV': (11.0, 31.0), # ~ 21 mV increase for middle trace
                'tau_rise_g': 12.0, # ~ 12 ms to 63% of peak
                'tau_decay_g': 30.0, # ~ 30 ms to 63% recovery to baseline
            }


            # ------------------------------------------------------------------
            # Parameters Hendrickson (2011) - Capabilities and Limitations ...
            # TODO: chech these, STN -> GPe should be excitatory, and STR->inhibitory

            # cp[Pop.STR][Rec.AMPA][Src.Hendrickson2011] = {
            #     'gbar':         0.25,
            #     'tau_rise_g':   1.0,
            #     'tau_decay_g':  3.0,
            #     'Erev':         0.0,
            # }

            # cp[Pop.STR][Rec.NMDA][Src.Hendrickson2011] = {
            #     'gbar':         0.25,
            #     'tau_rise_g':   10.0,
            #     'tau_decay_g':  30.0,
            #     'Erev':         0.0,
            # }

            # cp[Pop.STN][Rec.GABAA][Src.Hendrickson2011] = {
            #     'gbar':         0.25,
            #     'tau_rise_g':   1.0,
            #     'tau_decay_g':  12.0,
            #     'Erev':         -80.0,
            # }

        if post_pop == Pop.STN:

            cp[Pop.CTX] =  {
                Rec.AMPA: {},
                Rec.NMDA: {},
            }

            cp[Pop.GPE] = {
                Rec.GABAA: {},
                Rec.GABAB: {},
            }

            # TODO: correct both gmax/Imax for attenuation from dendrites to soma. 
            #       Do this separately for GPe and CTX inputs since they have to travel different path lengths.

            # TODO SETPARAM: make sure default delay, Vpre etc. is set for each (pre,post,rec) combination

            # gmax calculation:
            # gmax is in [uS] in POINT_PROCESS synapses
            #   1 [uS] * 1 [mV] = 1 [nA]
            #   we want ~ -300 [pA] = -0.3 [nA]
            #   gmax [uS] * (Ermp [mV] - Erev [mV]) = Isyn [nA]
            #
            #       => gmax [uS] = Isyn [nA] / (Ermp-Erev) [mV]

            #######################################################################
            # CTX -> STN parameters
            #######################################################################

            Erev = 0.

            # ---------------------------------------------------------------------
            # Default parameters
            for ntr in cp[Pop.CTX].keys():
                cp[Pop.CTX][ntr][Src.Default] = {
                    'delay': 1.0,
                    'Vpre_threshold': 0.0,
                }

            cp[Pop.CTX][Rec.AMPA][Src.Default].update({
                'f_med_PSP_single': 16.87,
                'f_med_PSP_burst': 16.95,
            })

            cp[Pop.CTX][Rec.NMDA][Src.Default].update({
                'f_med_PSP_single': 0.58,
                'f_med_PSP_burst': 2.14,
            })

            # ---------------------------------------------------------------------
            # Parameters Chu (2015) - Heterosynaptic Regulation of External Globus
            # Pallidus Inputs to the Subthalamic Nucleus by the Motor Cortex

            # AMPA from Chu (2015)
            Ermp = -80.
            cp[Pop.CTX][Rec.AMPA][Src.Chu2015] = {
                'Ipeak': -275.,     # peak synaptic current (pA)
                'gbar': -275.*1e-3 / (Ermp - Erev), # gbar calculation (~ 3e-3 uS)
                'tau_rise_g': 1.0,
                'tau_decay_g': 4.0,
                'Erev': Erev,
            }

            # Changes in parkinsonian state
            if physio_state == PhysioState.PARKINSONIAN:
                Ermp = -80. # Fig. 2G
                cp[Pop.CTX][Rec.AMPA][Src.Chu2015].update({
                    'Ipeak': -390., # pA
                    'gbar': -390.*1e-3 / (Ermp - Erev), 
                })
            park_gain = 390./275. # increase in Parkinsonian condition

            # NMDA from Chu (2015)
            Ermp = 30. # NMDA EPSC measured at 30 mV rmp to remove Mg2+ block
            I_peak = 210. # opposite sign: see graph
            cp[Pop.CTX][Rec.NMDA][Src.Chu2015] = {
                'Ipeak': I_peak,    # peak synaptic current (pA)
                'gbar': I_peak * 1e-3 / (Ermp - Erev), # gbar calculation
                'tau_rise_g': 3.7,
                'tau_decay_g': 80.0,
                'Erev': Erev,
            }

            # Changes in parkinsonian state
            if physio_state == PhysioState.PARKINSONIAN:
                Ermp = -80. # Fig. 2G
                # NOTE: increase peak conductance by same factor as AMPA
                cp[Pop.CTX][Rec.NMDA][Src.Chu2015].update({
                    'Ipeak': park_gain * cp[Pop.CTX][Rec.NMDA][Src.Chu2015]['Ipeak'],
                    'gbar':  park_gain * cp[Pop.CTX][Rec.NMDA][Src.Chu2015]['Ipeak'],   
                })

            # ---------------------------------------------------------------------
            # Parameters Gradinaru (2009)

            cp[Pop.CTX][Rec.AMPA][Src.Gradinaru2009] = {
                'tau_rec_STD': 200.,
                'tau_rec_STP': 1., # no facilitation
                'P_release_base': 0.7,
            }


            # ---------------------------------------------------------------------
            # Default parameters

            # Set params Chu (2015) as default parameters
            cp[Pop.CTX][Rec.AMPA][Src.Default].update(cp[Pop.CTX][Rec.AMPA][Src.Chu2015])
            cp[Pop.CTX][Rec.NMDA][Src.Default].update(cp[Pop.CTX][Rec.NMDA][Src.Chu2015])
            
            # Modification to NMDA conductance
            #   -> NMDA conductance is typically 70% of that of AMPA (see EPFL MOOC)
            cp[Pop.CTX][Rec.NMDA][Src.Default]['gbar'] = 0.7 * cp[Pop.CTX][Rec.AMPA][Src.Default]['gbar']
            
            
            
            #######################################################################
            # GPe -> STN parameters
            #######################################################################

            # ---------------------------------------------------------------------
            # Default parameters
            Erev_GABAA = -85.
            Erev_GABAB = -93.
            Ermp = -70.

            cp[Pop.GPE][Rec.GABAA][Src.Default] = {
                'Erev': Erev_GABAA,
                'delay': 1.0,
                'Vpre_threshold': 0.0,
                'f_med_PSP_single': 0.19, # median frequency of PSP triggered by burst
                'f_med_PSP_burst': 0.05, # median frequency of PSP triggered by single spike
            }

            cp[Pop.GPE][Rec.GABAB][Src.Default] = {
                'Erev': Erev_GABAB,
                'gbar': 350.*1e-3 / (Ermp - Erev_GABAB), # Chu2015 in healthy state
                'delay': 1.0,
                'Vpre_threshold': 0.0,
                'f_med_PSP_single': 4.22,
                'f_med_PSP_burst': 8.72,
            }

            # ---------------------------------------------------------------------
            # Parameters from Chu (2015)

            # GABA-A ISPSc examples in Figures 2,3,4,5,6.
            cp[Pop.GPE][Rec.GABAA][Src.Chu2015] = {
                'Ipeak': 350.,  # peak synaptic current (pA)
                'gbar': 350.*1e-3 / (Ermp - Erev_GABAA), # gbar calculation
                'tau_rise_g': 2.6,
                'tau_decay_g': 5.0,
            }
            if physio_state == PhysioState.PARKINSONIAN:
                cp[Pop.GPE][Rec.GABAA][Src.Chu2015].update({
                    'Ipeak': 450.,
                    'gbar': 450.*1e-3 / (Ermp - Erev_GABAA),
                    'tau_rise_g': 3.15,
                    'tau_decay_g': 6.5,
                })

            # ---------------------------------------------------------------------
            # Parameters from Fan (2012)

            # Figure 2 shows magnitudes for spike-invoked IPSCs
            cp[Pop.GPE][Rec.GABAA][Src.Fan2012] = {
                'gbar': {
                    'mean': 7.03e-3,
                    'deviation': 3.10e-3,
                    'units': 'uS'
                    },
                'tau_rise_g': 1.0,
                'tau_decay_g': 6.0,
            }
            if physio_state == PhysioState.PARKINSONIAN:
                cp[Pop.GPE][Rec.GABAA][Src.Fan2012].update({
                    'gbar': {
                        'mean': 11.17e-3,
                        'deviation': 5.41e-3,
                        'units': 'uS'
                        },
                    'tau_rise_g': 1.0,
                    'tau_decay_g': 8.0,
                })

            # ---------------------------------------------------------------------
            # Parameters from Atherton (2013)

            # Figure 2 shows short-term depression
            cp[Pop.GPE][Rec.GABAA][Src.Atherton2013] = {
                'tau_rec_STD': 17300., # See Fig. 2 & legend
                'tau_rec_STP': 1.0, # no STP: very fast recovery
                'P_release_base': 0.5, # Fitted
            }

            # ---------------------------------------------------------------------
            # Parameters from Baufreton (2009)
            cp[Pop.GPE][Rec.GABAA][Src.Baufreton2009] = {
                'Ipeak': {
                        'min': 375.,
                        'max': 680.,
                        'units': 'pA',
                    },
                'gbar': 530.*1e-3 / (Ermp - Erev_GABAA), # gbar calculation
                'tau_rise_g': 2.6,
                'tau_decay_g': 5.0,
            }

        # Make final dict with only (receptor -> params)
        cp_final = {}
        if use_sources is None:
            use_sources = self.preferred_sources
        use_sources = list(use_sources)
        use_sources.append(Src.Default) # this source was used for default parameters
        use_sources.append(Src.CommonUse) # source used for unreferenced parameters from other models

        for receptor in cp[pre_pop].keys(): # all receptors involved in this connection
            cp_final[receptor] = {}

            # Use preferred sources to update final parameters
            for citation in reversed(use_sources):
                
                if ((citation == Src.Custom) and (custom_params is not None)
                    and (receptor in custom_params)):
                    ntr_params = custom_params[receptor]

                elif (citation in cp[pre_pop][receptor]):
                    ntr_params = cp[pre_pop][receptor][citation]

                else:
                    continue # citation doesn't provide parameters about this receptor

                # Successively overwrite with params of each source
                cp_final[receptor].update(ntr_params)

            # Adjust for multi synapse rule
            if (adjust_gsyn_msr) and ('gbar' in cp_final[receptor]):
                cp_final[receptor]['gbar'] = cp_final[receptor]['gbar']/adjust_gsyn_msr

        return cp_final

    # alias
    getPhysioConParams = get_physiological_parameters


    def get_nrn_from_physio_params(self, nrn_mech_name, physio_params):
        """
        Get parameters to assign to NEURON objects.

        @return     dictionary {nrn_obj_type: {attr_name: value} }
                    
                    where nrn_obj_type is one of:

                        - 'pointprocess',
                        - 'netcon'
                        - 'netstim'

                    and the inner dictionary are parameter names and values for 
                    these object types.

        USAGE:

            physio_params = cc.get_physiological_parameters(pre_pop, post_pop, use_sources)
            mech = 'Exp2Syn'
            nrn_params = cc.get_nrn_from_physio_params(mech, physio_params)
        
        """
        # Get NT Receptors used in connection
        receptors = physio_params.keys()

        # Get mapping from physiological to NEURON mechanism parameters
        nrn_param_map = synmechs.get_physio2mech_map(nrn_mech_name)

        # keep track of parameters that are assigned
        nrn_param_values = {
            'pointprocess': {},
            'netcon':       {},
            'netstim':      {},
        }

        # For each NT receptor, look up the physiological connection parameters,
        # and translate them to a parameter of the synaptic mechanism or NetCon
        for rec in receptors:

            physio_to_nrn = nrn_param_map[rec] 
            physio_to_values = physio_params[rec]

            # Translate each physiological parameter to NEURON parameter
            for physio_name, nrn_param_spec in physio_to_nrn.iteritems():

                # Check if parameters is available from given sources
                if physio_name not in physio_to_values:
                    logger.anal("Parameter {dictpar} not found for connection. "
                                "This means that parameter {mechpar} will not be set\n".format(
                                dictpar=physio_name, mechpar=nrn_param_spec))
                    continue

                mech_type, param_name, param_index = interpretParamSpec(nrn_param_spec)

                # Get the actual parameter value
                value_spec = physio_to_values[physio_name]
                value = synmechs.evalValueSpec(value_spec, self._rng)

                # Determine target (synapse or NetCon)
                if mech_type == 'syn':
                    nrn_dict = nrn_param_values['pointprocess']
                else:
                    nrn_dict = nrn_param_values[mech_type]

                # Set attribute
                if param_index is None:
                    nrn_dict[param_name] = value
                else:
                    values = nrn_dict.get(param_name, None)
                    if values is None:
                        values = [0.0] * (param_index+1)
                        nrn_dict[param_name] = values
                    elif len(values) <= param_index:
                        values.extend([0.0] * (param_index+1-len(values)))
                    values[param_index] = value

        # Apply possible corrections to synaptic parameters
        correct_syn = synmechs.syn_mech_correctors.get(nrn_mech_name, None)
        if correct_syn is not None:
            correct_syn(nrn_param_values)

        return nrn_param_values


    # alias
    getNrnObjParams = get_nrn_from_physio_params


    def get_synaptic_mechanism(self, receptors):
        """
        Get a suitable mechanism that implements the given receptors

        @param      receptors : iterable(str) / iterable(NTReceptor)

        @return     mechanism_name : str
        """
        # Get the synaptic mechanism name
        suitable_mechanism = None
        for mech in self.preferred_mechanisms:
            supported_receptors = synmechs.get_mechanism_receptors(mech)
            required_receptors = (NTReceptors.from_str(rec) for rec in receptors)
            
            # Check that the mechanism implements all receptors
            if all(((rec in supported_receptors) for rec in required_receptors)):
                suitable_mechanism = mech

        if suitable_mechanism is None:
            raise ValueError("Could not find synaptic mechanism that "
                    "implements all receptors in {}".format(receptors))

        return suitable_mechanism



    def make_empty_synaptic_connection(self, pre_post_pop, pre_post_obj, syn_mech):
        """
        Create the synaptic POINT_PROCESS and NetCon instance without
        setting their parameters.


        ARGUMENTS
        ---------
        
        @param      pre_post_pop : tuple(str, str)
                    Identifiers for pre- and post-synaptic populations, e.g.
                    (Populations.GPE, Populations.STN)

        @param      pre_post_obj : tuple(object, object)
                    Source and target object for synaptic connection. 
                    Synapse will be inserted into target object
        """
        pre_pop, post_pop = pre_post_pop
        pre_obj, post_obj = pre_post_obj

        if not isinstance(post_obj, nrn.Segment):
            raise ValueError("Post-synaptic object {} is not of type nrn.Segment".format(repr(post_obj)))

        # Get synapse type constructor and make it
        syn_ctor = synmechs.syn_wrapper_classes.get(syn_mech, getattr(h, syn_mech))
        syn = syn_ctor(post_obj)

        # Store some data on the synapse
        if hasattr(syn, 'htype'): # check if wrapper class (see nrn/lib/python/hclass.py)
            syn.pre_pop = pre_pop.name
            syn.post_pop = post_pop.name

        # Make NetCon connection
        if isinstance(pre_obj, nrn.Section):
            # Biophysical cells need threshold detection to generate events
            nc = synmechs.NetConWrapper(pre_obj(0.5)._ref_v, syn, sec=pre_obj)
        else:
            # Source object is POINT_PROCESS or other event-generating objcet
            nc = synmechs.NetConWrapper(pre_obj, syn)
        nc.pre_pop = pre_pop.name
        nc.post_pop = post_pop.name

        return syn, nc


    def make_synaptic_connection_from_DB(
            self, pre_post_pop, pre_post_obj, syn_type, receptors, 
            use_sources=None, custom_conpar=None, custom_synpar=None,
            con_par_data=None, weight_scales=None, weight_times=None):
        """
        Insert synapse POINT_PROCESS in given section.


        @param      custom_conpar : dict(NTR, dict)
                    
                    Custom physiological parameters to override those from DB,
                    in the form of a dict {NTR_0: {params_0}, NTR_1: {params_1}}


        @see        set_connection_params_from_physio()
        
        """

        pre_pop, post_pop = pre_post_pop

        # Initialize an empty synaptic connection
        syn, nc = self.make_empty_synaptic_connection(
                    pre_post_pop, pre_post_obj, syn_type)

        # Get physiological parameter descriptions
        if con_par_data is None:
            con_par_data = self.get_physiological_parameters(
                            pre_pop, post_pop, use_sources, custom_conpar)

        # Set connection parameters from physiological parameters + possible
        # overriding mechanism parameters
        weight_vecs = self.set_connection_params_from_physio(
                        syn, nc, syn_type, receptors, con_par_data, 
                        custom_synpar, weight_scales, weight_times)


        # Return refs to objects that need to stay alive
        return syn, nc, weight_vecs

    # Alias for backward compatibility
    make_synapse = make_synaptic_connection_from_DB


    def set_connection_params_from_physio(
            self, syn, nc, syn_type, receptors, 
            con_par_data=None, custom_synpar=None,
            weight_scales=None, weight_times=None):
        """
        Set parameters of synapse POINT_PROCESS and NetCon from database.

        USAGE
        -----

        - either provide argument 'use_sources' to fetch connection parameters
          or provide them yourself in argument 'con_par_data'

        
        ARGUMENTS
        ---------

        @param      con_par_data : dict
                    
                    The physiological parameters.
                    Either dict<NTReceptors, dict<str, object>> if multiple receptors
                    are used on the given synapse mechanism, or a dict<str, object>
                    with synapse parameters if only a single receptor is used.


        @param      custom_synpar : dict 

                    Custom mechanism parameters for the synaptic mechanism,
                    in the form of a dict {param_name: param_value, ...}


        @param      weight_scales : iterable

                    Scale factors for the weights in range (0,1). 
                    This must be an iterable: NetCon.weight[i] is scaled by the i-th value. 
                    If a vector is given, it is scaled and played into the weight.

        
        """
        
        # Get mapping from physiological to NEURON mechanism parameters
        syn_par_map = synmechs.get_physio2mech_map(syn_type)

        # keep track of parameters that are assigned
        syn_assigned_pars = {}
        netcon_assigned_pars = {}
        assigned_params = {
            'syn': syn_assigned_pars,
            'netcon': netcon_assigned_pars
        }

        # For each NT receptor, look up the physiological connection parameters,
        # and translate them to a parameter of the synaptic mechanism or NetCon
        num_receptors = len(receptors)
        for receptor_descr in receptors:
            rec = NTReceptors.get(receptor_descr)
            parname_map = syn_par_map[rec] # how each connection parameter is mapped to mechanism parameter

            # Get dict with synapse param values
            if num_receptors==1 and rec not in con_par_data.keys():
                phys_params = con_par_data
            else:
                phys_params = con_par_data[rec] # physiological parameters from given sources

            # Translate each param value to a mechanism parameter
            for phys_parname, mech_parspec in parname_map.iteritems():

                # Check if parameters is available from given sources
                if phys_parname not in phys_params:
                    # logger.anal(dedent("""\
                    #        Parameter {dictpar} not found for connection ({pre},{post},{rec}).
                    #        This means that parameter {mechpar} will not be set\n""").format(
                    #            dictpar=phys_parname, pre=pre_pop, post=post_pop, 
                    #            rec=rec, mechpar=mech_parspec))
                    continue

                # Interpret parameter specification (what is the parameter?)
                mechtype, mech_parname, paridx = interpretParamSpec(mech_parspec)

                # Get the actual parameter value
                value_spec = phys_params[phys_parname]
                par_val = synmechs.evalValueSpec(value_spec, self._rng)

                # Determine target (synapse or NetCon)
                if mechtype == 'syn':
                    target = syn
                
                elif mechtype == 'netcon':
                    target = nc
                else:
                    raise ValueError("Cannot set attribute of unknown mechanism type {}".format(mechtype))

                # Set attribute
                if paridx is None:
                    setattr(target, mech_parname, par_val)
                else:
                    getattr(target, mech_parname)[int(paridx)] = par_val

                # Save which parameters were assigned
                assigned_params[mechtype][mech_parname] = par_val

        # Custom synaptic mechanism parameters
        if custom_synpar is not None:
            for pname, pval in custom_synpar.iteritems():
                setattr(syn, pname, pval)

        # Apply possible corrections to synaptic parameters
        correct_syn = synmechs.syn_mech_correctors.get(syn_type, None)
        if correct_syn is not None:
            correct_syn(syn)

        # Set weights
        weight_vecs = []
        if (weight_scales is None) and ('weight' not in netcon_assigned_pars):
            # Weight was never assigned, assume synapse has gmax attribute
            nc.weight[0] = 1.0

        elif (weight_scales is not None):
            for i_w, weight in enumerate(weight_scales):

                if isinstance(weight, (int, float)):
                    # If weight was assigned as part of synapse parameter, scale it
                    if ('weight' in netcon_assigned_pars):
                        nc.weight[i_w] = nc.weight[i_w] * weight
                    else:
                        nc.weight[i_w] = weight

                else: # weight is h.Vector
                    # If weight was assigned as part of synapse parameter, scale it
                    if ('weight' in netcon_assigned_pars):
                        logger.anal("Weight was assigned as part of synapse parameters. New weight Vector will be allocated.")
                        wvec = h.Vector()
                        wvec.copy(weight)
                        wvec.mul(nc.weight[i_w]) # Vector.mul() : multiply in-place
                    else:
                        wvec = weight

                    # Play Vector into weight
                    timevec = weight_times[i_w]
                    wvec.play(nc._ref_weight[i_w], timevec)
                    weight_vecs.append(wvec)

        # Return refs to objects that need to stay alive
        return syn, nc, weight_vecs