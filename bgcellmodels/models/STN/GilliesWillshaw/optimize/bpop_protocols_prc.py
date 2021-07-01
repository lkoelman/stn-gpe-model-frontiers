"""
Stimulation protocols for measuring Phase Response Curves (PRC)
with BluePyOpt.

@author Lucas Koelman
@date   16/10/2019


USAGE
-----

- see comments marked 'SETPARAM' for free parameters

"""

# Standard library
import logging

# Third party
import numpy as np
import bluepyopt.ephys as ephys

# Custom modules
from bgcellmodels.extensions.bluepyopt import bpop_recordings, bpop_stimuli
from bgcellmodels.extensions.bluepyopt.bpop_protocol_ext import (
    SelfContainedProtocol, PhysioProtocol, BpopProtocolWrapper,
    rng_getter, connector_getter, PROTOCOL_WRAPPERS
)
from bgcellmodels.common import stimprotocols

# Module globals
logger = logging.getLogger('prc_protos')
StimProtocol = stimprotocols.StimProtocol


def init_stn_physiology(sim, model):
    """
    Initialize simulator to run plateau protocol

    NOTE: function must be declared at top-level of module in order to be pickled
    """
    h = sim.neuron.h
    h.celsius = 30
    h.v_init = -60
    h.set_aCSF(4)
    h.init()


loc_soma_center = ephys.locations.NrnSeclistCompLocation(
        name            = 'soma_center',
        seclist_name    = 'somatic',
        sec_index       = 0,
        comp_x          = 0.5)




class PhaseResponseProtocol(BpopProtocolWrapper):
    """
    Protocol for measuring phase response curves.

    Implementation using standard BluePyOpt classes: SweepProtocol,
    NrnNetStimStimulus, NrndPoinrProcessLocation.

    Based on example BluePyOpt/examples/expsyn/expsyn.py
    """

    def make_synapse_location(self, proximity):
        """
        Get location in dendritic tree based on proximity (prox/mid/dist)
        """
        dist_um = {
            'prox': 50.0,
            'mid'   : 150.0,
            'dist'  : 300.0,
        }[proximity]

        return ephys.locations.NrnSomaDistanceCompLocation(
                                name            = 'dend_' + proximity,
                                seclist_name    = 'dendritic',
                                soma_distance   = dist_um)

    def make_synapse_mechanism(
            self,
            syn_mech_name,
            syn_comp_loc,
            **syn_mod_params):
        """
        Make an Ephys synapse mechanism, including its parameters and
        PointProcessLocation to link them.
        """
        pp_mech = ephys.mechanisms.NrnMODPointProcessMechanism(
                        name='syn1',
                        suffix=syn_mech_name,
                        locations=[syn_comp_loc])

        pp_loc = ephys.locations.NrnPointProcessLocation(
                        name='pploc_syn1', pprocess_mech=pp_mech)

        pp_params = [
            ephys.parameters.NrnPointProcessParameter(
                param_name=k, value=v, frozen=True,
                locations=[pp_loc], name='syn_' + k)
            for k, v in syn_mod_params.items()
        ]


        return pp_mech, pp_params, pp_loc


    def __init__(
            self,
            stn_model_type=None,
            syn_comp_loc=None,
            syn_pp_mech=None,
            syn_pp_params=None,
            syn_pp_loc=None,
            expected_rate=None,
            **kwargs):
        """
        Initialize all protocol variables for given model type

        @param  stn_model_type : cellpopdata.StnModel
                Type of model to use

        @param  syn_comp_loc : ephys.locations.Location
                Location for synapse used to adminisiter PRC pulses
        """
        # protocol parameters
        syn_gmax = kwargs.get('syn_gmax', 1e-3)
        bias_rate = kwargs.get('bias_rate', 20.0)

        # Current-firing rate curves
        fullmodel_fI_data = [(0.01, 16.36), (0.015, 20.73), (.02, 25.1),
            (.03, 33.5), (0.032, 34.9), (0.035, 37.1), (.05, 47.3), (.1, 76.2)]
        Idata, fdata = zip(*fullmodel_fI_data)
        bias_current = np.interp(bias_rate, fdata, Idata)


        # Parameters for PRC
        cell_T = 1e3 / expected_rate    # (ms)
        prc_sampling_T = 2.0        # (ms) sampling period within ISI

        # Sample each phase X times, and randomize order
        prc_sampling_repeats = 10
        prc_delays = np.arange(0, cell_T*1.2, prc_sampling_T)
        prc_delays = np.tile(prc_delays, prc_sampling_repeats)
        np.random.shuffle(prc_delays) # modifies in-place
        logger.debug('Delays for sampling ISI of {} ms (including repeats): {}'.format(
                        cell_T, sorted(prc_delays)))

        # Global parameters
        stim_start = 200
        stim_stop = stim_start + (1.1 * cell_T * len(prc_delays))
        sim_dur = stim_stop
        if stim_stop > 5000.0:
            logger.debug('Long simulation time ({} ms) for combination of cell '
                         'firing rate, PRC sampling interval, repeats.'.format(
                             sim_dur))

        self.response_interval = (stim_start, sim_dur)


        # Current stimulus -----------------------------------------------------
        # To get neuron firing at target rate

        # STN firing rates in control condition are quite low < 17 Hz, so we
        # choose firing rate in DD condition

        stim_bias = ephys.stimuli.NrnSquarePulse(
                        step_amplitude  = bias_current,
                        step_delay      = 0.0,
                        step_duration   = sim_dur,
                        location        = loc_soma_center,
                        total_duration  = sim_dur)

        # Synaptic stimulus ----------------------------------------------------

        # Synaptic mechanism


        # TODO: set weight dynamically? Or calibrate once based on passive Ztransfer
        # NOTE: mechs and params passed to cellmodel in our code

        # Stimulate using NetVarDelay (adaptive, feeds back spikes with delay)
        stim_syn_prox = bpop_stimuli.NetVarDelayStimulus(
                        delays=prc_delays,
                        start_time=stim_start,
                        target_locations=[syn_pp_loc],
                        source_location=loc_soma_center,
                        source_threshold=-20.0,
                        source_delay=0.1,
                        target_delay=0.1,
                        target_weight=syn_gmax,
                        total_duration=sim_dur)

        rec_stim_syn1 = bpop_recordings.NetStimRecording(
                        name='PRC.stim_times', # name required by feature calculator
                        netstim=stim_syn_prox)


        rec_soma_v = ephys.recordings.CompRecording(
                        name            = '{}.soma.v'.format(self.IMPL_PROTO.name),
                        location        = loc_soma_center,
                        variable        = 'v')

        # Final protocol -------------------------------------------------------

        proto_stimuli = [stim_syn_prox, stim_bias]
        proto_recordings = [rec_soma_v, rec_stim_syn1]

        self.ephys_protocol = PhysioProtocol(
                        name        = self.IMPL_PROTO.name,
                        stimuli     = proto_stimuli,
                        recordings  = proto_recordings,
                        init_func   = init_stn_physiology)

        # TODO: fill in all params
        self.proto_vars = {
            'pp_mechs'      : [syn_pp_mech],
            'pp_comp_locs'  : [loc_soma_center, syn_comp_loc],
            'pp_target_locs': [syn_pp_loc],
            'pp_mech_params': syn_pp_params,
            'stims'         : proto_stimuli,
            'recordings'    : proto_recordings,
            # 'range_mechs' : [],
        }

        # Characterizing features and parameters for protocol
        # NOTE: these are feat_params used in feature_factory/make_features
        self.characterizing_feats = {
            'PRC_traditional': {
                'weight'        : 1.0,  # TODO: PRC weight
                'norm_factor'   : 1.0,  # TODO: PRC norm factor
                'int': {
                    'smooth_method'     : 'polyfit',
                    'distance_metric'   : 'sumsquared',
                    'poly_order'        : 3,
                },
                'traces'        : {
                    ''              : rec_soma_v.name,
                    syn_pp_loc.name : rec_stim_syn1.name,
                },
            },
            # TODO: other features to regularize cost (spike rate, ...)
        }

class PhaseResponseSynExcDist(PhaseResponseProtocol):
    """
    Protocol for measuring PRC with distally located, excitatory stimulus.
    """

    IMPL_PROTO = StimProtocol.PRC_SYN_EXC_DIST

    def __init__(
            self,
            stn_model_type=None,
            **kwargs):

        syn_loc = self.make_synapse_location(proximity='dist')

        pp_mech, pp_params, pp_loc = self.make_synapse_mechanism(
            'Exp2Syn', syn_loc,  e=0.0, tau1=1.0, tau2=3.0)

        return PhaseResponseProtocol.__init__(self,
                    stn_model_type  = stn_model_type,
                    syn_comp_loc    = syn_loc,
                    syn_pp_mech     = pp_mech,
                    syn_pp_params   = pp_params,
                    syn_pp_loc      = pp_loc,
                    syn_gmax        = 1e-3,
                    expected_rate   = 20.0,
                    **kwargs)


class PhaseResponseSynExcProx(PhaseResponseProtocol):
    """
    Protocol for measuring PRC with proximally located, excitatory stimulus.
    """

    IMPL_PROTO = StimProtocol.PRC_SYN_EXC_PROX

    def __init__(
            self,
            stn_model_type=None,
            **kwargs):

        syn_loc = self.make_synapse_location(proximity='prox')

        pp_mech, pp_params, pp_loc = self.make_synapse_mechanism(
            'Exp2Syn', syn_loc,  e=0.0, tau1=1.0, tau2=3.0)

        return PhaseResponseProtocol.__init__(self,
                    stn_model_type  = stn_model_type,
                    syn_comp_loc    = syn_loc,
                    syn_pp_mech     = pp_mech,
                    syn_pp_params   = pp_params,
                    syn_pp_loc      = pp_loc,
                    syn_gmax        = 1e-3,
                    expected_rate   = 20.0,
                    **kwargs)


class PhaseResponseSynInhProx(PhaseResponseProtocol):
    """
    Protocol for measuring PRC with proximally located, inhibitory stimulus.
    """

    IMPL_PROTO = StimProtocol.PRC_SYN_INH_PROX

    def __init__(
            self,
            stn_model_type=None,
            **kwargs):

        syn_loc = self.make_synapse_location(proximity='prox')

        pp_mech, pp_params, pp_loc = self.make_synapse_mechanism(
            'Exp2Syn', syn_loc,  e=-80.0, tau1=1.0, tau2=12.0)

        return PhaseResponseProtocol.__init__(self,
                    stn_model_type  = stn_model_type,
                    syn_loc         = syn_loc,
                    syn_comp_loc    = syn_loc,
                    syn_pp_mech     = pp_mech,
                    syn_pp_params   = pp_params,
                    syn_pp_loc      = pp_loc,
                    syn_gmax        = 1e-3,
                    expected_rate   = 15.0,
                    **kwargs)


# Register protocols implemented here
wrappers = [
    PhaseResponseSynExcDist, PhaseResponseSynExcProx,
    PhaseResponseSynInhProx
]

PROTOCOL_WRAPPERS.update({
    wrapper.IMPL_PROTO: wrapper for wrapper in wrappers
})
