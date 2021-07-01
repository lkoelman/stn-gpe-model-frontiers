"""
Protocols for optimization of reduced neuron model using BluePyOpt.

@author Lucas Koelman

@date   25/09/2017


NOTES

- see SETPARAM comment for parameters that should be set by user

"""

# Standard library
import collections

# Third party
import bluepyopt.ephys as ephys

# Custom modules
from bgcellmodels.extensions.bluepyopt.bpop_stimuli import NrnSpaceClamp
from bgcellmodels.extensions.bluepyopt.bpop_protocol_ext import (
    SelfContainedProtocol, PhysioProtocol, BpopProtocolWrapper,
    rng_getter, connector_getter, PROTOCOL_WRAPPERS
)
from bgcellmodels.models.STN.GilliesWillshaw.evalmodel import (
    proto_background, proto_simple_syn
)
from bgcellmodels.common import stimprotocols
import bgcellmodels.cellpopdata as cpd

StimProtocol = stimprotocols.StimProtocol

################################################################################
# PROTOCOLS
################################################################################

# Location for current clamp
soma_center_loc = ephys.locations.NrnSeclistCompLocation(
                name            = 'soma_center_loc',
                seclist_name    = 'somatic',
                sec_index       = 0,        # index in SectionList
                comp_x          = 0.5)  # x-location in Section

def init_spont(sim, model):
    """
    Initialize simulator to run plateau protocol

    NOTE: function must be declared at top-level of module in order to be pickled
    """
    h = sim.neuron.h
    h.celsius = 35
    h.v_init = -60
    h.set_aCSF(4) # gillies_model.set_aCSF(4) # if called before model.instantiate()
    h.init()


class BpopSpontProtocol(BpopProtocolWrapper):
    """
    Functions for setting up plateau protocol.
    """

    IMPL_PROTO = StimProtocol.SPONTANEOUS

    def __init__(self, stn_model_type=None):
        """
        Initialize all protocol variables for given model type

        @param stn_model_type       cellpopdata.StnModel enum instance
        """

        # stimulus interval for eFEL features
        sim_dur = 2000.0
        self.response_interval = (200.0, sim_dur)

        stim_dummy = ephys.stimuli.NrnSquarePulse(
                        step_amplitude  = 0,
                        step_delay      = 0,
                        step_duration   = 0,
                        location        = soma_center_loc,
                        total_duration  = sim_dur)


        rec_soma_V = ephys.recordings.CompRecording(
                        name            = '{}.soma.v'.format(self.IMPL_PROTO.name),
                        location        = soma_center_loc,
                        variable        = 'v')

        self.ephys_protocol = PhysioProtocol(
                        name        = self.IMPL_PROTO.name, 
                        stimuli     = [stim_dummy],
                        recordings  = [rec_soma_V],
                        init_func   = init_spont)

        self.proto_vars = {
            # 'pp_mechs':           [],
            # 'pp_comp_locs':       [],
            # 'pp_target_locs': [],
            # 'pp_mech_params': [],
            'stims':            [stim_dummy],
            'recordings':       [rec_soma_V],
            # 'range_mechs':        [],
        }

    # Characterizing features and parameters for protocol
    characterizing_feats = {
        ### SPIKE TIMING RELATED ###
        'Spikecount': {         # (int) The number of peaks during stimulus
            'weight':   2.0,
        },
        'ISI_CV': {             # (float) coefficient of variation of ISI durations
            'weight':   1.0,
        },
        ### SPIKE SHAPE RELATED ###
        'AP_rise_rate': {       # (array) Voltage change rate during the rising phase of the action potential
            'weight':   1.0,
        },
        'spike_half_width': {   # (array) The FWHM of each peak
            'weight':   1.0,
        },
        'min_AHP_values': {     # (array) Voltage values at the AHP
            'weight':   1.0,    # this feature recognizes that there is an elevated plateau
        },
    }

################################################################################
# PLATEAU protocol

def init_plateau(sim, model):
    """
    Initialize simulator to run plateau protocol

    NOTE: function must be declared at top-level of module in order to be pickled
    """
    h = sim.neuron.h
    h.celsius = 30
    h.v_init = -60
    h.set_aCSF(4) # gillies_model.set_aCSF(4) # if called before model.instantiate()
    h.init()


class BpopPlateauProtocol(BpopProtocolWrapper):
    """
    Functions for setting up plateau protocol.
    """

    IMPL_PROTO = StimProtocol.CLAMP_PLATEAU

    def __init__(self, stn_model_type=None):
        """
        Initialize all protocol variables for given model type

        @param stn_model_type       cellpopdata.StnModel enum instance

        @post                       following attributes will be available on this object:
                                    - ephys_protocol: PhysioProtocol instance
                                    - proto_vars: dict with protocol variables
                                    - response_interval: expected time interval of response
        """
        # stimulus parameters
        I_hyper = -0.17         # hyperpolarize to -70 mV (see fig. 10C)
        I_depol = I_hyper + 0.2 # see fig. 10D: 0.2 nA (=stim.amp) over hyperpolarizing current

        del_depol = 1000
        dur_depol = 50          # see fig. 10D, top right
        dur_total = 2000

        # stimulus interval for eFEL features
        plat_start = del_depol - 50
        plat_stop = del_depol + 200
        self.response_interval = (plat_start, plat_stop)

        stim1_hyp = ephys.stimuli.NrnSquarePulse(
                        step_amplitude  = I_hyper,
                        step_delay      = 0,
                        step_duration   = del_depol,
                        location        = soma_center_loc,
                        total_duration  = del_depol)

        stim2_dep = ephys.stimuli.NrnSquarePulse(
                        step_amplitude  = I_depol,
                        step_delay      = del_depol,
                        step_duration   = dur_depol,
                        location        = soma_center_loc,
                        total_duration  = del_depol + dur_depol)

        stim3_hyp = ephys.stimuli.NrnSquarePulse(
                        step_amplitude  = I_hyper,
                        step_delay      = del_depol + dur_depol,
                        step_duration   = del_depol,
                        location        = soma_center_loc,
                        total_duration  = dur_total)

        plat_rec1 = ephys.recordings.CompRecording(
                        name            = '{}.soma.v'.format(self.IMPL_PROTO.name),
                        location        = soma_center_loc,
                        variable        = 'v')

        self.ephys_protocol = PhysioProtocol(
                        name        = self.IMPL_PROTO.name, 
                        stimuli     = [stim1_hyp, stim2_dep, stim3_hyp],
                        recordings  = [plat_rec1],
                        init_func   = init_plateau)

        self.proto_vars = {
            # 'pp_mechs':           [],
            # 'pp_comp_locs':       [],
            # 'pp_target_locs': [],
            # 'pp_mech_params': [],
            'stims':            [stim1_hyp, stim2_dep, stim3_hyp],
            'recordings':       [plat_rec1],
            # 'range_mechs':        [],
        }

    # Characterizing features and parameters for protocol
    characterizing_feats = {
        ### SPIKE TIMING RELATED ###
        # 'Spikecount': {           # (int) The number of peaks during stimulus
        #   'weight':   2.0,
        # },
        # 'mean_frequency',     # (float) the mean frequency of the firing rate
        # 'burst_mean_freq',    # (array) The mean frequency during a burst for each burst
        'adaptation_index': {   # (float) Normalized average difference of two consecutive ISIs
            'weight':   1.0,
            'double':   {'spike_skipf': 0.0},
            'int':      {'max_spike_skip': 0},
        },
        'ISI_CV': {             # (float) coefficient of variation of ISI durations
            'weight':   1.0,
        },
        # 'ISI_log',            # no documentation
        'Victor_Purpura_distance': {
            'weight':   2.0,
            'exp_std': 50.0,
            'double': { 'spike_shift_cost_ms' : 20.0 }, # 20 ms is kernel quarter width
        },
        ### SPIKE SHAPE RELATED ###
        # 'AP_duration_change': {   # (array) Difference of the durations of the second and the first action potential divided by the duration of the first action potential
        #   'weight':   1.0,
        # },
        # 'AP_duration_half_width_change': { # (array) Difference of the FWHM of the second and the first action potential divided by the FWHM of the first action potential
        #   'weight':   1.0,
        # },
        # 'AP_rise_time': {     # (array) Time from action potential onset to the maximum
        #   'weight':   1.0,
        # },
        'AP_rise_rate': {       # (array) Voltage change rate during the rising phase of the action potential
            'weight':   1.0,
        },
        'AP_height': {          # (array) The voltages at the maxima of the peak
            'weight':   1.0,
        },
        # 'AP_amplitude': {     # (array) The relative height of the action potential
        #   'weight':   1.0,
        # },
        'spike_half_width': {   # (array) The FWHM of each peak
            'weight':   1.0,
        },
        # 'AHP_time_from_peak': {   # (array) Time between AP peaks and AHP depths
        #   'weight':   1.0,
        # },
        # 'AHP_depth': {            # (array) relative voltage values at the AHP
        #   'weight':   1.0,
        # },
        'min_AHP_values': {     # (array) Voltage values at the AHP
            'weight':   2.0,    # this feature recognizes that there is an elevated plateau
        },
    }

################################################################################
# REBOUND protocol

def init_rebound(sim, model):
    """
    Initialize simulator to run rebound protocol

    NOTE: function must be declared at top-level of module in order to be pickled
    """
    h = sim.neuron.h
    h.celsius = 35
    h.v_init = -60
    h.set_aCSF(4)
    h.init()


class BpopReboundProtocol(BpopProtocolWrapper):
    """
    Functions for setting up rebound protocol.
    """

    IMPL_PROTO = StimProtocol.CLAMP_REBOUND

    def __init__(self, **kwargs):
        """
        Initialize all protocol variables for given model type
        """

        dur_hyper = 500
        dur_burst = 150 # SETPARAM: approximate duration of burst
        reb_start, reb_stop = dur_hyper-50, dur_hyper+dur_burst # SETPARAM: stimulus interval for eFEL features
        self.response_interval = (reb_start, reb_stop)

        reb_clmp1 = NrnSpaceClamp(
                        step_amplitudes = [0, 0, -75],
                        step_durations  = [0, 0, dur_hyper],
                        total_duration  = reb_stop,
                        location        = soma_center_loc)

        reb_rec1 = ephys.recordings.CompRecording(
                        name            = '{}.soma.v'.format(self.IMPL_PROTO.name),
                        location        = soma_center_loc,
                        variable        = 'v')

        self.ephys_protocol = PhysioProtocol(
                        name        = self.IMPL_PROTO.name, 
                        stimuli     = [reb_clmp1],
                        recordings  = [reb_rec1],
                        init_func   = init_rebound)

        self.proto_vars = {
            # 'pp_mechs':           [],
            # 'pp_comp_locs':       [],
            # 'pp_target_locs': [],
            # 'pp_mech_params': [],
            'stims':            [reb_clmp1],
            'recordings':       [reb_rec1],
            # 'range_mechs':        [],
        }

    # Characterizing features and parameters for protocol
    characterizing_feats = {
        ### SPIKE TIMING RELATED ###
        'Spikecount': {         # (int) The number of peaks during stimulus
            'weight':   1.0,
        },
        # 'mean_frequency': {   # (float) the mean frequency of the firing rate
        #   'weight':   2.0,
        # },
        # 'burst_mean_freq',    # (array) The mean frequency during a burst for each burst
        'adaptation_index': {   # (float) Normalized average difference of two consecutive ISIs
            'weight':   1.0,
            'double':   {'spike_skipf': 0.0},
            'int':      {'max_spike_skip': 0},
        },
        'ISI_CV': {             # (float) coefficient of variation of ISI durations
            'weight':   1.0,
        },
        # 'ISI_log',            # no documentation
        # 'Victor_Purpura_distance': {
        #   'weight':   2.0,
        #   'exp_std': 50.0,
        #   'double': { 'spike_shift_cost_ms' : 20.0 }, # 20 ms is kernel quarter width
        # },
        ### SPIKE SHAPE RELATED ###
        # 'AP_duration_change': {   # (array) Difference of the durations of the second and the first action potential divided by the duration of the first action potential
        #   'weight':   1.0,
        # },
        # 'AP_duration_half_width_change': { # (array) Difference of the FWHM of the second and the first action potential divided by the FWHM of the first action potential
        #   'weight':   1.0,
        # },
        # 'AP_rise_time': {     # (array) Time from action potential onset to the maximum
        #   'weight':   1.0,
        # },
        'AP_rise_rate': {       # (array) Voltage change rate during the rising phase of the action potential
            'weight':   1.0,
        },
        'AP_height': {          # (array) The voltages at the maxima of the peak
            'weight':   2.0,
        },
        'AP_amplitude': {       # (array) The relative height of the action potential
            'weight':   2.0,
        },
        'spike_half_width': {   # (array) The FWHM of each peak
            'weight':   1.0,
        },
        # 'AHP_time_from_peak': {   # (array) Time between AP peaks and AHP depths
        #   'weight':   1.0,
        # },
        # 'AHP_depth': {            # (array) relative voltage values at the AHP
        #   'weight':   1.0,
        # },
        'min_AHP_values': {     # (array) Voltage values at the AHP
            'weight':   3.0,    # this feature recognizes that there is an elevated plateau
        },
    }

################################################################################
# SYNAPTIC protocol

class BpopSynBurstProtocol(BpopProtocolWrapper):
    """
    Functions for setting up rebound protocol.
    """

    IMPL_PROTO = StimProtocol.MIN_SYN_BURST

    def __init__(self, **kwargs):
        """
        Initialize all protocol variables for given model type

        @post                       following attributes will be available on this object:
                                    - ephys_protocol: PhysioProtocol instance
                                    - proto_vars: dict with protocol variables
                                    - response_interval: expected time interval of response
        """

        transient_delay = 350.0
        stim_delay = 700.0
        sim_end = stim_delay + 500.0
        self.spont_interval = (transient_delay, stim_delay)
        self.burst_interval = (stim_delay+5.0, 1010.0)
        self.total_interval = (transient_delay, sim_end)
        self.response_interval = self.total_interval

        rec_soma_V = ephys.recordings.CompRecording(
                        name            = '{}.soma.v'.format(self.IMPL_PROTO.name),
                        location        = soma_center_loc,
                        variable        = 'v')

        # Protocol setup
        proto_setup_funcs = [
            proto_simple_syn.make_inputs_BURST_impl
        ]

        proto_init_funcs = [
            proto_simple_syn.init_sim_BURST_impl
        ]

        proto_setup_kwargs_const = {
            'do_map_synapses': True,
            'base_seed':  25031989, # SETPARAM: same as StnModelEvaluator
            'gid': 1, # same as StnModelEvaluator
            'physio_state': cpd.PhysioState.NORMAL.name,
            'n_gpe_syn': 1,
            'n_ctx_syn': 4,
            'delay': stim_delay,
        }

        proto_setup_kwargs_getters = collections.OrderedDict([
            ('rng', rng_getter),
            ('connector', connector_getter),
        ])

        # Recording and plotting traces
        proto_rec_funcs = [
            # proto_background.rec_spikes,
        ]

        proto_plot_funcs = [
            # stimprotocols.plot_all_spikes,
            # stimprotocols.report_spikes,
        ]

        self.ephys_protocol = SelfContainedProtocol(
                    name                = self.IMPL_PROTO.name, 
                    recordings          = [rec_soma_V],
                    total_duration      = sim_end,
                    proto_init_funcs            = proto_init_funcs,
                    proto_setup_funcs_pre       = proto_setup_funcs,
                    proto_setup_kwargs_const    = proto_setup_kwargs_const,
                    proto_setup_kwargs_getters  = proto_setup_kwargs_getters,
                    rec_traces_funcs    = proto_rec_funcs,
                    plot_traces_funcs   = proto_plot_funcs)

        self.proto_vars = {
            'pp_mechs':         [],
            'pp_comp_locs':     [],
            'pp_target_locs':   [],
            'pp_mech_params':   [],
            'stims':            [],
            'recordings':       [rec_soma_V],
            'range_mechs':      [],
        }

        # Characterizing features and parameters for protocol
        self.characterizing_feats = {
            ### SPIKE TIMING RELATED ###
            'Spikecount': {         # (int) The number of peaks during stimulus
                'weight':   4.0,
                'response_interval': self.burst_interval,
            },
            'mean_frequency': { # (float) the mean frequency of the firing rate
                'weight':   4.0,
                'response_interval': self.spont_interval, # SPONTANEOUS period
            },
            # 'burst_mean_freq',    # (array) The mean frequency during a burst for each burst
            'adaptation_index': {   # (float) Normalized average difference of two consecutive ISIs
                'weight':   1.0,
                'double':   {'spike_skipf': 0.0},
                'int':      {'max_spike_skip': 0},
                'response_interval': self.burst_interval,
            },
            'ISI_CV': {             # (float) coefficient of variation of ISI durations
                'weight':   1.0,
                'response_interval': self.burst_interval,
            },
            'time_to_first_spike': {
                'weight':   2.0,
                'response_interval': self.burst_interval,
            },
            # 'ISI_log',            # no documentation
            'Victor_Purpura_distance': {
                'weight':   2.0,
                'exp_std': 50.0,
                'double': { 'spike_shift_cost_ms' : 20.0 }, # 20 ms is kernel quarter width
                'response_interval': self.total_interval,
            },
            ### SPIKE SHAPE RELATED ###
            # 'AP_duration_change': {   # (array) Difference of the durations of the second and the first action potential divided by the duration of the first action potential
            #   'weight':   1.0,
            # },
            # 'AP_duration_half_width_change': { # (array) Difference of the FWHM of the second and the first action potential divided by the FWHM of the first action potential
            #   'weight':   1.0,
            # },
            # 'AP_rise_time': {     # (array) Time from action potential onset to the maximum
            #   'weight':   1.0,
            # },
            'AP_rise_rate': {       # (array) Voltage change rate during the rising phase of the action potential
                'weight':   1.0,
                'response_interval': self.spont_interval,
            },
            'AP_height': {          # (array) The voltages at the maxima of the peak
                'weight':   2.0,
                'response_interval': self.burst_interval,
            },
            'AP_amplitude': {       # (array) The relative height of the action potential
                'weight':   2.0,
                'response_interval': self.burst_interval,
            },
            'spike_half_width': {   # (array) The FWHM of each peak
                'weight':   1.0,
                'response_interval': self.spont_interval,
            },
            # 'AHP_time_from_peak': {   # (array) Time between AP peaks and AHP depths
            #   'weight':   1.0,
            # },
            # 'AHP_depth': {            # (array) relative voltage values at the AHP
            #   'weight':   1.0,
            # },
            'min_AHP_values': {     # (array) Voltage values at the AHP
                'weight':   2.0,    # this feature recognizes that there is an elevated plateau
                'response_interval': self.burst_interval,
            },
        }

################################################################################
# BACKGROUND protocol


class BpopBackgroundProtocol(BpopProtocolWrapper):
    """
    Functions for setting up rebound protocol.
    """

    IMPL_PROTO = StimProtocol.SYN_BACKGROUND_HIGH

    def __init__(
            self,
            impl_proto=None,
            **kwargs):
        """
        Initialize all protocol variables for given model type

        @param    kwargs            any other keyword arguments will be passed to 
                                    SelfContainedProtocol as 'proto_setup_kwargs_const'

        @post                       following attributes will be available on this object:
                                    - ephys_protocol: PhysioProtocol instance
                                    - proto_vars: dict with protocol variables
                                    - response_interval: expected time interval of response
        """
        if impl_proto is not None:
            self.IMPL_PROTO = impl_proto # quick fix

        sim_end = 3000.0
        self.response_interval = (300.0, sim_end)

        bg_recV = ephys.recordings.CompRecording(
                        name            = '{}.soma.v'.format(self.IMPL_PROTO.name),
                        location        = soma_center_loc,
                        variable        = 'v')

        # Protocol setup
        proto_setup_funcs = [
            proto_background.make_inputs_ctx_impl, 
            proto_background.make_inputs_gpe_impl
        ]

        proto_init_funcs = [
            proto_background.init_sim_impl
        ]

        proto_setup_kwargs_const = {
            'base_seed': 8, # SETPARAM: 25031989 in StnModelEvaluator,
            'gid': 1,
            'do_map_synapses': True,
            'physio_state': cpd.PhysioState.NORMAL.name,
        }
        # Pass user arguments
        proto_setup_kwargs_const.update(kwargs)

        proto_setup_kwargs_getters = collections.OrderedDict([
            ('rng', rng_getter),
            ('connector', connector_getter),
        ])

        # Recording and plotting traces
        proto_rec_funcs = [
            proto_background.rec_spikes,
        ]

        proto_plot_funcs = [
            stimprotocols.plot_all_spikes,
            stimprotocols.report_spikes,
        ]

        self.ephys_protocol = SelfContainedProtocol(
                            name                = self.IMPL_PROTO.name, 
                            recordings          = [bg_recV],
                            total_duration      = sim_end,
                            proto_init_funcs            = proto_init_funcs,
                            proto_setup_funcs_pre       = proto_setup_funcs,
                            proto_setup_kwargs_const    = proto_setup_kwargs_const,
                            proto_setup_kwargs_getters  = proto_setup_kwargs_getters,
                            rec_traces_funcs    = proto_rec_funcs,
                            plot_traces_funcs   = proto_plot_funcs)

        self.proto_vars = {
            'pp_mechs':         [],
            'pp_comp_locs':     [],
            'pp_target_locs':   [],
            'pp_mech_params':   [],
            'stims':            [],
            'recordings':       [bg_recV],
            'range_mechs':      [],
        }

    # Characterizing features and parameters for protocol
    characterizing_feats = {
        'Victor_Purpura_distance': {
            'weight':   8.0,
            'exp_std': 50.0,
            'double': { 'spike_shift_cost_ms' : 20.0 }, # 20 ms is kernel quarter width
        },
        # 'burst_mean_freq': {
        #   'weight':   2.0,
        # }
        'ISI_voltage_distance': {
            'weight':   5.0,
            'exp_std': 150.0,
        },
        'instantaneous_rate': {
            'weight':   5.0,
            'exp_std': 700.0,
            'int':      {'min_AP': 2.0},
            'double':   {'bin_width': 50.0},
        },
    }


# Register protocols implemented here
wrappers = [
    BpopSpontProtocol, BpopPlateauProtocol, BpopReboundProtocol,
    BpopSynBurstProtocol, BpopBackgroundProtocol
]

PROTOCOL_WRAPPERS.update({
    wrapper.IMPL_PROTO: wrapper for wrapper in wrappers
})

# One class used for two protocols
PROTOCOL_WRAPPERS[StimProtocol.SYN_BACKGROUND_LOW] = BpopBackgroundProtocol
