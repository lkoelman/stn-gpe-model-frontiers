"""
Custom EFeatures to compare spike trains using different similarity mesures / 
distance metrics.

@author Lucas Koelman

@date   4/10/2017


NOTES

- the returned score of efeature.calculate_score() is the distance to the target
  normalized by the standard devation
- i.e. distance is sum_{i..N}(feat[i] - exp_mean) / N / exp_std
- see https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/ephys/efeatures.py
- this meands the standard deviation is the inverse 'weight' of the score
    - weight = 1 / efeat.exp_std



ARCHITECTURE:

   - make custom =EFeature= class that allows to set a target (in the form of spike train or histogram, e.g. a neo spiketrain)

       - this can use eFEL to get the spike times using =efel.getFeatureValues('peak_time')=

       - then call *elephant* or other library to compute distance with target spike train

       - as long as this offers =calculate_score(responses)=, you don't need a custom Objective class

   
   - *PROS* flexible, easy to integrate and easily call into any module

"""

# Third party libraries
import numpy as np

import bluepyopt.ephys as ephys
from elephant import spike_train_dissimilarity as stds
import quantities as pq
import neo

# Custom modules
from bgcellmodels.common import signal, curvedist

import logging
logger = logging.getLogger(__name__)

DEBUG_TESTRUN = True
debug_data = None


def getFeatureNames():
    """
    Return list with names of available features.
    """
    return SpikeTrainFeature.CALC_SCORE_FUNCS.keys()


################################################################################
# Feature input value calculation
################################################################################


def calc_feat_values_spiketimes(self, efel_trace, raise_warnings):
    """
    Calculate feature values for spike-timing based metrics.
    """

    self._setup_efel()
    import efel

    efel_feats = ['peak_time']

    # Calculate spike times from response
    values = efel.getFeatureValues(
        [efel_trace],
        efel_feats,
        raise_warnings=raise_warnings
    )
    efeat_values = {
        feat_name: values[0][feat_name] for feat_name in efel_feats
    }

    efel.reset()

    return efeat_values


def calc_feat_values_ISI_voltage(self, efel_trace, raise_warnings):
    """
    Calculate feature values required for ISI voltage distance.
    """

    self._setup_efel()
    import efel

    # Calculate required features / dependencies
    efel_feats = ['AP_begin_indices', 'AP_end_indices']
    values = efel.getFeatureValues(
        [efel_trace],
        efel_feats,
        raise_warnings=raise_warnings
    )
    efeat_values = {
        feat_name: values[0][feat_name] for feat_name in efel_feats
    }

    # Voltage itself is required as well
    # TODO: only save values where stim_start <= T <= stim_end
    efeat_values['V'] = efel_trace['V'].values
    efeat_values['dt'] = efel_trace['T'][1] - efel_trace['T'][0]

    efel.reset()

    return efeat_values


################################################################################
# Feature scores calculation
################################################################################

def calc_data_PRC(self, efel_trace, trace_check):
    """
    Calculate feature score for new response with regard to target values

    @pre    efel_trace must contain recording named 'PRC.stim_times'

    @pre    Algorithm parameters set in 'double_settings':
            'smooth_kernel_sigma'

    @see    signal.compute_PRC
    """

    # Get spike times
    self._setup_efel()
    import efel

    feat_vals = efel.getFeatureValues([efel_trace], ['peak_time'],
                                      raise_warnings = True)
    
    spike_times = feat_vals[0]['peak_time']

    # stimulus times (are just recordings)
    stim_times = efel_trace['PRC.stim_times']

    # Phase response curve
    phi0, delta_phi0 = signal.compute_PRC(spike_times, stim_times, sort_by='phi')

    # Bin phi values
    bin_dphi = self.double_settings.get('bin_dphi', 0.025)
    phi, delta_phi = curvedist.bin_curve_samples(phi0, delta_phi0, dx=bin_dphi,
                                                 bin_y_func='median')


    if DEBUG_TESTRUN:
        import matplotlib.pyplot as plt
        # from scipy.ndimage.filters import gaussian_filter
        # sigma_phi = self.double_settings['smooth_kernel_sigma']
        coeff = np.polyfit(phi, delta_phi, 3)
        poly = np.poly1d(coeff)
        dphi_polyfit = [poly(x) for x in phi]
        dphi_smooth = curvedist.smooth(delta_phi, window_len=11, window='hanning')

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(phi0, delta_phi0, '.', label='original')
        ax.plot(phi, dphi_smooth, '-', marker='.', label='smooth')
        ax.plot(phi, dphi_polyfit, label='polynomial')
        ax.legend()

        global debug_data
        debug_data = fig, ax # modifiable
        

    # Save feature data
    efeat_values = {
        'spike_times'   : spike_times,
        'stim_times'    : stim_times,
        'prc_phi'       : phi,
        'prc_delta_phi' : delta_phi,
    }

    return efeat_values


def calc_score_PRC(self, efel_trace, trace_check):
    """
    Calculate feature score for new response with regard to target values

    @see    signal.compute_PRC
    """

    # Get spike and stimulus times
    self._setup_efel()
    import efel
    feat_vals = efel.getFeatureValues([efel_trace], ['peak_time'],
                                      raise_warnings = True)
    spike_times = feat_vals[0]['peak_time']
    stim_times = efel_trace['PRC.stim_times']

    # Calculate PRC for candidate response
    phi_cand, dphi_cand = signal.compute_PRC(spike_times, stim_times, sort_by='phi')

    # Get PRC for original model/experiments
    phi_targ = self.target_value_data['prc_phi']
    dphi_targ = self.target_value_data['prc_delta_phi']

    # Compare phase response curves
    # NOTE: by choosing the target PRC as reference, you can ensure it has
    #       enough sample points. If the candidate PRC would be the reference,
    #       you could get very small distance for unequally sampled PRC.
    dist_params = dict(self.int_settings)
    dist_params.update(self.double_settings)
    dist_params['reference'] = 0

    score = curvedist.curve_distance((phi_targ, dphi_targ),
                                     (phi_cand, dphi_cand), **dist_params)

    score /= self.exp_std
    if self.force_max_score:
            score = max(score, self.max_score)

    return score


def calc_score_VP_dist(self, efel_trace, trace_check):
    """
    Calculate feature score for new response with regard to target values
    """

    self._setup_efel()
    import efel

    # Calculate spike times from response
    efel_feat = 'peak_time'
    feat_vals = efel.getFeatureValues(
        [efel_trace],
        [efel_feat],
        raise_warnings = True
    )
    resp_spike_times = feat_vals[0][efel_feat]
    resp_spike_train = self._construct_neo_spiketrain(resp_spike_times)
    target_spike_times = self.target_value_data['peak_time']
    target_spike_train = self._construct_neo_spiketrain(target_spike_times)


    # Set spike shift cost parameter
    if 'spike_shift_cost_ms' in self.double_settings:
        cost_ms = self.double_settings['spike_shift_cost_ms']
        spike_shift_cost = 1.0/(cost_ms*1e-3) * pq.Hz
    
    elif 'spike_shift_cost_hz' in self.double_settings:
        spike_shift_cost = self.double_settings['spike_shift_cost_hz'] * pq.Hz
    
    else:
        spike_shift_cost = 1.0/(20e-3) * pq.Hz # 20 ms is kernel quarter width

    # Compute distance function
    dist_mat = stds.victor_purpura_dist(
                    [target_spike_train, resp_spike_train], 
                    q = spike_shift_cost, 
                    algorithm = 'fast')

    score = dist_mat[0, 1]
    

    # Need to take into account std (can be used as weight too)
    efel.reset()
    score /= self.exp_std
    if self.force_max_score:
            score = max(score, self.max_score)

    return score


def calc_score_instantaneous_rate(self, efel_trace, trace_check):
    """
    Calculate feature score for new response with regard to target values
    """

    self._setup_efel()
    import efel
    from bgcellmodels.common.analysis import numpy_avg_rate_simple
    import numpy as np

    # Calculate spike times from response
    efel_feat = 'peak_time'
    feat_vals = efel.getFeatureValues(
        [efel_trace],
        [efel_feat],
        raise_warnings = True
    )
    resp_spike_times = feat_vals[0][efel_feat]
    target_spike_times = self.target_value_data['peak_time']

    # min_spk = self.int_settings.get('min_AP', 2)
    bin_width = self.double_settings.get('bin_width', 50.0)

    resp_psth = numpy_avg_rate_simple(
                        [resp_spike_times], 
                        self.stim_start, self.stim_end,
                        bin_width)

    tar_psth  = numpy_avg_rate_simple(
                        [target_spike_times], 
                        self.stim_start, self.stim_end,
                        bin_width)

    # Sum of squared differences, averaged
    score = np.sum((tar_psth-resp_psth)**2) / tar_psth.size

    # Need to take into account std (can be used as weight too)
    efel.reset()
    score /= self.exp_std
    if self.force_max_score:
            score = max(score, self.max_score)

    return score


def calc_score_Kreuz_ISI_dist(self, efel_trace, trace_check):
    """
    Calculate feature score for new response with regard to target values
    """

    self._setup_efel()
    import efel
    import pyspike

    # Calculate spike times from response
    efel_feat = 'peak_time'
    feat_vals = efel.getFeatureValues(
        [efel_trace],
        [efel_feat],
        raise_warnings = True
    )

    # Construct spike trains
    resp_spike_times = feat_vals[0][efel_feat]
    st_resp = pyspike.SpikeTrain(
                        resp_spike_times,
                        [self.stim_start, self.stim_end],
                        is_sorted=True)

    target_spike_times = self.target_value_data['peak_time']
    st_targ = pyspike.SpikeTrain(
                        target_spike_times,
                        [self.stim_start, self.stim_end],
                        is_sorted=True)

    # Sum of squared differences, averaged
    score = pyspike.isi_distance(st_targ, st_resp)

    # Need to take into account std (can be used as weight too)
    efel.reset()
    score /= self.exp_std
    if self.force_max_score:
            score = max(score, self.max_score)

    return score


def calc_score_ISI_voltage(self, efel_trace, trace_check):
    """
    Calculate feature score for new response with regard to target values
    """

    self._setup_efel()
    import efel
    import numpy as np

    # import pyximport; pyximport.install()
    # from efeatures_fast_cython import calc_ISI_voltage_distance_dt_equal
    from efeatures_fast_numba import calc_ISI_voltage_distance_dt_equal

    # Calculate required features / dependencies
    efel_feats = ['AP_begin_indices', 'AP_end_indices']
    feat_vals = efel.getFeatureValues(
        [efel_trace],
        efel_feats,
        raise_warnings=True
    )

    # Compute distance function
    tar_AP_begin    = self.target_value_data['AP_begin_indices']
    tar_AP_end      = self.target_value_data['AP_end_indices']
    tar_Vm          = self.target_value_data['V']
    tar_dt          = self.target_value_data['dt']

    cur_AP_begin    = feat_vals[0]['AP_begin_indices']
    cur_AP_end      = feat_vals[0]['AP_end_indices']
    cur_Vm          = efel_trace['V'].values # pandas.Series to numpy.ndarray
    cur_dt          = efel_trace['T'][1] - efel_trace['T'][0]

    dt_equal = abs(tar_dt-cur_dt) <= 0.00001
    if not dt_equal:
        raise Exception("ISI voltage distance only implemented for traces calculated with equal time step (dt_old={}, dt_new={}).".format(tar_dt, cur_dt))

    for vec in (tar_AP_begin, tar_AP_end, cur_AP_begin, cur_AP_end):
        if not isinstance(vec, np.ndarray):
            return self.max_score
        if not np.issubdtype(vec.dtype, int):
            return self.max_score
    
    score = calc_ISI_voltage_distance_dt_equal(
                            tar_Vm, cur_Vm, 
                            tar_AP_begin, cur_AP_begin,
                            tar_AP_end, cur_AP_end,
                            self.stim_start, self.stim_end, tar_dt)

    # Need to take into account std (can be used as weight too)
    efel.reset()
    score /= self.exp_std
    if self.force_max_score:
        score = max(score, self.max_score)

    return score


################################################################################
# BluePyOpt Ephys classes
################################################################################


class SpikeTrainFeature(ephys.efeatures.EFeature, ephys.serializer.DictMixin):
    """
    Feature representing spike times in a particular interval. The distance is
    computed usinng a given similarity metric.
    """

    # Fields to be serialized for data transfer between threads via pickle/unpickle
    SERIALIZED_FIELDS = ('name', 'metric_name', 'recording_names',
                         'stim_start', 'stim_end', 'target_value_data',
                         'threshold', 'comment')

    # Functions for calculating feature input data from traces
    CALC_FEAT_FUNCS = {
        'Victor_Purpura_distance'   : calc_feat_values_spiketimes,
        'instantaneous_rate'        : calc_feat_values_spiketimes,
        'ISI_voltage_distance'      : calc_feat_values_ISI_voltage,
        'Kreuz_ISI_distance'        : calc_feat_values_spiketimes,
        'PRC_traditional'           : calc_data_PRC, 
    }

    # Functions for calculating score (comparing feature values)
    CALC_SCORE_FUNCS = {
        'Victor_Purpura_distance'   : calc_score_VP_dist,
        'instantaneous_rate'        : calc_score_instantaneous_rate,
        'ISI_voltage_distance'      : calc_score_ISI_voltage,
        'Kreuz_ISI_distance'        : calc_score_Kreuz_ISI_dist,
        'PRC_traditional'           : calc_score_PRC, 
    }

    def __init__(
            self,
            name,
            metric_name=None,
            recording_names=None,
            stim_start=None,
            stim_end=None,
            threshold=None,
            interp_step=None,
            comment='',
            double_settings=None,
            int_settings=None,
            force_max_score=False,
            max_score=250
    ):
        """
        Constructor

        Args:

            name (str):             name of the SpikeTrainFeature object
            
            metric_name (str):      name of similarity measure / distance metric to use
            
            recording_names (dict): eFEL features can accept several recordings
                                    as input. These are pairs 'location_name' :
                                    'recording_name'.
            
            stim_start (float):     stimulation start time (ms)
            
            stim_end (float):       stimulation end time (ms)
            
            threshold(float):       spike detection threshold (mV)
            
            comment (str):          comment
            
            interp_step(float):     interpolation step (ms)
            


            target_trace:           dict {
                                        T:np.array , V:np.array , 
                                        stim_start:[float] , stim_end:[float]
                                    }
        """

        super(SpikeTrainFeature, self).__init__(name, comment)

        if double_settings is None:
            double_settings = {}
        if int_settings is None:
            int_settings = {}

        self.recording_names = recording_names
        self.metric_name = metric_name # function of var 'efeal_feature_name' in original class

        self.stim_start = stim_start
        self.stim_end = stim_end
        
        self.threshold = threshold
        self.interp_step = interp_step
        self.double_settings = double_settings
        self.int_settings = int_settings
        
        self.force_max_score = force_max_score
        self.max_score = max_score
        
        self.exp_std = 1.0
        self.exp_mean = None

        self.target_value_data = {} # data for computing distance


    def set_target_values(self, value_dict):
        """
        Set data needed for computing distance.

        @param  value_dict : dict[str, object]
                result of a call to self.calculate_feature()
        """
        self.target_value_data = value_dict


    def _construct_neo_spiketrain(self, spike_times):
        """
        Construct Neo spike train for use with Elephant library.
        """
        if len(spike_times)>0 and (spike_times[0] < self.stim_start or spike_times[-1] > self.stim_end):
            spike_times = [t for t in spike_times if 
                                        (t >= self.stim_start and t <= self.stim_end)]
        
        spike_train = neo.SpikeTrain(
                            spike_times,
                            units ='ms',
                            t_start = self.stim_start,
                            t_stop = self.stim_end)

        return spike_train


    def _construct_efel_trace(self, responses):
        """
        Construct trace that can be passed to eFEL

        @note   data type of time and voltage series is pandas.core.series.Series

        DEVNOTES
        --------

        The original function with very convoluted logic is in
        https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/ephys/efeatures.py
        
        The minimum data required for all features are "T", "V", 
        "stim_start", "stim_end" and these are set using efel.setDoubleSetting(...)
        internally.

        So the required steps for a Feature (i.e. this class) to use eFel
        is to construct a dict using these four values and pass it to eFel.
        How it chooses one of the recordings to dump into the "T" and "V"
        entries is up to this class.

        The way it's done in the BluePyOpt example is that eFEL features
        (ephys.efeatures.eFELFeature) are always passed the argument
        `recording_names = {'': 'soma.v'}`, and this recording (marked
        with empty location '') is used for setting 'T' and 'V'.
        """

        trace = {}
        if '' not in self.recording_names:
            raise Exception('Recording with location \'\' needs to be in recording_names.\n'
                            'This is the default recording used for T,V.')

        trace['stim_start'] = [self.stim_start]
        trace['stim_end']   = [self.stim_end]

        # Use recording with location '' for T and V
        default_rec_name = self.recording_names['']
        if responses[default_rec_name] is None:
            return None # same as original behavior

        trace['T']  = responses[default_rec_name]['time']
        trace['V'] = responses[default_rec_name]['voltage']
        
        for location_name, recording_name in self.recording_names.items():
            if location_name == '':
                continue

            if recording_name not in responses:
                raise ValueError('Expected recording "{}" in responses.'.format(
                                    recording_name))

            if responses[recording_name] is None:
                return None

            if recording_name.endswith('.v'):
                # assume it is CompRecording -> original behavior
                trace['V;' + location_name] = responses[recording_name]['voltage']
                trace['T;' + location_name] = responses[recording_name]['time']
            
            elif recording_name.endswith('_times'):
                # assume it is spike times, e.g. NetStimRecording
                # TODO: for now, NetStimResponse is saved as TimeVoltageResponse, but could use Neo SpikeTrain
                trace[recording_name] = responses[recording_name]['time']
            

        return trace

    def _setup_efel(self):
        """
        Set up efel before extracting the feature
        """

        import efel
        efel.reset()

        if self.threshold is not None:
            efel.setThreshold(self.threshold)

        if self.interp_step is not None:
            efel.setDoubleSetting('interp_step', self.interp_step)


    def calculate_feature(self, responses, raise_warnings=True):
        """
        Calculate feature value for comparison with a candidate model. 
        The result of the comparison is the feature score.

        The resulting feature data is assigned as a target by calling
        set_target_value(), which assigns self.target_value_data that is
        is used in the score functions.

        @return     feat_vals : dict
                    All feature data from the target model that is required
                    to calculate a feature score for a candidate model.
        """

        efel_trace = self._construct_efel_trace(responses)

        if efel_trace is None:
            feat_vals = None
        else:
            feat_name = self.metric_name
            feat_func = self.CALC_FEAT_FUNCS[feat_name]
            feat_vals = feat_func(self, efel_trace, raise_warnings)

        # logger.debug('Calculated feature value for %s: %s', self.name, feat_vals)
        
        return feat_vals


    def calculate_score(self, responses, trace_check=False):
        """
        Calculate score: spike train distance between given response and target
        spike train.

        Call graph
        ----------
        
        ephys.evaluators.CellEvaluator.evaluate_with_dicts(param_dict)
        |- responses = evaluator.run_protocols
        |- evaluator.objectivescalculator.calculate_scores(responses)
           `- objective.calculate_feature_scores(responses)
           `- feature.calculate_score(responses)
        """

        efel_trace = self._construct_efel_trace(responses)

        if efel_trace is None:
            score = self.max_score
        else:
            feat_name = self.metric_name
            score_func = self.CALC_SCORE_FUNCS[feat_name]
            score = score_func(self, efel_trace, trace_check)

        logger.debug('Calculated score for %s: %f', self.name, score)

        return score


    def __str__(self):
        """
        String representation
        """

        return "%s for %s with stim start %s and end %s, " \
            "and AP threshold override %s" % \
            (self.metric_name,
             self.recording_names,
             self.stim_start,
             self.stim_end,
             self.threshold)

