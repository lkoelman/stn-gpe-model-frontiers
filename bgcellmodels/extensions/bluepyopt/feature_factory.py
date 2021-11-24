"""
Creation of EFEL features for STN model optimization.

@author Lucas Koelman

@date   3/10/2017

"""
import math

import bluepyopt.ephys as ephys
import bgcellmodels.extensions.bluepyopt.efeatures_spiketrain as espk
import efel

import logging
logger = logging.getLogger('bluepyopt.ephys.efeatures')

def make_features(proto_wrapper):
    """
    Make eFEL features based on a protocol definition and its characterizing 
    features.

    @return     a dictionary {feature_names: feature_objects}
    """
    protocol = proto_wrapper.ephys_protocol

    candidate_feats = {}
    default_trace = {'': protocol.recordings[0].name} # location : str -> recording : str

    eFEL_available_features = efel.getFeatureNames()
    custom_spiketrain_features = espk.getFeatureNames()

    for feat_name, feat_params in proto_wrapper.characterizing_feats.items():

        # Get interval of response trace from which feature is calculated
        response_interval   = feat_params.get('response_interval',
                                             proto_wrapper.response_interval)
        max_score           = feat_params.get('max_score', None)
        force_max_score     = max_score is not None

        # Identify library providing feature and make it
        if feat_name in eFEL_available_features:

            feature = ephys.efeatures.eFELFeature(
                        name                ='{}.{}'.format(protocol.name, feat_name),
                        efel_feature_name   = feat_name,
                        recording_names     = feat_params.get('traces', default_trace),
                        stim_start          = response_interval[0],
                        stim_end            = response_interval[1],
                        double_settings     = feat_params.get('double', None),
                        int_settings        = feat_params.get('int', None),
                        threshold           = proto_wrapper.spike_threshold,
                        max_score           = max_score,
                        force_max_score     = force_max_score,
                        )

        elif feat_name in custom_spiketrain_features:

            feature = espk.SpikeTrainFeature(
                        name                ='{}.{}'.format(protocol.name, feat_name),
                        metric_name         = feat_name,
                        recording_names     = feat_params.get('traces', default_trace),
                        stim_start          = response_interval[0],
                        stim_end            = response_interval[1],
                        double_settings     = feat_params.get('double', None),
                        int_settings        = feat_params.get('int', None),
                        threshold           = proto_wrapper.spike_threshold,
                        max_score           = max_score,
                        force_max_score     = force_max_score,
                        )

        else:
            raise ValueError('Unknown feature: <{}>'.format(feat_name))

        candidate_feats[feat_name] = feature

    return candidate_feats


def make_opt_features(proto_wrappers):
    """
    Make features that associate each protocol with some relevant metrics.

    @return     dict(StimProtocol : dict(feature_name : tuple(feature, weight)))

                    I.e. a dictionary that maps ephys.protocol objects to another
                    dictionary, that maps feature names to a feature object and
                    its weight for the optimization.

    @note   available features:
                - import efel; efel.getFeatureNames()
                - see http://efel.readthedocs.io/en/latest/eFeatures.html
                - see pdf linked there (ctrl+f: feature name with underscores as spaces)
                - each feature has specific (required_features, required_trace_data, required_parameters)

    """

    proto_feat_dict = {} # StimProtocol -> dict[feat_name : str, tuple[feature, weight]]

    # For each protocol used in optimization: make the Feature objects
    for proto in proto_wrappers:

        proto_feat_dict[proto.IMPL_PROTO] = {} # feature_name -> (feature, weight)

        # Make eFEL features
        candidate_feats = make_features(proto)

        # Save each feature and its weight
        for feat_name in candidate_feats.keys():
            # Weight and normalization factor for feature score (distance)
            feat_weight = proto.characterizing_feats[feat_name]['weight']
            feat_std = proto.characterizing_feats[feat_name].get('exp_std', float('nan'))
            
            # Add to feature dict if nonzero weight
            # NOTE: feature score = sum(feat[i] - exp_mean) / N / exp_std  => so exp_std determines weight (in case of SingletonObjective)
            if feat_weight > 0.0:
                feat_obj = candidate_feats[feat_name]
                feat_obj.exp_std = feat_std
                proto_feat_dict[proto.IMPL_PROTO][feat_name] = feat_obj, feat_weight

    return proto_feat_dict


def calc_feature_targets(protos_feats, protos_responses, remove_problematic=True,
                         raise_check=True):
    """
    Calculate target values for features used in optimization (using full model).


    @param      protocols           BpopProtocolFactory objects


    @param      saved_responses     file path to pickled responses dictionary

    @post       for each EFelFeature in given dictionary, the feature.exp_mean
                and feature.exp_std will be set
    """

    eFEL_available_features = efel.getFeatureNames()
    custom_spiketrain_features = espk.getFeatureNames()

    if raise_check:
        def check_warn(msg): raise ValueError(msg)
    else:
        check_warn = logger.warning

    # Run each protocol and get its responses
    for stim_proto, feat_dict in protos_feats.items():

        # Get response traces
        responses = protos_responses[stim_proto.name]

        # Mark features that were not calculated correctly
        problem_feat_names = []

        # Use response to calculate target value for each features
        for feat_name, feat_data in feat_dict.items():

            # Calculate feature value from full model response
            e_feature, weight = feat_data
            # weight = math.sqrt(weight) # because feature scores will be squared

            target_value = e_feature.calculate_feature(responses, raise_warnings=True)

            if feat_name in eFEL_available_features:

                # Check if value is sane
                if (target_value is None) or math.isinf(target_value) or math.isnan(target_value):
                    logger.warning('Feature {} value {} not sane for protocol {}'.format(
                                feat_name, target_value, stim_proto.name))
                    problem_feat_names.append(feat_name)

                # Now we can set the target value
                # NOTE: distance is sum_i^N(feat[i] - exp_mean) / N / exp_std
                # NOTE: exp_std will have as much influence as weight! This needs to be set in protocols file
                e_feature.exp_mean = target_value

                # if a normalization factor was set, use it, otherwise distance will be relative to initial target
                if math.isnan(e_feature.exp_std) or e_feature.exp_std is None:
                    check_warn(
                        "No standard deviation 'exp_std' specified for efeature {}.".format(e_feature.name))

                # OLD SOLUTION (no weighted objective -> exp_std = 1 / weight)
                #     e_feature.exp_std = abs(target_value / weight)
                # else:
                #     e_feature.exp_std /= abs(weight)

            elif feat_name in custom_spiketrain_features:

                # Set target spike train
                e_feature.set_target_values(target_value)

                # target is not a number so need to decide normalization factor in other way
                # exp_std was set to normalization factor or NaN
                if math.isnan(e_feature.exp_std) or e_feature.exp_std is None:
                    check_warn(
                        'No standard deviation "exp_std" given for eFeature {}. '
                        'Since custom eFeatures may not have a value for a single response '
                        'it is highly recommended to provide a normalization factor.'.format(
                            e_feature.name))

                # OLD SOLUTION (no weighted objective -> exp_std = 1 / weight)
                #     e_feature.exp_std = abs(1.0 / weight)
                # else:
                #     e_feature.exp_std /= abs(weight) # exp_std already set to norm factor

            else:
                raise ValueError('Unknown feature: <{}>'.format(feat_name))

        # Remove problematic features
        if remove_problematic:
            for feat_name in problem_feat_names:
                feat_dict.pop(feat_name)
                logger.debug('Removed feature {}:{}'.format(stim_proto.name, feat_name))


def all_features_weights(feature_dicts):
    """
    Concatenate all EFelFeature objects and all their corresponding weights
    for the objective calculation into two lists.

    @param  feature_dicts       an iterable of dictionaries with feature names as
                                keys and tuples (EFelFeature, weight) as values.
    """
    all_features = []
    all_weights = []
    for featdict in feature_dicts:
        # Add features and weights for this protocol
        feats, weights = zip(*featdict.values()) # values ist list of (feature, weight)
        all_features.extend(feats)
        all_weights.extend(weights)

    return all_features, all_weights