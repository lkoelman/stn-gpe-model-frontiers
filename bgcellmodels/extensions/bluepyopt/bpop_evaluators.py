"""
Extensions to BluePyOpt optimization-related classes

@author Lucas Koelman

@date   13/09/2017


"""


import bluepyopt.ephys as ephys

import math
import logging
logger = logging.getLogger('bpop_ext')




class SumOfSquaresObjective(ephys.objectives.EFeatureObjective):
    """
    Objective that calculates sum of squares of EFeature scores
    """

    def __init__(self, name, features):
        """
        Constructor

        Args:
            name (str): name of this object
            features (list of EFeatures): eFeatures in the objective
        """

        super(SumOfSquaresObjective, self).__init__(name, features)


    def calculate_score(self, responses):
        """
        Objective score
        """

        feature_scores = self.calculate_feature_scores(responses)

        sumsq = sum((score**2 for score in feature_scores))
        return sumsq

class WeightedSumOfSquaresObjective(ephys.objectives.EFeatureObjective):
    """
    Objective that calculates sum of squares of EFeature scores
    """

    def __init__(self, name, features, weights):
        """
        Constructor

        Args:
            name (str): name of this object
            features (list of EFeatures): eFeatures in the objective
            weights (list of float): weights of the eFeatures, in same order
        """

        super(WeightedSumOfSquaresObjective, self).__init__(name, features)
        self.weights = weights


    def calculate_score(self, responses):
        """
        Objective score
        """

        feature_scores = self.calculate_feature_scores(responses)

        score = 0.0
        for feature_score, weight in zip(feature_scores, self.weights):
            score += weight * (feature_score**2)

        return score


class RootMeanSquareObjective(ephys.objectives.EFeatureObjective):
    """
    Objective that calculates sum of squares of EFeature scores
    """

    def __init__(self, name, features):
        """
        Constructor

        Args:
            name (str): name of this object
            features (list of EFeatures): eFeatures in the objective
        """

        super(RootMeanSquareObjective, self).__init__(name, features)


    def calculate_score(self, responses):
        """
        Objective score
        """

        feature_scores = self.calculate_feature_scores(responses)

        rms = math.sqrt(
                    sum((score**2 for score in feature_scores)) / len(feature_scores))
        return rms


class CellEvaluatorCaching(ephys.evaluators.CellEvaluator):
    """
    CellEvaluator extension that can save responses after
    distributed evaluation.
    """

    def set_responses_filename(self, folder, prefix):
        """
        Set directory to save responses to and prefix
        for responses filename.
        """
        self.responses_folder = folder
        self.responses_prefix = prefix


    def evaluate_with_dicts(self, param_dict=None, ind_suffix=None):
        """Run evaluation with dict as input and output"""

        if self.fitness_calculator is None:
            raise Exception(
                'CellEvaluator: need fitness_calculator to evaluate')

        logger.debug('Evaluating %s', self.cell_model.name)

        responses = self.run_protocols(
            self.fitness_protocols.values(),
            param_dict)

        # Save responses if folder was set
        folder = getattr(self, 'responses_folder', None)
        prefix = getattr(self, 'responses_prefix', '')

        if folder is not None and ind_suffix is not None:
            import os.path, pickle
            fname = os.path.join(folder, prefix+str(ind_suffix))
            with open(fname, 'wb') as f:
                pickle.dump(responses, f)

        return self.fitness_calculator.calculate_scores(responses)


    def evaluate_with_lists(self, param_list=None, ind_suffix=None):
        """
        Run evaluation with lists as input and outputs

        @param  param_list  list of parameter values in order corresponding to
                            self.param_names
        """

        param_dict = self.param_dict(param_list)

        obj_dict = self.evaluate_with_dicts(
                            param_dict=param_dict, 
                            ind_suffix=ind_suffix)

        return self.objective_list(obj_dict)