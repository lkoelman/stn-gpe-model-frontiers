"""
Modification of bluepyopt.deapext.optimisations and bluepyopt.deapext.algorithms
modules to record different statistics.

@author Lucas Koelman

@date   02-10-2017
"""

import numpy as np
import random
import logging
import functools

import deap
import deap.base
import deap.algorithms
import deap.tools


from bluepyopt.deapext import tools
# from bluepyopt.deapext.optimisations import WeightedSumFitness, WSListIndividual

import bluepyopt.optimisations

logger = logging.getLogger('__main__')


import sys
from operator import truediv


class SumOfSquaredDistFitness(object):
    """
    This is a reimplementation of deap.base.Fitness object
    (https://github.com/BlueBrain/deap/blob/bbp-specific/deap/base.py#L121)
    when the values that are assigned represent distances that can be positive
    or negative and should be close to zero for a good fitness

    @note   self.values are assigned from call to toolbox.evaluate()


    @note   The HallOfFame uses fitness comparison operators implemented here

    @note   The genetic algorithm selector implements its own fitness comparisons

                - IBEA selector uses self.wvalues

                - NSGA2 selector uses self.wvalues and self.dominates

    @note   Feature scores are distances that can be positive or negative, 
            of which the absolute values represent a COST rather than
            a FITNESS, hence a fitness value is -1.0*abs(dist)

    @note   self.fitness.values is a python property: when it is assigned a list of values,
            they are multiplied with self.weights and stored in self.wvalues


    TODO: currently weights are squared when evaluating squared sum,
              (weights are incorporated int efeat scores assigned via setValues).
              This can be fixed by dividing by weight before squaring (pass
              obj_weights to __init__)

    TODO: weights are squared in SumOfSquaredDistFitness which is used in HOF, 
          but not in selectors
    """

    weights = None
    wvalues = ()

    def __init__(self, values=(), obj_size=None):
        """
        @param  obj_size     number of objectives
        """

        self.weights = [-1.0] * obj_size if obj_size is not None else [-1]

        if len(values) > 0:
            self.values = values


    def getValues(self):
        """
        Returns absolute values of distances
        """
        return tuple(map(truediv, self.wvalues, self.weights))


    def setValues(self, dist_values):
        """
        Transforms distance values into wvalues that can be interpreted as fitnesses
        (i.e. higher is better).
        """
        try:
            self.wvalues = tuple((w*abs(d) for w,d in zip(self.weights, dist_values)))
        except TypeError:
            _, _, traceback = sys.exc_info()
            raise TypeError, ("Both weights and assigned values must be a "
                              "sequence of numbers when assigning to values of "
                              "%r. Currently assigning value(s) %r of %r to a fitness with "
                              "weights %s."
                              % (self.__class__, dist_values, type(dist_values), self.weights)), traceback

    def delValues(self):
        self.wvalues = ()

    values = property(getValues, setValues, delValues,
                      ("Fitness values. Use directly ``individual.fitness.values = values`` "
                       "in order to set the fitness and ``del individual.fitness.values`` "
                       "in order to clear (invalidate) the fitness. The (unweighted) fitness "
                       "can be directly accessed via ``individual.fitness.values``."))


    @property
    def neg_squared_sum(self):
        """
        Sum of squares of weighted values, negated.
        """
        return -1.0 * sum((v**2 for v in self.wvalues)) # wvalues are negative!


    @property
    def sum(self):
        """
        Return sum of (positive) distances

        @note   usef for keeping statistics / logs
        """
        return -sum(self.wvalues)


    def dominates(self, other, obj=slice(None)):
        """
        Return true if each objective of *self* is not strictly worse than
        the corresponding objective of *other* and at least one objective is
        strictly better.

        :param obj: Slice indicating on which objectives the domination is
                    tested. The default value is `slice(None)`, representing
                    every objectives.
        """
        not_equal = False
        for self_wvalue, other_wvalue in zip(self.wvalues[obj], other.wvalues[obj]):
            if self_wvalue > other_wvalue:
                not_equal = True
            elif self_wvalue < other_wvalue:
                return False
        return not_equal


    @property
    def valid(self):
        """Assess if a fitness is valid or not."""
        return len(self.wvalues) != 0

    def __hash__(self):
        return hash(self.wvalues)

    def __le__(self, other):
        """
        This fitness is _worse_ than other if the negated (squared) sum
        of distances is smaller (more negative) than the other fitness.
        """
        return self.neg_squared_sum <= other.neg_squared_sum

    def __lt__(self, other):
        return self.neg_squared_sum < other.neg_squared_sum

    def __gt__(self, other):
        return not self.__le__(other)

    def __ge__(self, other):
        return not self.__lt__(other)

    def __eq__(self, other):
        return self.wvalues == other.wvalues

    def __ne__(self, other):
        return not self.__eq__(other)


    def __deepcopy__(self, _):
        """
        Override deepcopy
        """

        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __str__(self):
        """Return the values of the Fitness object."""
        return str(self.values if self.valid else tuple())

    def __repr__(self):
        """Return the Python code to build a copy of the object."""
        return "%s.%s(%r)" % (self.__module__, self.__class__.__name__,
                              self.values if self.valid else tuple())


class ScoresDictIndividual(dict):
    """
    Individual consisting of fields

        fitness: SumOfSquaredDistFitness

        genes:  list()
    """

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, *args, **kwargs):
        """
        Constructor

        @param  args             list of parameter values

        @param  'obj_names'      list of objective names

        @note   any object passed in kwargs will be available 
                using dot notation

        REQUIREMENTS

            - must implement Sequence interface for DEAP toolbox
                - must be able to get/set param values (genes) using indexing
                - len() must correspond to number of parameters
            
            - must be compatible with _evaluate_invalid_fitness()

        """

        # This sets it as dictionary item, which works with serialization
        kwargs['fitness'] = SumOfSquaredDistFitness(obj_size=kwargs['obj_size'])
        # del kwargs['obj_size']
        if 'fit_vals' in kwargs:
            kwargs['fitness'].values = kwargs['fit_vals']
            # del kwargs['fit_vals']
        
        param_values = args[0]
        kwargs['param_size'] = len(param_values)
        for i, pval in enumerate(param_values):
            kwargs[i] = pval
        
        super(ScoresDictIndividual, self).__init__(**kwargs)


    def get_param_values(self):
        return [self[i] for i in xrange(self['param_size'])]


    # Implement Sequence interface for DEAP toolbox
    def __len__(self):
        """
        Return number of parameters (length of 'DNA')
        """
        return self['param_size']

    def __iter__(self):
        """
        Iterate over param values (genes)
        """
        return iter(self.get_param_values())

    def __reduce__(self):
        """
        Support serialization using pickle module.
        """
        return (
            deserialize_ind, 
            ([self.get_param_values()], {
                'obj_size': self['obj_size'],
                'fit_vals': self.fitness.values}))

    def __eq__(self, other):
        return np.allclose(
                    self.get_param_values(), 
                    other.get_param_values(),
                    atol=0.0, rtol=1e-6)


def deserialize_ind(args, kwargs):
    return ScoresDictIndividual(*args, **kwargs)


class DEAPOptimisation(bluepyopt.optimisations.Optimisation):

    """DEAP Optimisation class"""

    def __init__(self, evaluator=None,
                 use_scoop=False,
                 seed=1,
                 offspring_size=10,
                 eta=10,
                 selector_name=None,
                 mutpb=1.0,
                 cxpb=1.0,
                 hofs=None,
                 map_function=None,
                 inc_default_params=True
                 ):
        """Constructor

        Args:
            
            evaluator (Evaluator): Evaluator object
            
            seed (float): Random number generator seed
            
            offspring_size (int): Number of offspring individuals in each
                generation
            
            eta (float): Parameter that controls how far the crossover and
                mutation operator disturbe the original individuals
            
            mutpb (float): Mutation probability
            
            cxpb (float): Crossover probability
            
            map_function (function): Function used to map (parallelise) the
                evaluation function calls
            
            hof (hof): Hall of Fame object
            
            selector_name (str): The selector used in the evolutionary
                algorithm, possible values are 'IBEA' or 'NSGA2'

            inc_default_params (bool): whether to include individual with default
                parameter values in population
        """

        super(DEAPOptimisation, self).__init__(evaluator=evaluator)

        self.use_scoop = use_scoop
        self.seed = seed
        self.offspring_size = offspring_size
        self.eta = eta
        self.cxpb = cxpb
        self.mutpb = mutpb
        self.map_function = map_function
        self.include_default_params = inc_default_params

        self.selector_name = selector_name
        if self.selector_name is None:
            self.selector_name = 'IBEA'

        if hofs is None:
            # HallOfFame that uses ind.fitness.__le__ operator for selection
            hof_leastsquares = deap.tools.HallOfFame(50)
            # hallOfFame that uses deap.base.Fitness.dominates operator for selection
            hof_pareto = deap.tools.ParetoFront()
            hofs = [hof_leastsquares, hof_pareto]
        self.hofs = hofs

        # Create a DEAP toolbox
        self.toolbox = deap.base.Toolbox() # see https://github.com/BlueBrain/deap/blob/bbp-specific/deap/base.py

        self.setup_deap()


    def setup_deap(self):
        """Set up optimisation"""

        # Number of objectives
        OBJ_SIZE = len(self.evaluator.objectives)

        # Set random seed
        random.seed(self.seed)

        # Eta parameter of crossover / mutation parameters
        # Basically defines how much they 'spread' solution around
        # The lower this value, the more spread
        ETA = self.eta

        # Number of parameters
        IND_SIZE = len(self.evaluator.params)

        # Bounds for the parameters

        LOWER = []
        UPPER = []

        for parameter in self.evaluator.params:
            LOWER.append(parameter.lower_bound)
            UPPER.append(parameter.upper_bound)

        # ordered_param_names = [p.name for p in self.evaluator.params]
        # ordered_objective_names = [o.name for o in self.evaluator.objectives]

        # Define a function that will uniformly pick an individual
        # def uniform(lower_list, upper_list, dimensions):
            
        #     Sample params to fill initial population list

        #     @param  lower_list  lower bound for each parameter

        #     @param  upper_list  upper bound for each parameter

        #     @return     param values in order returned by CellEvaluator.params,
        #                 which is the same as its param_names constructor argument
            

        #     if hasattr(lower_list, '__iter__'):
        #         return [random.uniform(lower, upper) for lower, upper in
        #                 zip(lower_list, upper_list)]
        #     else:
        #         return [random.uniform(lower_list, upper_list)
        #                 for _ in range(dimensions)]

        # Register the 'uniform' function
        # self.toolbox.register("uniformparams", uniform, LOWER, UPPER, IND_SIZE)

        # Register the individual format
        # An indiviual is create by WSListIndividual and parameters
        # are initially picked by 'uniform'


        # Function that generates initial population
        def gen_initial_pop(
                n=None, 
                init_param_vectors=None, 
                param_sampler=None,
                obj_size=None):
            """
            Generate initial population

            @param  popsize     population size

            @param  
            """
            new_samples = [param_sampler() for _ in range(n-len(init_param_vectors))]
            param_samples = init_param_vectors + new_samples
            pop = [
                ScoresDictIndividual(pvals, obj_size=obj_size)
                    for pvals in param_samples
            ]
            return pop

        # Function that picks parameters for an individual in initial population
        def uniform_sampler():
            """
            Sample parameters uniformly between lower and upper bound.
            """
            return np.random.uniform(LOWER, UPPER)

        # Should default parameters be included as an individual in initial pop?
        if self.include_default_params:
            initial_param_vectors = [
                [p.value for p in self.evaluator.params]
            ]
        else:
            initial_param_vectors = []

        # Register the population generator function
        self.toolbox.register(
            "population",
            gen_initial_pop,
            init_param_vectors=initial_param_vectors,
            param_sampler=uniform_sampler,
            obj_size=OBJ_SIZE)
        

        # Register the evaluation function for the individuals
        # import deap_efel_eval1
        # NOTE: the result of this evaluation is assigned to individual.fitness.values
        self.toolbox.register(
            "evaluate", 
            self.evaluator.evaluate_with_lists)

        # Register the mate operator
        self.toolbox.register(
            "mate",
            deap.tools.cxSimulatedBinaryBounded,
            eta=ETA,
            low=LOWER,
            up=UPPER)

        # Register the mutation operator
        self.toolbox.register(
            "mutate",
            deap.tools.mutPolynomialBounded,
            eta=ETA,
            low=LOWER,
            up=UPPER,
            indpb=0.5)

        # Register the variate operator
        self.toolbox.register("variate", deap.algorithms.varAnd)

        # Register the selector (picks parents from population)
        if self.selector_name == 'IBEA':
            self.toolbox.register("select", tools.selIBEA)
        elif self.selector_name == 'NSGA2':
            self.toolbox.register("select", deap.tools.emo.selNSGA2)
        else:
            raise ValueError('DEAPOptimisation: Constructor selector_name '
                             'argument only accepts "IBEA" or "NSGA2"')

        def _reduce_method(meth):
            """Overwrite reduce"""
            return (getattr, (meth.__self__, meth.__func__.__name__))
        import copyreg
        import types
        copyreg.pickle(types.MethodType, _reduce_method)

        if self.use_scoop:
            if self.map_function:
                raise Exception(
                    'Impossible to use scoop is providing self '
                    'defined map function: %s' %
                    self.map_function)

            from scoop import futures
            self.toolbox.register("map", futures.map)

        elif self.map_function:
            self.toolbox.register("map", self.map_function)


    def run(self,
            max_ngen=10,
            offspring_size=None,
            continue_cp=False,
            cp_filename=None,
            cp_frequency=1):
        """
        Run optimisation
        """
        # Allow run function to override offspring_size
        # TODO probably in the future this should not be an object field anymore
        # keeping for backward compatibility
        if offspring_size is None:
            offspring_size = self.offspring_size

        # Generate the population object
        pop = self.toolbox.population(n=offspring_size)

        stats = deap.tools.Statistics(key=lambda ind: ind.fitness.sum)
        import numpy
        stats.register("avg", numpy.mean)
        stats.register("std", numpy.std)
        stats.register("min", numpy.min)
        stats.register("max", numpy.max)

        pop, log, history = eaAlphaMuPlusLambdaCheckpoint(
            pop,
            self.toolbox,
            offspring_size,
            self.cxpb,
            self.mutpb,
            max_ngen,
            stats=stats,
            halloffames=self.hofs,
            cp_frequency=cp_frequency,
            continue_cp=continue_cp,
            cp_filename=cp_filename)

        return pop, self.hofs, log, history


################################################################################
"""
Algorithm class

Modified from bluepyopt.deapext.algorithms
"""

import pickle


def _evaluate_invalid_fitness(toolbox, population):
    '''
    Evaluate the individuals with an invalid fitness
    Returns the count of individuals with invalid fitness

    @note   fitnesses are invalidated after applying crossover + mutation
            to populations
    '''
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    ind_genes = [ind.get_param_values() for ind in invalid_ind]
    fitnesses = toolbox.map(toolbox.evaluate, ind_genes)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    return len(invalid_ind)


def _update_history_and_hof(halloffames, history, population):
    '''
    Update the hall of fame with the generated individuals

    Note: History and Hall-of-Fame behave like dictionaries

    @param  halloffames     HallOfFame objects with different selection criteria
    '''
    for hof in halloffames:
        hof.update(population)

    history.update(population)


def _record_stats(stats, logbook, gen, population, invalid_count):
    '''
    Update the statistics with the new population
    '''
    record = stats.compile(population) if stats is not None else {}
    logbook.record(gen=gen, nevals=invalid_count, **record)


def _get_offspring(parents, toolbox, cxpb, mutpb):
    '''
    return the offsprint, use toolbox.variate if possible
    '''
    if hasattr(toolbox, 'variate'):
        return toolbox.variate(parents, toolbox, cxpb, mutpb)
    return deap.algorithms.varAnd(parents, toolbox, cxpb, mutpb)


def eaAlphaMuPlusLambdaCheckpoint(
        population,
        toolbox,
        mu,
        cxpb,
        mutpb,
        ngen,
        stats=None,
        halloffames=None,
        cp_frequency=1,
        cp_filename=None,
        continue_cp=False):
    r"""This is the :math:`(~\alpha,\mu~,~\lambda)` evolutionary algorithm

    Args:
        population(list of deap Individuals)
        
        toolbox(deap Toolbox)
        
        mu(int): Total parent population size of EA
        
        cxpb(float): Crossover probability
        
        mutpb(float): Mutation probability
        
        ngen(int): Total number of generation to run
        
        stats(deap.tools.Statistics): generation of statistics
        
        halloffame(deap.tools.HallOfFame): hall of fame
        
        cp_frequency(int): generations between checkpoints
        
        cp_filename(string): path to checkpoint filename
        
        continue_cp(bool): whether to continue
    """

    if continue_cp:
        # A file name has been given, then load the data from the file
        import cPickle
        cp = None
        with open(cp_filename, "rb") as f:
            # Load last saved checkpoint
            while True:
                try:
                    cp = cPickle.load(f)
                except EOFError:
                    break
        population = cp["population"]
        parents = cp["parents"]
        start_gen = cp["generation"]
        halloffame = cp["halloffame"]
        paretofront = cp["paretofront"]
        halloffames = [halloffame, paretofront]
        logbook = cp["logbook"]
        history = cp["history"]
        random.setstate(cp["rndstate"])
    else:
        # Start a new evolution
        start_gen = 1
        parents = population[:]
        logbook = deap.tools.Logbook()
        logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
        history = deap.tools.History()

        # TODO this first loop should be not be repeated !
        invalid_count = _evaluate_invalid_fitness(toolbox, population)
        _update_history_and_hof(halloffames, history, population)
        _record_stats(stats, logbook, start_gen, population, invalid_count)

    # Begin the generational process
    for gen in range(start_gen + 1, ngen + 1):
        offspring = _get_offspring(parents, toolbox, cxpb, mutpb)

        population = parents + offspring

        invalid_count = _evaluate_invalid_fitness(toolbox, offspring)
        _update_history_and_hof(halloffames, history, population)
        _record_stats(stats, logbook, gen, population, invalid_count)

        # Select the next generation parents
        parents = toolbox.select(population, mu)

        logger.info(logbook.stream)

        if(cp_filename and cp_frequency and
           gen % cp_frequency == 0):
            cp = dict(population=population,
                      generation=gen,
                      parents=parents,
                      halloffame=halloffames[0],
                      paretofront=halloffames[1],
                      history=history,
                      logbook=logbook,
                      rndstate=random.getstate())

            # Append checkpoint to end of pickle file (after previous checkpoints)
            with open(cp_filename, "ab") as f:
                pickle.dump(cp, f)
            
            logger.debug('Wrote checkpoint to %s', cp_filename)

    return population, logbook, history
