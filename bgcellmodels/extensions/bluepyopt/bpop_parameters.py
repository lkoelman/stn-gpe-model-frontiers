"""
Extensions to Ephys.parameters module

@author Lucas Koelman

@date   8/05/2018
"""

import bluepyopt.ephys as ephys

import logging
logger = logging.getLogger('bpop_ext')


class NrnSegmentParameter(ephys.parameters.NrnParameter, ephys.serializer.DictMixin):
    """
    Parameter of a section
    """
    SERIALIZED_FIELDS = ('name', 'value', 'frozen', 'locations', )

    def __init__(
            self,
            name,
            value=None,
            frozen=False,
            param_name=None,
            locations=None):
        """
        Contructor
        Args:
            name (str): name of the Parameter
            value (float): Value for the parameter, required if Frozen=True
            frozen (bool): Whether the parameter can be varied, or its values
            is permently set
            param_name (str): name used within NEURON
            locations (list of ephys.locations.Location): locations on which
                to instantiate the parameter
        """

        super(NrnSegmentParameter, self).__init__(
            name,
            value=value,
            frozen=frozen)

        self.locations = locations
        self.param_name = param_name


    def instantiate(self, sim=None, icell=None):
        """
        Instantiate
        """
        if self.value is None:
            raise Exception(
                'NrnSegmentParameter: impossible to instantiate parameter "%s" '
                'without value' % self.name)

        for location in self.locations:
            isegments = location.instantiate(sim=sim, icell=icell)
            for segment in isegments:
                setattr(segment, self.param_name, self.value)
            logger.debug(
                'Set %s in %s to %s',
                self.param_name,
                location,
                self.value)

    def __str__(self):
        """String representation"""
        return '%s: %s %s = %s' % (self.name,
                                   [str(location)
                                    for location in self.locations],
                                   self.param_name,
                                   self.value if self.frozen else self.bounds)

class NrnScaleRangeParameter(ephys.parameters.NrnParameter, ephys.serializer.DictMixin):
    """
    Parameter that scales a NEURON RANGE parameter 
    (pre-existing spatial distribution) in target region.
    """

    SERIALIZED_FIELDS = ('name', 'value', 'frozen', 'bounds', 'param_name',
                         'value_scaler', 'locations', )

    def __init__(
            self,
            name,
            value=None,
            frozen=False,
            bounds=None,
            param_name=None,
            locations=None,
            segment_filter=None):
        """
        Contructor

        Args:
            name (str): name of the Parameter
            
            value (float): Value for the parameter, required if Frozen=True
            
            frozen (bool): Whether the parameter can be varied, or its values
            is permently set
            
            bounds (indexable): two elements; the lower and upper bounds
                (Optional)
            
            param_name (str): name used within NEURON
            
            locations (list of ephys.locations.Location): locations on which
                to instantiate the parameter
        """

        super(NrnScaleRangeParameter, self).__init__(
            name,
            value=value,
            frozen=frozen,
            bounds=bounds)

        self.locations = locations
        self.param_name = param_name
        self.segment_filter = segment_filter


    def instantiate(self, sim=None, icell=None):
        """
        Instantiate (i.e. apply the parameter)
        """

        if self.value is None:
            raise Exception(
                'NrnRangeParameter: impossible to instantiate parameter "%s" '
                'without value' % self.name)

        for location in self.locations:
            for isection in location.instantiate(sim=sim, icell=icell):
                for seg in isection:
                    # Skip segment if doesn't match filter
                    if (self.segment_filter is not None) and (not self.segment_filter(seg)):
                        continue
                    # Scale the parameter
                    old_val = getattr(seg, '%s' % self.param_name)
                    new_val = old_val * self.value
                    setattr(seg, '%s' % self.param_name, new_val)
        
        logger.debug(
                'Scaled %s in %s by factor %s', self.param_name,
                [str(location) for location in self.locations],
                self.value)


    def __str__(self):
        """String representation"""
        return '%s: %s %s = %s' % (self.name,
                                   [str(location)
                                    for location in self.locations],
                                   self.param_name,
                                   self.value if self.frozen else self.bounds)


class NrnOffsetRangeParameter(ephys.parameters.NrnParameter, ephys.serializer.DictMixin):
    """
    Parameter that offsets a NEURON RANGE parameter
    (pre-existing spatial distribution) in target region.
    """

    SERIALIZED_FIELDS = ('name', 'value', 'frozen', 'bounds', 'param_name',
                         'value_scaler', 'locations', )

    def __init__(
            self,
            name,
            value=None,
            frozen=False,
            bounds=None,
            param_name=None,
            locations=None,
            threshold=None):
        """
        Contructor

        Args:
            name (str):         name of the Parameter
            
            value (float):      Value for the parameter, required if Frozen=True
            
            frozen (bool):      Whether the parameter can be varied, or its values
            is permently set
            
            bounds (indexable): two elements; the lower and upper bounds
                                (Optional)
            
            param_name (str):   name used within NEURON
            
            locations (list of ephys.locations.Location):
                                locations on which to instantiate the parameter

            threshold:          threshold on parameter value: only offset if parameter
                                is larger than this value.
        """

        super(NrnOffsetRangeParameter, self).__init__(
            name,
            value=value,
            frozen=frozen,
            bounds=bounds)

        self.locations = locations
        self.param_name = param_name
        self.threshold = threshold


    def instantiate(self, sim=None, icell=None):
        """
        Instantiate (i.e. apply the parameter)
        """

        if self.value is None:
            raise Exception(
                'NrnRangeParameter: impossible to instantiate parameter "%s" '
                'without value' % self.name)

        for location in self.locations:
            for isection in location.instantiate(sim=sim, icell=icell):
                for seg in isection:
                    # Scale the parameter
                    old_val = getattr(seg, '%s' % self.param_name)
                    if old_val > self.threshold:
                        new_val = old_val + self.value
                        setattr(seg, '%s' % self.param_name, new_val)
        
        logger.debug(
            'Offset %s in %s by value %s', self.param_name,
            [str(location)
             for location in self.locations],
            self.value)


    def __str__(self):
        """String representation"""
        return '%s: %s %s = %s' % (self.name,
                                   [str(location)
                                    for location in self.locations],
                                   self.param_name,
                                   self.value if self.frozen else self.bounds)


class DistanceParamaterizedRangeVar(ephys.parameters.NrnParameter, ephys.serializer.DictMixin):
    """
    Parameter that sets a NEURON RANGE parameter according to a parameterized
    function/distribution.

    @see    Similar functionality can be obtained using:
            - classes in bluepyopt/ephys/parameterscalers.py
            - example in bluepyopt/tests/test_ephys/test_parameters.py#L60
    """

    SERIALIZED_FIELDS = ('name', 'value', 'frozen', 'bounds',
                        'param_name', 'dist_func', 'dist_func_params',
                        'locations')

    def __init__(
            self,
            name,
            param_name=None,
            locations=None,
            dist_func=None,
            **dist_func_params):
        """
        Contructor

        @param  param_name : str
                Name of NEURON RANGE parameter, e.g. 'gbar_NaF'

        @param  dist_func : function
                Function defined at top-level of module. Its signature should
                be the distance from soma as first, ordered parameter, and
                keyword arguments as subsequent parameters. I.e.:
                `dist_func(x, param1=None, param2=None, ...)`

        @param  dist_func_params : **kwargs
                Keyword arguments for the distribution function

        @pre    parameter names in 'dist_func_params' must not clash with
                attribute names of this class

        EXAMPLE
        -------

        def gaussian(x, mu=0, sigma=2):
            return np.exp(-(x-mu)**2/(2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)

        gkbar_distr = ParameterizedRangeVar(
                        name        = 'gmax_sKCa_distribution',
                        param_name  = 'gk_sKCa'
                        dist_func   = gaussian,
                        mu          = 0.0,
                        sigma       = 1.0)

        mu_param = ephys.parameters.MetaParameter(
                        name        = 'gkbar_mu',
                        obj         = gkbar_distr,
                        attr_name   = 'mu',
                        value       = 0.0,
                        bounds      = [-0.5, 0.5],
                        frozen      = False)

        """

        super(DistanceParamaterizedRangeVar, self).__init__(
            name,
            value=1,        # unused
            frozen=True,    # unused, freeze metaparameters
            bounds=[0, 2])  # unused

        self.param_name = param_name
        self.locations = locations
        self.dist_func = dist_func
        self.dist_func_params = dist_func_params

        # Attributes set by the metaparameters
        for k,v in self.dist_func_params.items():
            setattr(self, k, v)


    def instantiate(self, sim=None, icell=None):
        """
        Instantiate (i.e. apply the parameter)
        """

        soma = icell.soma[0]
        sim.neuron.h.distance(0, 0.5, sec=soma)

        for k in self.dist_func_params.keys():
            self.dist_func_params[k] = getattr(self, k)

        for location in self.locations:
            for sec in location.instantiate(sim=sim, icell=icell):
                for segment in sec:
                    distance = sim.neuron.h.distance(1, segment.x, sec=segment.sec)
                    value = self.dist_func(distance, **self.dist_func_params)
                    setattr(segment, self.param_name, value)


    def __str__(self):
        """String representation"""
        return "{} with distribution function {} and parameters: {}".format(
                self.name, self.dist_func, self.dist_func_params)