# Module Description

Contains all optimization routines, after morphology has been reduced.

# Module Contents

Modules prefixed by `bpop_` contain functions and classes for doing
optimizations using [BluePyOpt](https://github.com/BlueBrain/BluePyOpt)


Modules prefixed by `neurotune_` are for optimization using 
[Neurotune](https://github.com/NeuralEnsemble/neurotune).


Modules prefixed by `praxis_` use the NEURON built-in 'praxis' method.


--------------------------------------------------------------------------------
# Notes on extending BluePyOpt

Different approaches to creating a cell model

- create own ephys.Model or ephys.models.CellModel

    + see CellModel class
        * https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/ephys/models.py
    
    + see external use cases
        * https://github.com/BlueBrain/BluePyOpt/wiki
        * particularly the example at https://github.com/apdavison/BluePyOpt/blob/pynn-models/bluepyopt/ephys_pyNN/models.py

    + see example of complex model at https://github.com/BlueBrain/BluePyOpt/tree/master/examples/l5pc

    + see dummy cell model at https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/tests/test_ephys/testmodels/dummycells.py