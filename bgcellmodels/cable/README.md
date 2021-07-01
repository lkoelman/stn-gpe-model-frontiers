Numerical methods, and mathematical analysis tools relating to the cable equation.


# References

## Online Resources

Example simulators:

- [BRIAN algorithm notes](https://github.com/brian-team/brian/blob/master/brian/experimental/morphology/algorithms.txt)

- [J. Rieke's simple SciPy solver](https://github.com/jrieke/NeuroSim)


Reading out NEURON state in Python:

- https://github.com/WillemWybo/SGF_formalism

Example of converting NMODL to Python code: [BlueBrain NMODL tutorial](https://bluebrain.github.io/nmodl/html/notebooks/nmodl-python-tutorial.html)

## Articles: Matrix Equations Branched Cable


+ Hines (1984) - Efficient Computation of Branched Nerve Equations

+ Hines (1989) - A Program for Simulation of Nerve Equations with Branching Geometries

+ Mascagni (1991) - A Parallelizing Algorithm for Computing Solutions to Arbitrarily Branched Cable Neuron Models

+ Book Koch (1999) - The Biophysics of Computation - Appendix C - Sparse Matrix Methods for Modeling Single Neurons (p. 487)

+ Book Koch, Segev (1998) - Methods in Neural Modeling - Ch. 14 - Numerical Methods for Neuronal Modeling

## Extending NEURON

- see `KSChan.cpp`

- see Neuron's Python wrapper code

- use Python as glue code using `nonvint_block_supervisor.py`
    + explained by hines [here](https://www.neuron.yale.edu/phpBB/viewtopic.php?f=8&t=3392)
    + analogous to how RxD module is implemented (see `nrn/share/lib/python/neuron/rxd`)
    + en example of implementing a channel like this can be found [here](https://github.com/nrnhines/nrntest/blob/master/nrniv/rxd/kchanarray.py)
        * `nonvint_block_supervisor.register(<list of predefined functions>)`
