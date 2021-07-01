# Plan

## TODO: write code that assembles matrix & does integration step

- see xmind notes Neuron + Arbor -> refer to source code

- first write simple implementation in Python

- then potentially write one in Rust/C/C++ that gets passed cfunc/ufunc/gufunc

## TODO: write Python mechanism class

Features:

- parse .mod files and convert into Python representation
    + see [Arbor: modcc](https://github.com/eth-cscs/arbor/blob/master/modcc/modcc.cpp)
        * Possibly reuse `Module` and `Parser` object, then use it's state to generate own object instead of generating Cpp code
    + see [CoreNEURON: mod2c](https://github.com/BlueBrain/mod2c)
    + google terms "python domain specific language" "python parser"

- get steady-state values as tables

- get equations as Numba functions (ufunc/lambdify/LLVMjitprinter)

## TODO: write Numba cfunc + integrator backend

See Numba documentation: [Creating C callbacks with @cfunc](http://numba.pydata.org/numba-doc/dev/user/cfunc.html).

Numba provides functions for jit-compiling python methods and turning them into C-callable functions. These can be passed to the simulator object implemented using C ABI interface (exported functions).

See [Numba developer documentation](http://numba.pydata.org/numba-doc/dev/index.html), 
particularly sections 'cfuncs', compiled classes, compiling ahead of time

Alternatively: Sympy expressions can be compiled using Numba backend. See following links:

- [Compiling Sympy expressions](http://docs.sympy.org/latest/modules/numeric-computation.html)
    + different backends for printing/compiling can be found [on GitHub](https://github.com/sympy/sympy/tree/master/sympy/printing)
    + [Blogpost about compiling Sympy expressions](http://matthewrocklin.com/blog/work/2013/03/19/SymPy-Theano-part-1)
    + [LLVM JIT compilation with Numba backend](see https://github.com/sympy/sympy/blob/master/sympy/printing/llvmjitcode.py)
    + [Disscussion of Numba compilation of Sympy expressions](https://github.com/sympy/sympy/issues/8406)


## TODO: write Numba cuda backend

See Numba documentation: [writing CUDA kernels](http://numba.pydata.org/numba-doc/dev/cuda/kernels.html).
