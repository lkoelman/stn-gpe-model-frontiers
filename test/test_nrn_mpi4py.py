"""
Test running NEURON with MPI.

If NEURON is compiled correctly with MPI support, using the same MPI
libs as mpi4py, both mpi4py and NEURON's ParallelContext should be on the same
rank ("I am X" should be the same on each output line of script).

@see    For more tests of MPI with neuron see NEURON folder /nrn/src/parallel/
"""
from mpi4py import MPI
from neuron import h
pc = h.ParallelContext()
s = "mpi4py thinks I am %d of %d,\
 NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, \
           pc.id(),pc.nhost())
pc.done()