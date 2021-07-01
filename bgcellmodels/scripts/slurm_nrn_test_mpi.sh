#!/bin/bash -l
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks-per-node=4   # Number of processors
#SBATCH -t 00:01:00           # Time requested (1 min)

module load anaconda

# GCC + OpenMPI toolchain
# module load gcc openmpi/3.1.4

# Intel toolchain
export FI_PROVIDER=verbs
module load intel/intel-cc intel/intel-mkl intel/intel-mpi

# Python environment with NEURON installed
conda activate --stack localpy27

# Execute Python script with MPI
mpirun -n 4 python nrn_test_mpi.py
