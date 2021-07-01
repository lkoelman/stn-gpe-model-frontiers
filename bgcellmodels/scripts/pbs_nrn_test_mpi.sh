#!/bin/bash -l
#PBS -l nodes=1:ppn=16
#PBS -l walltime=00:01:00
#PBS -j oe
#PBS -o ./output_testmpi.o
#PBS -N myMPItest
#PBS -V

mpirun -n 4 python nrn_test_mpi.py

