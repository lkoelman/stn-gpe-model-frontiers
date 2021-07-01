#!/bin/bash -l

################################################################################
# QSUB CONFIGURATION
################################################################################

# All command line options to `sbatchz and their equivalent environment
# variables can be found at https://slurm.schedmd.com/sbatch.html

# Working directory where submitted script will be executed ($SLURM_SUBMIT_DIR)
# This must be a full path.
#SBATCH -D /home/people/15203008/workspace/bgcellmodels/bgcellmodels/models/network/LuNetStnGpe

# Log files for standard output and standard error
#SBATCH -o LuNetStnGpe-%j-stdout.log.o
#SBATCH -e LuNetStnGpe-%j-stderr.log.e

# E-mail on begin (b), abort (a) and end (e) of job
#SBATCH --mail-type=BEGIN,END,FAIL

# E-mail address of recipient
#SBATCH --mail-user lucas.koelman@ucdconnect.ie


################################################################################
# JOB SCRIPT
################################################################################

echo -e "
Job id is                               $SLURM_JOB_ID
Job name is                             $SLURM_JOB_NAME
The Sonic node is                       $SLURM_NODEID
Number of tasks is                      $SLURM_NTASKS
Working directory for submitted script: $SLURM_SUBMIT_DIR
"

# Setup environment
module purge
module load anaconda

## GCC toolchain:
module load gcc openmpi/3.1.4

## Python environment
conda activate --stack localpy27

# Get all the paths
if [ -z "$outdir" ]; then
    outdir="~/scratch"
fi

# Simulation script
model_dir="${HOME}/workspace/bgcellmodels/bgcellmodels/models/network/LuNetStnGpe"
model_filename=model_parameterized.py
model_filepath="${model_dir}/${model_filename}"

# Command with minimum required arguments
# numproc = $((SLURM_JOB_NUM_NODES*SLURM_NTASKS_PER_NODE))
mpi_command="mpirun -n ${numproc} python ${model_filepath} -id ${SLURM_JOB_ID}"

################################################################################
# Arguments for python script

# ARGUMENTS: short arguments
opt_names_short=("d" "o" "wi" "tp" "p" "dc" "cs" "cc" "ca" "dm")
for optname in "${opt_names_short[@]}"; do
    if [ -n "${!optname}" ]; then
        mpi_command="${mpi_command} -${optname} ${!optname}"
    fi
done


# ARGUMENTS: long arguments
opt_names_long=("dur" "simdt" "scale" "seed" \
    "writeinterval" "transientperiod" "reportinterval" \
    "outdir" "config")
for optname in "${opt_names_long[@]}"; do
    if [ -n "${!optname}" ]; then
        mpi_command="${mpi_command} --${optname} ${!optname}"
    fi
done


# ARGUMENTS: flags
flag_names=("lfp" "nolfp" "dbs" "nodbs" "dd" "dnorm" "progress" "nogui" "gui")
for flagname in "${flag_names[@]}"; do
    if [ -n "${!flagname}" ]; then
        mpi_command="${mpi_command} --${flagname}"
    fi
done

cd $model_dir

# Sanity check
echo -e "
--------------------------------------------------------------------------------
Executing script with following inputs:

- ncell = ${ncell}
- dur = ${dur}
- outdir = ${outdir}

--------------------------------------------------------------------------------
The final command (following '> ') is:

> $mpi_command
"

echo -e "
################################################################################
REPRODUCIBILITY INFO

--------------------------------------------------------------------------------
Model repository version:
"
git log -1

echo -e "
--------------------------------------------------------------------------------
Python package versions:
"
pip freeze

echo -e "
--------------------------------------------------------------------------------
The contents of job generation file is:
"
cat "${model_dir}/batchjobs/slurm_generate_jobs.sh"

echo -e "
--------------------------------------------------------------------------------
The contents of the model configuration file is:
"
cat $model_config

echo -e "
Model output follows below line:
================================================================================
"
# Evaluate our command
eval $mpi_command
