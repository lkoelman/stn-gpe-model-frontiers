#!/bin/bash -l

################################################################################
# QSUB CONFIGURATION
################################################################################

# Set the number of nodes & processors per node
# 24 cores (threads) available per node, so max ppn=24
## PBS -l nodes=1:ppn=8


# Set the walltime of the job to 1 hour (format is hh:mm:ss)
# - 1 jobs/node: 11 [s walltime] / 50 [ms simtime] for LuNetStnGpe with 50 STN + 100 GPe cells
# - 3 jobs/node is (19-23)s / 50 ms for 50 STN + 100 GPe cells
# - to get hh and mm (for hh:mm:ss) do hh = runtime // 3600; mm = runtime % 3600
## PBS -l walltime=00:45:00

# Declares that all environment variables in the qsub command's
# environment are to be exported to the batch job.
#PBS -V

# Working directory where submitted script will be executed ($PBS_O_INITDIR)
# This must be a full path.
#PBS -d /home/people/15203008/workspace/bgcellmodels/bgcellmodels/models/network/LuNetStnGpe

# Specifies the jobname. The default name is the script name (basename)
## PBS -N BG_network_model

# Defines the path to be used for the standard output stream of the batch job.
# The default output is <job_name.job_number> which is pretty clear
## PBS -o ./output_bgnetmodel.o.log
# Same for standard error stream:
## PBS -e ./output_bgnetmodel.e.log

# E-mail on begin (b), abort (a) and end (e) of job
#PBS -m bae

# E-mail address of recipient
#PBS -M lucas.koelman@ucdconnect.ie


################################################################################
# JOB SCRIPT
################################################################################

# Example shell command:
# >>> qsub runjob_lucastest.sh -l walltime=00:45:00 \
# >>> -v ncell=100,dur=10000,seed=15203008,\
# >>> config=~/workspace/bgcellmodels/models/KlmnNetMorpho/configs/simple_config.json,
# >>> outdir=~/storage,id=with_DA_depleted_1


echo -e "
Directory where qsub was called:        $PBS_O_WORKDIR
Working directory for submitted script: $PBS_O_INITDIR
The job id is                           $PBS_JOBID
The job name is                         $PBS_JOBNAME
"

# Setup environment
module load intel-mpi gcc anaconda
source activate localpy27

# Get all the paths
if [ -z "$outdir" ]; then
    outdir="~/storage"
fi

model_dir="${HOME}/workspace/bgcellmodels/bgcellmodels/models/network/LuNetStnGpe"
model_script=model_parameterized.py
model="${model_dir}/${model_script}"
config_file="configs/${config}"
model_config="${model_dir}/${config_file}"

# Command with minimum required arguments
mpi_command="mpirun -n ${numproc} python ${model} \
--scale 1.0 --dur ${dur} \
--nogui --progress --config ${model_config} \
--outdir ${outdir} -id ${PBS_JOBID}"

# Optional arguments passed to python script
opt_names_long=("seed" "writeinterval" "transientperiod")
for optname in "${opt_names_long[@]}"; do
    if [ -n "${!optname}" ]; then
        mpi_command="${mpi_command} --${optname} ${!optname}"
    fi
done

opt_names_short=("wi" "tp")
for optname in "${opt_names_short[@]}"; do
    if [ -n "${!optname}" ]; then
        mpi_command="${mpi_command} -${optname} ${!optname}"
    fi
done


# Optional flags passed to python script
flag_names=("lfp" "no-lfp" "dd" "dnorm")
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
cat "${model_dir}/batchjobs/generate_jobs.sh"

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
