#!/bin/bash -l

model_dir="${HOME}/workspace/bgcellmodels/bgcellmodels/models/network/LuNetStnGpe"
job_script="${model_dir}/batchjobs/runjob_bgnetmodel.sh"

# Number of cluster nodes and processes per nodes
num_nodes=1
num_ppn=8
num_proc=$((num_nodes*num_ppn))

# Config files you want to repeat with different seeds
# configs=( \
#     "config.json" \
# )
outputs_clipboard="sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-3.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-6.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-9.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-12.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-15.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-18.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-21.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-24.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-27.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-30.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-33.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-36.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-39.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-42.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-45.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-48.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-51.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-54.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-57.0-Hz.json
sweeps_f-burst-input/syn-V18_severe_gpe2stn_f-burst-60.0-Hz.json"
readarray -t configs <<< "${outputs_clipboard}"

start_seed=888

for conf in "${configs[@]}"; do
    for seed in {0..0}; do
        qsub_command="qsub ${job_script} \
-l walltime=4:00:00,nodes=${num_nodes}:ppn=${num_ppn} \
-v seed=$((start_seed+seed)),config=${conf},numproc=${num_proc},\
dur=7e3,transientperiod=0.0,writeinterval=10e3"

        echo -e "Submitting qsub command:\n> $qsub_command"
        eval $qsub_command
    done
# for fburst in {5,10,20,25,50,100}; do
#         qsub_command="qsub ${job_script} \
# -l walltime=02:30:00 \
# -v ncell=100,dur=26000,burst=${fburst},config=${conf}"

#         echo -e "Submitting qsub command:\n> $qsub_command"
#         eval $qsub_command
#     done
done

echo -n "Waiting for 5 seconds to check job status"
for i in {1..5}; do
    echo -n "."
    sleep 1
done
echo ""

qstat
