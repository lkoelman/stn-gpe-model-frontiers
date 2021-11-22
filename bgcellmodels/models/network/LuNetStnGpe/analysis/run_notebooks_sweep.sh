#!/bin/bash -l

model_dir="${HOME}/workspace/bgcellmodels/bgcellmodels/models/network/KlmnNetMorpho"
nb_filename="synchrony_analysis_auto.ipynb"
nb_dir="${model_dir}/analysis"
notebook="${nb_dir}/${nb_filename}"
logfile="${nb_dir}/nb_exec_list.log" # change for copies of this script
conffile="${nb_dir}/nb_exec_conf.py" # change for copies of this script
export NB_CONF_FILE=${conffile} # for subprocesses

# List simulation output directories to analyze
outputs_clipboard="/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/q14_sweep-gpe-surrogate-frac/2018.08.22_job-783224.sonic-head_DA-depleted-v3_CTX-f0_GPE-surrogates-frac0.0
/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/q14_sweep-gpe-surrogate-frac/2018.08.22_job-783359.sonic-head_DA-depleted-v3_CTX-f2.5-GPE-surrogates-frac0.0
/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/q14_sweep-gpe-surrogate-frac/2018.08.22_job-783360.sonic-head_DA-depleted-v3_CTX-f2.5-GPE-surrogates-frac0.1
/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/q14_sweep-gpe-surrogate-frac/2018.08.22_job-783361.sonic-head_DA-depleted-v3_CTX-f2.5-GPE-surrogates-frac0.2
/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/q14_sweep-gpe-surrogate-frac/2018.08.22_job-783362.sonic-head_DA-depleted-v3_CTX-f2.5-GPE-surrogates-frac0.5
/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/q14_sweep-gpe-surrogate-frac/2018.08.22_job-783363.sonic-head_DA-depleted-v3_CTX-f2.5-GPE-surrogates-frac0.75
/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/q14_sweep-gpe-surrogate-frac/2018.08.22_job-783364.sonic-head_DA-depleted-v3_CTX-f2.5-GPE-surrogates-frac0.9"
readarray -t output_dirs <<< "${outputs_clipboard}"
# output_dirs=(${outputs_clipboard//$'n'// }) # substitute newlines and make array
# output_dirs=( \
#     "/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/syn-v2_q7_sweep-transmission-delay/2018.08.03_job-780788.sonic-head_DA-depleted-v3_CTX-f0_StnGpe-d4.0_GpeStn-d6.0" \
# )

# Generate output dirs from parent directory
# parent_dir="/run/media/luye/Windows7_OS/Users/lkoelman/simdata-win/syn-v2_q3_const-freq_vary-burst-frac"
# dir_contents=$(ls -d -1 ${parent_dir}/**)
# output_dirs=($(echo ${dir_contents})) # echo without quote turns newlines into spaces

num_outputs=${#output_dirs[@]}
# hpfreqs=(3 5 5 5 5 5 5)
hpfreqs=($(for ((i=0;i<num_outputs;i++)); do echo 5; done))
# lpfreqs=(20 20 20 20 20 20 20 20 20)
lpfreqs=($(for ((i=0;i<num_outputs;i++)); do echo 20; done))

# Parameter sweep
sweep_name="gmax_gpe_gpe"
sweep=(0.0 0.0 0.1 0.2 0.5 0.75 0.9)

cd $nb_dir

# for outdir in "${output_dirs[@]}"; do
for ((i=0;i<num_outputs;i++)); do
    outdir=${output_dirs[i]}
    hpfreq=${hpfreqs[i]}
    lpfreq=${lpfreqs[i]}
    outfile="${outdir}/${nb_filename}"

    # Don't run notebook if exists
    # if [ -e ${outfile} ]; then
    #     echo -e "Skipping ${outdir}"
    #     continue
    # fi

    # notebook will read variables from Python configuration script
    echo $outdir >> $logfile
    pyscript="outputs = \"${outdir}\"
hpfreq = ${hpfreq}
lpfreq = ${lpfreq}
sweep_var_name = \"${sweep_name}\"
sweep_var_value = ${sweep[i]}"
    echo -e "${pyscript}" > ${conffile} # (quote preserves newlines)

    # Execute notebook
    # (timeout option is for notebook cell execution, default is 30 seconds)
    nb_exec="jupyter nbconvert --ExecutePreprocessor.timeout=None --execute --to notebook ${notebook} --output-dir=${outdir}"
    echo -e "Executing notebook for output: ${outdir}"
    eval $nb_exec

    echo -e "Finished! Saving notebook to output directory ..."
    nb_save="jupyter nbconvert ${outfile} --template=toc2 --output-dir=${outdir}"
    eval $nb_save
done

echo "Finished executing notebooks!"