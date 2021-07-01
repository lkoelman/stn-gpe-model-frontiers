#!/bin/bash -l

target_dir=/home/luye/storage

# Tip: copy source file names using ctrl+i, ctrl+c in Nautilus
source_files=( \
    "2018.06.24_job-777379.sonic-head_DA-control-v2_CTX-beta" \
    "2018.06.24_job-777375.sonic-head_DA-control-v2_CTX-beta" \
    "2018.06.24_job-777378.sonic-head_DA-control-v2_CTX-beta" \
    "2018.06.24_job-777377.sonic-head_DA-control-v2_CTX-beta" \
    "2018.06.24_job-777376.sonic-head_DA-control-v2_CTX-beta" \
)

# Script will be interactive and request for password
rsync_cmd="rsync -avz -e ssh 15203008@sonic.ucd.ie:/home/people/15203008/storage/{"


for simfile in "${source_files[@]}"; do
    rsync_cmd="${rsync_cmd}${simfile},"
done

# We don't need the last comma
rsync_cmd="${rsync_cmd::-1}} ${target_dir}"
eval $rsync_cmd
