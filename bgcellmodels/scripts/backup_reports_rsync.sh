#!/bin/bash -l


target_dir="~/cloudstore_m/simdata/LuNetSGS"
source_dir="/home/luye/storage"

inc_patterns=("*.html" "*.ipynb" "*.pkl")


rsync_cmd='rsync -vrm --include="*/" --exclude="*"'

for pattern in "${inc_patterns[@]}"; do
    rsync_cmd='${rsync_cmd} --include="${pattern}"'
done

rsync_cmd="${rsync_cmd} ${source_dir} ${target_dir}"
eval $rsync_cmd

# Example command:
# rsync -vrm --include="*/" --include="*.html" --include="*AUTO.pkl" --exclude="*" simdata_newsonic/ ~/cloudstore_m/simdata/LuNetDBS/