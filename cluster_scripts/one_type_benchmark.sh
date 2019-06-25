#!/bin/bash

# location for home directory to work with
folder_dir="/michorlab/jroney/bpinference"

cd $folder_dir

while read -r args
do
 export args
 # Command to send to nodes - will use sbatch here
 sbatch --job-name=one_type_benchmark \
    $folder_dir/cluster_scripts/one_type_benchmark.sbatch
echo "${args}"
sleep 0.3
done < $folder_dir/benchmarks/one_type_benchmark_params.txt

