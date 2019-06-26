#!/bin/bash

# location for home directory to work with
folder_dir="/michorlab/jroney/bpinference"

cd $folder_dir

while read -r args
do
 # Command to send to nodes - will use sbatch here
 nsims=$(echo ${args} | grep -o '^[0-9]*') 
 nreps=$(echo ${args} | grep -o '[0-9]*$')
 export nsims
 for ((c=1; c<=$nreps; c++))
 do
   sbatch -o /dev/null --job-name=one_type_benchmark_unif \
   $folder_dir/cluster_scripts/one_type_benchmark_unif.sbatch 
   echo "${nsims}"
 done
done < $folder_dir/benchmarks/one_type_benchmark_params.txt

