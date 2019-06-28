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
   sbatch -o /dev/null --job-name=logistic_benchmark \
   $folder_dir/cluster_scripts/logistic_benchmark.sbatch 
   echo "${nsims}"
 done
done < $folder_dir/benchmarks/logistic_benchmark_params.txt

