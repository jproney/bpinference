#!/bin/bash

# location for home directory to work with
folder_dir="/michorlab/jroney/bpinference"

cd $folder_dir

while read -r args
do
 # Command to send to nodes - will use sbatch here
 noise =$(echo ${args}) 
 export noise
 for ((c=1; c<=$nreps; c++))
 do
   sbatch --job-name=one_type_noise \
   $folder_dir/noise-experiments/one_type_noise_benchmark.sbatch 
   echo "${noise}"
 done
done < $folder_dir/noise-experiments/noise_params.txt