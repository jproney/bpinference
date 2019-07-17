#!/bin/bash

folder_dir="/michorlab/jroney/bpinference"

cd $folder_dir

while read -r args
do
 # Command to send to nodes - will use sbatch here
 cell_name=$(echo ${args} | grep -o '^[a-z,0-9,A-Z,_,-]*') 
 drug=$(echo ${args} | grep -o '[a-z,0-9,A-Z,_,-]*$')
 export cell_name
 export drug
 sbatch --job-name=lincs \
 $folder_dir/lincs-code/lincs.sbatch 
 echo "${cell_name}"
 echo "${drug}"
done < $folder_dir/lincs-code/lincs_inputs.txt
