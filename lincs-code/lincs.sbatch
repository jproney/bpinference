#!/bin/bash
##lincs.sbatch
##SBATCH -n 1
##SBATCH -N 1
##SBATCH --mem 2000


module load R/3.5.1
Rscript "/michorlab/jroney/bpinference/lincs-code/lincs_logistic_birth_parallel.R" "${cell_name}" "${drug}"
