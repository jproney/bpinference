#!/bin/bash
#logistic_benchmark.sbatch
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 2000

module load R/3.5.1
Rscript "/michorlab/jroney/bpinference/benchmarks/double_logistic_benchmark.R" "${nsims}"

