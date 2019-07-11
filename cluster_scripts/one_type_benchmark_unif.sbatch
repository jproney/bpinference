#!/bin/bash
#one_type_benchmark.sbatch
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 2000

module load R/3.5.1
Rscript "/michorlab/jroney/bpinference/benchmarks/one_type_benchmark_unif.R" "${nsims}"

