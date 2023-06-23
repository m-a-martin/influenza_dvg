#!/bin/bash

#SBATCH --job-name=simulatedSEIR_prop60_early45
#SBATCH --partition=week-long-cpu
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4



beast -threads 4 SEIR_simple_prop60_sampled_during_early45_1218.xml
