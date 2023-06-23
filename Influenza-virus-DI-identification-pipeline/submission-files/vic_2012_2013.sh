#!/bin/bash

#SBATCH --job-name=vic_2012_2013
#SBATCH --partition=day-long-cpu
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=128G


nextflow run -c Influenza-virus-DI-identification-pipeline/new-config-files/vic_2012_2013.conf  Influenza-virus-DI-identification-pipeline/full-pipeline-customizable-mm-final.nf
