#!/bin/bash

#SBATCH --job-name=perth_2010_2012
#SBATCH --partition=week-long-cpu
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=128G


nextflow run -c Influenza-virus-DI-identification-pipeline/new-config-files/perth_2010_2012.conf  Influenza-virus-DI-identification-pipeline/full-pipeline-customizable-mm-final.nf -resume

