#!/bin/bash

#SBATCH --job-name=hk_2014_2015
#SBATCH --partition=week-long-cpu
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=128G

nextflow run -c Influenza-virus-DI-identification-pipeline/new-config-files/hk_2014_2015.conf  Influenza-virus-DI-identification-pipeline/full-pipeline-customizable-mm-final.nf -resume
