#!/bin/bash

#SBATCH --job-name=microbiome
#SBATCH --mail-type=ALL
#SBATCH --partition=amd
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%j.log

module load openjdk/11.0.12_7
module load singularity/3.8.5


./nextflow -log logs/nextflow.log run \
    workflows/main.nf \
    -params-file inputs.json \
    -resume

#ealloc_1f962_microbiome
