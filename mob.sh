#!/bin/bash
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --job-name=shovill
#SBATCH --error=mob_pipe.err
#SBATCH --output=mob_pipe.out

../nextflow run download_databases.nf -profile singularity --download_all -c  nextflow_slurm.config
