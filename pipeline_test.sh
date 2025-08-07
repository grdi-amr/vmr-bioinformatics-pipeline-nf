#!/bin/bash
##SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=100:00:00
#SBATCH --job-name=rgimob
#SBATCH --error=mob_db.err
#SBATCH --output=mob_db.out
##SBATCH --cpus-per-task=1
ulimit -u 16384
pid=""
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1


# Run database stub test
~/nextflow run pipeline.nf \
  --contigs "test_data/*.fa" \
  -profile singularity,test -resume --dry-run
