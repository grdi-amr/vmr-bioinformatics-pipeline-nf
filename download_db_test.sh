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

# Set ETE3 env var
export ETE3_HOME="$PWD/ete3_data"
mkdir -p $ETE3_HOME

# Run database stub test
~/nextflow run download_databases.nf \
  --download_all true \
  --overwrite true \
  -profile singularity,test 
