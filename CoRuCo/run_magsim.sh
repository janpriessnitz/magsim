#!/bin/bash
#SBATCH -J CoRuCo
#SBATCH -N 1
#SBATCH --partition=qcpu
#SBATCH -t 5:00:00
##SBATCH --mem=200G
#SBATCH --exclusive --no-requeue
#SBATCH -A OPEN-27-59
#SBATCH --ntasks-per-node=36
#SBATCH -c 1

cd $SLURM_SUBMIT_DIR

time ~/magsim/magsim/magsim config.in > stdout
