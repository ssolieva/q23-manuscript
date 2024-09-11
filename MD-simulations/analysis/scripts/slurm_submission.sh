#!/bin/bash

# queue
#SBATCH --partition=amdcpu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH -t 0-12:00
#SBATCH --output=outputSlurm/%j.out

python axe_classification_every10frames.py $1 $2 $3 $4
