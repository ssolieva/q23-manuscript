#!/bin/bash
#SBATCH --job-name=CH01cf
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --export=ALL
#SBATCH --gres=gpu:1
#SBATCH --mem=40Gb
#SBATCH --time=7-00:00:00 # 7 day limit. 

singularity run --nv \
  -B /home/ssolieva/.cache:/cache -B $(pwd):/work \
  /wistar/kulp/users/cagostino/software/colabfold_1.5.3-cuda11.8.0.sif \
  colabfold_batch . . --templates --custom-template-path ./cass_path --num-seeds=1

