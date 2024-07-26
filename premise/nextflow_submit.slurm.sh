#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="nextflow"
#SBATCH --output=nextflow.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node117,node118,node105

cd ~/Nc_SPANDx_all_seqs/

module purge
module load anaconda/colsa

conda activate spandx
source ~/.bashrc

nextflow run ~/SPANDx_git_clone/
