#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="nextflow"
#SBATCH --output=nextflow.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node[101-112,117-125]
##https://premise.sr.unh.edu/slurm.html#srun-excludenode117node118-hostname
##excluding a bunch bc they are often used by priority users and we don't
## need the special resources

cd ~/Nc_SPANDx_all_seqs/

module purge
module load anaconda/colsa

conda activate spandx
source ~/.bashrc

nextflow run ~/SPANDx_git_clone/
