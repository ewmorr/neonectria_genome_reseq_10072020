#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="nextflow"
#SBATCH --output=nextflow.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node101,node102,node103,node104,node105,node106,node107,node108,node109,node110,node111,node112,node117,node118,node119,node120,node121,node122,node123,node124,node125
##SBATCH --exclude=node\[101-112,117-125\]
##https://premise.sr.unh.edu/slurm.html#srun-excludenode117node118-hostname
##excluding a bunch bc they are often used by priority users and we don't
## need the special resources
## Except the node range syntax isn't working so we will list it out

cd ~/Nc_SPANDx_all_seqs/

module purge
module load anaconda/colsa

conda activate spandx
source ~/.bashrc

nextflow run ~/SPANDx_git_clone/
