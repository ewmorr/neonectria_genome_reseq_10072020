#!/bin/sh

module purge
module load anaconda/colsa

#initial setup to add bioconda
#conda config --add channels defaults
#conda config --add channels conda-forge
#conda config --add channels bioconda

conda create --name spandx -c dsarov -c bioconda -c conda-forge spandx

conda activate spandx

nextflow pull dsarov/spandx

nextflow clone dsarov/spandx SPANDx_git_clone

#navigate to SPANDx_nextflow/nextflow.config and update the `executor` line to "slurm" OR run with --executor "slurm" flag

#The current install of anaconda on Premise (02/16/2021) does not allow calling conda activate from subshells (e.g., within nextflow).
#Add the following lines to .bashrc **outside** of the # >>>conda init <<< block https://github.com/conda/conda/issues/7980 and https://stackoverflow.com/questions/69001097/conda-4-10-3-and-snakemake-5-conda-exe-problem
export -f conda
export -f __conda_activate
export -f __conda_reactivate
export -f __conda_hashr
export -f __conda_exe

