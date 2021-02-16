#!/bin/sh

module purge
module load anaconda/colsa

#initial setup to add bioconda
#conda config --add channels defaults
#conda config --add channels conda-forge
#conda config --add channels bioconda

conda create --name spandx -c dsarov -c bioconda -c conda-forge spandx

conda activate spandx


nextflow clone dsarov/spandx SPANDx_nextflow

#navigate to SPANDx_nextflow/nextflow.config and update the `executor` line to "slurm" OR run with --executor "slurm" flag
#then run
nextflow pull dsarov/spandx

#ACTUALLY this doesn't seem to work as recommended, so try running with the --executor flag
#executor flag works

#The current install of anaconda on Premise (02/16/2021) does not allow calling conda activate from subshells (e.g., within nextflow).
#Add the following lines to .bashrc **outside** of the # >>>conda init <<< block https://github.com/conda/conda/issues/7980
export -f conda
export -f __conda_activate
export -f __conda_reactivate
export -f __conda_hashr

