#!/bin/sh

module purge
module load anaconda/colsa

#initial setup to add bioconda
#conda config --add channels defaults
#conda config --add channels conda-forge
#conda config --add channels bioconda

#mamba create --name spandx -c dsarov -c bioconda -c conda-forge spandx

mamba create --name spandx -c dsarov -c bioconda -c conda-forge spandx "python>=3" "nextflow<=22.10.6-0"

conda activate spandx

nextflow pull dsarov/spandx

nextflow clone dsarov/spandx SPANDx_git_clone

#navigate to SPANDx_nextflow/nextflow.config and update the `executor` line to "slurm" OR run with --executor "slurm" flag

#The current install of anaconda on Premise (02/16/2021) does not allow calling conda activate from subshells (e.g., within nextflow).
#Add the following lines to .bashrc **outside** of the # >>>conda init <<< block https://github.com/conda/conda/issues/7980 and https://stackoverflow.com/questions/69001097/conda-4-10-3-and-snakemake-5-conda-exe-problem
```
export -f conda
export -f __conda_activate
export -f __conda_reactivate
export -f __conda_hashr
export -f __conda_exe
```
Make sure that nextflow.config is updated if necessary (https://github.com/dsarov/SPANDx#usage)
The config file is where CPUs etc are denoted as well as the resource manager (e.g., SLURM). Newer versions of the package also have `notrim` set to `true`. This should be `false` in that case. Also, note that we have cloned the git repo after installing via conda, and made some modificaations to the `main.nf` script (i.e., changing the `gatk HaplotypeCaller` comand at line 865 to include `--ploidy 1` flag) and to the `./bin/Master_vcf.sh` script (i.e., removing `-ploidy 1` from the `gatk GenotypeGVCFs` command). Also, note that `.bashrc` may need to be updated as described [here](./SPANDx_conda_install.sh). Finally, the reference genome assembly must be loacted in the SPANDx working directory (along with the reads), and can be indicated by path in the config file. Then, to run SPANDx (note that nextflow is pointed to the cloned git repo)
