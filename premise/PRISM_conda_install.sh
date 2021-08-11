#!/bin/sh

module purge
module load anaconda/colsa

conda create -n PRISM r-essentials r-base
conda activate PRISM
#conda install -c bioconda r-prism

R --slave -e "install.packages('prism', repos='http://cran.us.r-project.org')"
