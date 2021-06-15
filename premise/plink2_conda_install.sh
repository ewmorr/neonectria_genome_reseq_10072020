#!/bin/sh

module purge
module load anaconda/colsa

conda create --name plink2 -c bioconda plink2

conda activate
