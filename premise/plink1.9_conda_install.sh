#!/bin/sh

module purge
module load anaconda/colsa

conda create --name plink1.9 -c bioconda plink

conda activate
