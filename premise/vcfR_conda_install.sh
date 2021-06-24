#!/bin/sh

module purge
module load anaconda/colsa

conda create -n vcfR r-essentials r-base
conda activate vcfR
conda install -c bioconda r-vcfr
