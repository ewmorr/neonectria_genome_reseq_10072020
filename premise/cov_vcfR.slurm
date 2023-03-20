#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcR"
#SBATCH --output=vcfR_depth.out
#SBATCH --partition=shared

cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/

module purge
module load anaconda/colsa

conda activate vcfR

srun Rscript ~/repo/neonectria_genome_reseq_10072020/R_scripts/depth_from_vcf.r
