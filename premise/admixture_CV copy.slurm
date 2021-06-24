#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcR"
#SBATCH --output=vcfR_depth.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

cd ~/SPANDx_Nf_run2/Outputs/Master_vcf/

module/purge
module load anaconda/colsa

srun Rscript ~/repo/neonectria_genome_reseq_10072020/R_scripts/depth_from_vcf.r
