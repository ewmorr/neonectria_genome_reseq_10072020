#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcf2ped"
#SBATCH --output=plink1.9.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate plink1.9

cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/

srun plink --vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf --recode 01 --missing-genotype 9 --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9 --allow-extra-chr


