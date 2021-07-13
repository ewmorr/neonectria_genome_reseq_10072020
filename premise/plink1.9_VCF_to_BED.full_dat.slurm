#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcf2bed"
#SBATCH --output=plink1.9.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate plink1.9

cd ~/SPANDx_Nf/Outputs/Master_vcf

srun plink --vcf out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf --out out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode --allow-extra-chr --make-bed
