#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcf2ped"
#SBATCH --output=plink1.9.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate plink1.9

cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/

srun plink --vcf out.filtered.LD_filtered_0.5_10Kb.vcf --recode12 --out out.filtered.LD_filtered_0.5_10Kb --allow-extra-chr


