#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bcftools"
#SBATCH --output=bcftools.out
#SBATCH --partition=shared
#SBATCH --exclude=node[117-118]

module purge
module load linuxbrew/colsa

cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/

srun bcftools view -i 'F_MISSING<0.25' out.filtered.PASS.DP_filtered.recode.vcf -Ov -o out.filtered.PASS.DP_filtered.lt25missing.vcf

