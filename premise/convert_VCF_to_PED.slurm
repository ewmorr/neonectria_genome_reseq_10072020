#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcf2ped"
#SBATCH --output=vcftools.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter

srun vcftools --vcf Nf.out.filtered.LD_filtered_0.5_10Kb.comma_rm.vcf --plink --chr 1 --out Nf.out.filtered.LD_filtered_0.5_10Kb
srun vcftools --vcf Nf.out.filtered.LD_filtered_0.5_50Kb.comma_rm.vcf --plink --chr 1 --out Nf.out.filtered.LD_filtered_0.5_50Kb

