#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcf2ped"
#SBATCH --output=plink1.9.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate plink1.9

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter

srun plink --vcf Nf.out.filtered.LD_filtered_0.5_10Kb.lt25missing.vcf --recode12 --out Nf.out.filtered.LD_filtered_0.5_10Kb.lt25missing --allow-extra-chr --make-bed --[output-]missing-genotype

srun plink --vcf Nf.out.filtered.LD_filtered_0.5_50Kb.lt25missing.vcf --recode12 --out Nf.out.filtered.LD_filtered_0.5_50Kb.lt25missing --allow-extra-chr --make-bed --[output-]missing-genotype

