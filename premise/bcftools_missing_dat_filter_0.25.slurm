#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bcftools"
#SBATCH --output=bcftools.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter

srun bcftools view -i 'F_MISSING<0.25' Nf.out.filtered.LD_filtered_0.5_10Kb.vcf -Ov -o Nf.out.filtered.LD_filtered_0.5_10Kb.lt25missing.vcf

srun bcftools view -i 'F_MISSING<0.25' Nf.out.filtered.LD_filtered_0.5_50Kb.vcf -Ov -o Nf.out.filtered.LD_filtered_0.5_50Kb.lt25missing.vcf
