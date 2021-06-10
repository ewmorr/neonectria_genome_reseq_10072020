#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bcftools_LD"
#SBATCH --output=bcftools_LD.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/

srun bcftools +prune -l 0.5 -w 50000 ~/SPANDx_Nf/Outputs/Master_vcf/out.filtered.vcf -Ov -o Nf.out.filtered.LD_filtered_0.5_50Kb.vcf

