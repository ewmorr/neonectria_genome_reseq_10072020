#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bcftools_LD"
#SBATCH --output=bcftools_LD.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter

srun bcftools +prune -l 0.5 -w 10000 ~/SPANDx_Nf/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.vcf -Ov -o Nf.out.filtered.LD_filtered_0.5_10Kb.vcf

