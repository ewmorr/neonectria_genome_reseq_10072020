#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="VCFvalidate"
#SBATCH --output=vcf_validator.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

echo "file one"
srun vcf-validator neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf
echo "file two"
srun srun vcf-validator neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_50Kb.vcf

