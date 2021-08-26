#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcftools"
#SBATCH --output=vcftools.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/

####This *should* set individual genotypes to "." rather than filtering out whole variants. Can then filter based on % missing genotypes. Will need to check number of rows PASS remaining and number NA genotypes

srun vcftools --vcf out.filtered.PASS.no_comma.vcf --out out.filtered.PASS.DP_filtered --recode --minDP 4 --maxDP 22


