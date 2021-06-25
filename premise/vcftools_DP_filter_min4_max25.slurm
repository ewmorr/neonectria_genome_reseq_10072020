#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcftools"
#SBATCH --output=vcftools.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/SPANDx_Nf/Outputs/Master_vcf/

####This *should* set individual genotypes to "." rather than filtering out whole variants. Can then filter based on % missing genotypes. Will need to check number of rows PASS remaining and number NA genotypes

srun vcftools --vcf out.filtered.PASS.vcf --out out.filtered.PASS.DP_filtered --minDP 4 --maxDP 22


