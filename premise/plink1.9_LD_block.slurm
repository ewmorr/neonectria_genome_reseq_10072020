#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="plinkLD"
#SBATCH --output=plink1.9.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate plink1.9

cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/

srun plink --bfile out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode --blocks 'no-pheno-req' --blocks-max-kb 200 --threads 24 --allow-extra-chr --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_block


