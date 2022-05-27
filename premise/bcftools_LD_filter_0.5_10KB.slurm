#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bcftools_LD"
#SBATCH --output=bcftools_LD.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/

srun bcftools +prune -l 0.5 -w 10000 out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf -Ov -o out.filtered.LD_filtered_0.5_10Kb.vcf

