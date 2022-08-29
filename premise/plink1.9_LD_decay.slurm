#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="plinkLD"
#SBATCH --output=plink1.9.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate plink1.9

cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/

#srun plink --bfile out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.bim \
#    --r2 \
#    --ld-window-r2 0 \#include all SNP correaltions in report
#    --ld-window-r2 999999 \#maximum number of SNPs between pairwise comparisons, set to more than number of SNPs
#    --ld-window-kb 6000 \#set to greater than size of largest contig (5.5 Mb)
#    --gz \
#    --parallel 24 \
#    --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_decay


srun plink --bfile out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.bim --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 6000 --gz --parallel 24 --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_decay
