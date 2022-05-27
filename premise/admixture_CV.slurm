#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="admixture"
#SBATCH --output=admixture.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/Nf_SPANDx_all_seqs_admixture

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do(
    srun admixture --haploid="*" --cv ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.ped $K | tee log$K.out
)
done

grep -h CV log*.out > CV_by_K.text
