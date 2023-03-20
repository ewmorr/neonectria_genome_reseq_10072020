#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="admixture"
#SBATCH --output=admixture.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/Nd_SPANDx_all_seqs_admixture
#prev had up to 15 for Nf
for K in 1 2 3 4 5 6 7 8
do(
    srun admixture --haploid="*" --cv ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.ped $K | tee log$K.out
)
done

grep -h CV log*.out > CV_by_K.text
