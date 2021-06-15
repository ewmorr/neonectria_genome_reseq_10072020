#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="admixture"
#SBATCH --output=admixture.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/admixture

for K in 1 2 3 4 5 6 7 8 9 10
do(
    srun admixture --cv ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.lt25missing.ped $K | tee log$K.out
)
done

grep -h CV log*.out > CV_by_K.text
