#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="faststructure"
#SBATCH --output=faststructure2.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate faststructure

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/

#for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
for K in 3 4 5 6 7 8 9 10 11 12 13 14 15
do(
    ~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K $K \
        --input=Nf.out.filtered.LD_filtered_0.5_10Kb \
        --output=faststructure/simple \
        --full --seed=100 &
)
done
wait
