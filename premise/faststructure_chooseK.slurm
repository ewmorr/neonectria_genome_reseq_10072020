#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="faststructure"
#SBATCH --output=faststructure.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate faststructure

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/chooseK.py \
--input=faststructure/simple > faststructure/simple.chooseK.txt


