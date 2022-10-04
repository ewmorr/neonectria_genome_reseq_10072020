#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="LDplots"
#SBATCH --output=LD_plots.out
#SBATCH --partition=shared
#SBATCH --mem=512G

module purge
module load anaconda/colsa

conda activate PRISM

cd ~/neonectria_genome_reseq_10072020
srun Rscript ~/repo/neonectria_genome_reseq_10072020/R_scripts/LD_decay_plots.R
