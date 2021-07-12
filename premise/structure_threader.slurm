#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="str_th"
#SBATCH --output=structure_threader.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

export PATH=~/.local/bin:$PATH

cd ~/neonectria_genome_reseq_10072020/

structure_threader run -K 15 -R 5 \
-i ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.recode.strct_in \
-o ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/structure_th/ \
--params ~/neonectria_genome_reseq_10072020/ \
-t 24 -st /mnt/home/garnas/ericm/.local/bin/structure

