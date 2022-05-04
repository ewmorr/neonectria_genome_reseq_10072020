#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24 
#SBATCH --output busco.out

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/maker2_run/


srun busco -i makerFINAL.transcripts.aed-1.0.fasta -l $HOME/augustus_config/lineage/sordariomyceta_odb9 -m tran -o busco_transcript_eval -f -c 24
