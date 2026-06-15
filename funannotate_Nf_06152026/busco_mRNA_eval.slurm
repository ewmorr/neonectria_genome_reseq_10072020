#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24 
#SBATCH --output busco.out
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

cd ~/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026/output_funann/predict_results/


srun busco -i Neonectria_faginata.mrna-transcripts.fa -l $HOME/augustus_config/lineage/sordariomyceta_odb9 -m tran -o busco_transcript_eval_all -f -c 24
