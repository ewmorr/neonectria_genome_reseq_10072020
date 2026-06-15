#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24 
#SBATCH --output busco.out
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

cd ~/neonectria_pop_gen_analyses/Nd_annotate/output_funann/predict_results/


srun busco -i Neonectria_ditissima.mrna-transcripts.fa -l $HOME/augustus_config/lineage/sordariomyceta_odb9 -m tran -o busco_transcript_eval_all -f -c 24
