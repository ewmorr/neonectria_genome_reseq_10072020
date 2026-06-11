#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J antismash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output antismash.out
#SBATCH --exclude=node117,node118


###########################
#antismash instead run on web server
#the conda env is broken and the online server runs fast
###########################

module purge
module load anaconda/colsa

conda activate antismash-7.1.0

cd ~/neonectria_pop_gen_analyses/Nd_annotate
mkdir output_antismash

gbkFile=$HOME/neonectria_pop_gen_analyses/Nd_annotate/output_funann/predict_results/Neonectria_ditissima.gbk

srun antismash --taxon fungi \
    --output-dir output_antismash \
    -c 24 \
    $gbkFile
