#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J interproscan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output interproscan.out
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa

conda activate interproscan-5.69.101

cd ~/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026
mkdir output_iprscan_prot

srun interproscan.sh -i output_funann/predict_results/Neonectria_faginata.proteins.fa \
    --cpu 24 \
    -f XML \
    -o output_iprscan_prot/ipr_out.xml
