#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vsearch"
#SBATCH --output=vsearch_MAT_map.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_03312022/

INDIR=reads/bbtrimmed
OUTDIR=reads/bbtrimmed/vsearch

srun vsearch --usearch_global $OUTDIR/all_seqs.derep.fasta \
    --db ~/neonectria_genome_reseq_10072020/Neonectria_MAT1_and_MAT2.fasta \
    --sizein \
    --id 1 \
    --otutabout all_seqs.derep.Nf_Nd_MAT.otu_tab.txt \
    --threads 16
