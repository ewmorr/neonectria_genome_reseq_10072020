#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vsearch"
#SBATCH --output=vsearch_ITS_map.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/

INDIR=reads/bbtrimmed
OUTDIR=reads/bbtrimmed/vsearch

srun vsearch --usearch_global $OUTDIR/al_seqs.derep.fasta \
    --db Nf_Nd_ITS2.fasta \
    --sizein \
    --id 0.97 \
    --otutabout all_seqs.derep.Nf_Nd_ITS.otu_tab.txt \
    --threads 16
