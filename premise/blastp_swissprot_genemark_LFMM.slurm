#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J blastp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24 
#SBATCH --output Nf_blastp.out

module purge
module load linuxbrew/colsa

srun blastp -task blastp -query ~/neonectria_genome_reseq_10072020/genemark_run/LFMM_query_seqs.faa \
-db ~/blast_dbs/swiss-prot_fungi/swiss-prot_fungi_reviewed.fasta \
-outfmt 6 -culling_limit 5 -num_threads 24 \
-out ~/neonectria_genome_reseq_10072020/genemark_LFMM_query_seqs-swissprot_fungi_reviewed.blastp
