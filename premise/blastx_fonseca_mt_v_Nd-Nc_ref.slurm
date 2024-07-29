#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J blastn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node117,node118
#SBATCH --output mt_blastn.out

module purge
module load linuxbrew/colsa

srun blastn -task blastn -db ~/reference_genomes/fonseca_mt/Fonseca_nectriaceae_mt_genomes.fasta \
-query ~/Nd_SPANDx_all_seqs/ref.fasta \
-evalue 0.00001 \
-outfmt 6 -culling_limit 1 -num_threads 24 \
-out ~/neonectria_genome_reseq_10072020/Nd_fonseca_mt.blast

srun blastn -task blastn -db ~/reference_genomes/fonseca_mt/Fonseca_nectriaceae_mt_genomes.fasta \
-query ~/Nc_SPANDx_all_seqs/ref.fasta \
-evalue 0.00001 \
-outfmt 6 -culling_limit 1 -num_threads 24 \
-out ~/neonectria_genome_reseq_10072020/Nc_fonseca_mt.blast
