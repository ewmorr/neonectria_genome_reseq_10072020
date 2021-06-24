#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="coverage"
#SBATCH --output=coverage.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/SPANDx_Nf_run2/Outputs/bams


for i in *.bam
do

SAMPLE=${i%.dedup.bam}

srun samtools depth -a $i > ${SAMPLE}.coverage.txt
srun perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/coverage_per_contig_samtools.pl ${SAMPLE}.coverage.txt > ${SAMPLE}.coverage_by_sequence.txt

done

