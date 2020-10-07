#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bbduk"
#SBATCH --output=trimmomatic.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/

INDIR=reads
OUTDIR=reads/trimmomatictrimmed
mkdir $OUTDIR

ADAPTERS=/mnt/lustre/software/linuxbrew/colsa/Cellar/trimmomatic/0.39/share/trimmomatic/adapters/TruSeq3-PE.fa

for i in reads/*R1*.fastq.gz
do

    FILE=${i##*/}
    BEFFILE=${FILE%R1*}
    AFTFILE=${FILE##*R1}
    FORWARD=$FILE
    REVERSE=${BEFFILE}R2${AFTFILE}

    srun trimmomatic PE -threads 24 $INDIR/$FORWARD $INDIR/$REVERSE\
        $OUTDIR/$FORWARD $OUTDIR/$FORWARD\
        $OUTDIR/$REVERSE $OUTDIR/$REVERSE\
        ILLUMINACLIP:$ADAPTERS:2:30:10\
        SLIDINGWINDOW:4:15 MINLEN:36

done
