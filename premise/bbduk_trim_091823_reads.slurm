#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bbduk"
#SBATCH --output=bbduk.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate qiime2-2019.4

cd ~/neonectria_genome_reseq_09182023/

INDIR=reads
OUTDIR=reads/bbtrimmed
mkdir $OUTDIR

ADAPTERS=/mnt/lustre/software/anaconda/colsa/envs/qiime2-2019.4/opt/bbmap-38.57-0/resources/adapters.fa

for i in reads/*R1*.fastq.gz
do

    FILE=${i##*/}
    BEFFILE=${FILE%R1*}
    AFTFILE=${FILE##*R1}
    FORWARD=$FILE
    REVERSE=${BEFFILE}R2${AFTFILE}

srun bbduk.sh in=$INDIR/$FORWARD in2=$INDIR/$REVERSE \
out=$OUTDIR/$FORWARD out2=$OUTDIR/$REVERSE \
ref=$ADAPTERS ktrim=r k=23 mink=11 hdist=1 tpe tbo \
qtrim=r trimq=10 -minlength=36

done
