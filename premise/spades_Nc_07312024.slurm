#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J spades
#SBATCH --output spades.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=512000

module purge
module load linuxbrew/colsa

cd ~/neonectria_coccinea_genome_reseq_07182024

indir=reads/bbtrimmed

for i in $indir/*R1*fastq.gz
do(
    fwd=${i##*/}
    pre=${fwd%R1*}
    post=${fwd##*R1}
    rev=${pre}R2${post}
    sample_id=${fwd%_S*_L00*}
    
    echo "Starting $sample_id"
    mkdir $sample_id

    spades.py -1 $indir/$fwd -2 $indir/$rev -o $sample_id -m 512 -t 24
)
done
