#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vsearch"
#SBATCH --output=vsearch_ITS_map.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_09182023/

INDIR=reads/bbtrimmed
OUTDIR=reads/bbtrimmed/vsearch
#mkdir $OUTDIR

#for i in reads/*R1*.fastq.gz
#do

#    FILE=${i##*/}
#    BEFFILE=${FILE%R1*}
#    AFTFILE=${FILE##*R1}
#    FORWARD=$FILE
#    REVERSE=${BEFFILE}R2${AFTFILE}
#    sample_id=${FILE%_S*_L001*}

#    srun vsearch --fastq_mergepairs $INDIR/$FORWARD \
        --threads 16 \
        --reverse $INDIR/$REVERSE \
        --fastaout $OUTDIR/${sample_id}.fasta

#    srun vsearch --derep_fulllength $OUTDIR/${sample_id}.fasta \
        --sizeout \
        --threads 16 \
        --output $OUTDIR/${sample_id}.derep.fasta

#    sed -i "s/^>.*/&;sample=${sample_id};/g" $OUTDIR/${sample_id}.derep.fasta
#done

#cat $OUTDIR/*.derep.fasta > $OUTDIR/all_seqs.derep.fasta

srun vsearch --usearch_global $OUTDIR/all_seqs.derep.fasta \
    --db ~/neonectria_genome_reseq_10072020/Nf_Nd_ITS2.fasta \
    --sizein \
    --id 0.97 \
    --otutabout all_seqs.derep.Nf_Nd_ITS.otu_tab.txt \
    --threads 16
