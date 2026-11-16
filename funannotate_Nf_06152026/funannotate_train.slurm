#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J funannotate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output funannotate_train.out
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa
conda activate funannotate_new

gunzip ~/Nf_maker_reference_genome/*fasta.gz

cd ~/neonectria_faginata_draft_genome/Nf_annotate/funannotate_11162026


funannotate train -i ../funannotate_06152026/RepeatMasker/ref.fasta.masked.clean -s "Neonectria faginata" \
    -o output_funann \
    --left ../funannotate_09022026/trinity_work/trimmed/*_R1_* \
    --right ../funannotate_09022026/trinity_work/trimmed/*_R2_* \
    --trinity ../funannotate_09022026/trinity_out/trinity_out.Trinity.fasta \
    --stranded RF --jaccard_clip --no_trimmomatic \
    --cpus 24
    
gzip ~/Nf_maker_reference_genome/*fasta
