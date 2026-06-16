#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J funannotate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output funannotate_annotate.out
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa
conda activate funannotate_new

cd ~/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026

funannDir=$HOME/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026/output_funann

eggnogPath=$HOME/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026/output_eggnog_prot/Nf_.emapper.annotations
iprscanPath=$HOME/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026/output_iprscan_prot/ipr_out.xml
# antismash is gbk format from web server
antismashPath=$HOME/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026/output_antismash/Neonectria_faginata.gbk
# phobius is short format from web server
phobiusPath=$HOME/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026/phobius_prediction.txt
# Signal P is eukarya short format from web server
signalpPath=$HOME/neonectria_faginata_draft_genome/Nf_annotate/funannotate_06152026/signalP.all.txt

srun funannotate annotate -i $funannDir \
    -o $funannDir \
    --rename RS324p \
    --cpus 24 \
    --eggnog $eggnogPath \
    --iprscan $iprscanPath \
    --antismash $antismashPath \
    --phobius $phobiusPath \
    --signalp $signalpPath
