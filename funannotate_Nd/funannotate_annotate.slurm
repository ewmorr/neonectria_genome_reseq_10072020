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

cd ~/Nf_annotate/funannotate_run/

funannDir=/mnt/home/garnas/ewj4/Nf_annotate/funannotate_run/output_funann

eggnogPath=$HOME/Nf_annotate/funannotate_run/output_eggnog_prot/Nf_.emapper.annotations
iprscanPath=$HOME/Nf_annotate/funannotate_run/output_iprscan_prot/ipr_out.xml
antismashPath=$HOME/Nf_annotate/funannotate_run/output_antismash/Neonectria_faginata.gbk
phobiusPath=$HOME/Nf_annotate/funannotate_run/phobius_prediction.txt
signalpPath=$HOME/Nf_annotate/funannotate_run/signalP.all.txt

srun funannotate annotate -i $funannDir \
    -o $funannDir \
    --rename Q2T73 \
    --cpus 24 \
    --eggnog $eggnogPath \
    --iprscan $iprscanPath \
    --antismash $antismashPath \
    --phobius $phobiusPath \
    --signalp $signalpPath \
    --sbt template.sbt
    
    
    

