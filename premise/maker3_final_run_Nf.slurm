#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J makerSNAP
#SBATCH --output maker_FINAL.out
#SBATCH --ntasks=24

module purge
module load anaconda/colsa
source activate maker-3.01.02

AUGUSTUS_CONFIG_PATH=/mnt/home/garnas/ericm/augustus_config/config/
#mkdir ~/neonectria_genome_reseq_10072020/maker_run
#
cd ~/neonectria_genome_reseq_10072020/maker_run

mv ~/neonectria_genome_reseq_10072020/maker_run/maker_opts_final.ctl ~/neonectria_genome_reseq_10072020/maker_run/maker_opts.ctl

mpiexec -n 24 maker -base makerFINAL -f -fix_nucleotides
wait

mv ~/neonectria_genome_reseq_10072020/maker_run/maker_opts.ctl ~/neonectria_genome_reseq_10072020/maker_run/maker_opts_final.ctl

base=makerFINAL

# maker creates an output folder with a datastore index. use gff3_merge to combine all the GFFs
srun --ntasks=1 gff3_merge -d $base.maker.output/$base\_master_datastore_index.log
srun --ntasks=1 fasta_merge -d $base.maker.output/$base\_master_datastore_index.log

#Extract AED<1
gff=$base.all.gff
#get IDs
perl -lne 'print $1 if /\tmRNA\t.+ID=([^;]+).+_AED=(.+?);/ and $2<1' $gff > $base.aed-1.0.ids
#extract from GFF
srun --ntasks=1 ~/perl_scripts/fgrep.pl -f $base.aed-1.0.ids $gff -d "=|;|,|:" > $base.aed-1.0.gff

#proteins
fasta=$base.proteins.fasta # or $base.transcripts.fasta
srun --ntasks=1 ~/perl_scripts/fastaqual_select.pl -f $fasta -inc $base.aed-1.0.ids > $base.proteins.aed-1.0.fasta

#transcripts
fasta=$base.transcripts.fasta 
srun --ntasks=1 ~/perl_scripts/fastaqual_select.pl -f $fasta -inc $base.aed-1.0.ids > $base.transcripts.aed-1.0.fasta
