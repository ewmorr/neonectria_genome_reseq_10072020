#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J makerSNAP
#SBATCH --output maker_SNAP.out
#SBATCH --ntasks=24

module purge
module load anaconda/colsa
source activate maker-3.01.02

AUGUSTUS_CONFIG_PATH=/mnt/home/garnas/ericm/augustus_config/config/
#mkdir ~/neonectria_genome_reseq_10072020/maker_run
#
cd ~/neonectria_genome_reseq_10072020/maker_run

mv ~/neonectria_genome_reseq_10072020/maker_run/maker_opts_1.ctl ~/neonectria_genome_reseq_10072020/maker_run/maker_opts.ctl

mpiexec -n 24 maker -base maker1 -f -fix_nucleotides
wait

#Combine and convert SNAP profiles
echo combining first SNAP
base=maker1

#mv maker CTL files for first run
mkdir $base.hmm
cd    $base.hmm


# maker creates an output folder with a datastore index. use gff3_merge to combine all the GFFs
srun --ntasks=1 gff3_merge -d ../$base.maker.output/$base\_master_datastore_index.log

# convert this GFF to a snap model:
srun --ntasks=1 ~/perl_scripts/maker2zff_mod.pl -n $base.all.gff
srun --ntasks=1 fathom genome.ann genome.dna -categorize 1000
srun --ntasks=1 fathom uni.ann uni.dna -export 1000 -plus
srun --ntasks=1 forge export.ann export.dna
srun --ntasks=1 hmm-assembler.pl $base . > $base.snap.hmm

#Move maker CTL files for second run
cd ~/neonectria_genome_reseq_10072020/maker_run
mv ~/neonectria_genome_reseq_10072020/maker_run/maker_opts.ctl ~/neonectria_genome_reseq_10072020/maker_run/maker_opts_1.ctl
mv ~/neonectria_genome_reseq_10072020/maker_run/maker_opts_2.ctl ~/neonectria_genome_reseq_10072020/maker_run/maker_opts.ctl



#Second maker run
echo second maker run
mpiexec -n 24 maker -base maker2 -f -fix_nucleotides
wait

#combine second SNAP run
echo combining second SNAP run

base=maker2

mkdir $base.hmm
cd    $base.hmm


# maker creates an output folder with a datastore index. use gff3_merge to combine all the GFFs
srun --ntasks=1 gff3_merge -d ../$base.maker.output/$base\_master_datastore_index.log

# convert this GFF to a snap model:
srun --ntasks=1 ~/perl_scripts/maker2zff_mod.pl -n $base.all.gff
srun --ntasks=1 fathom genome.ann genome.dna -categorize 1000
srun --ntasks=1 fathom uni.ann uni.dna -export 1000 -plus
srun --ntasks=1 forge export.ann export.dna
srun --ntasks=1 hmm-assembler.pl $base . > $base.snap.hmm


#Move maker CTL files for final run
cd ~/neonectria_genome_reseq_10072020/maker_run
mv ~/neonectria_genome_reseq_10072020/maker_run/maker_opts.ctl ~/neonectria_genome_reseq_10072020/maker_run/maker_opts_2.ctl
mv ~/neonectria_genome_reseq_10072020/maker_run/maker_opts_3.ctl ~/neonectria_genome_reseq_10072020/maker_run/maker_opts.ctl

#final maker run
echo do final maker run
mpiexec -n 24 maker -base maker3 -f -fix_nucleotides
wait

mv ~/neonectria_genome_reseq_10072020/maker_run/maker_opts.ctl ~/neonectria_genome_reseq_10072020/maker_run/maker_opts_3.ctl

