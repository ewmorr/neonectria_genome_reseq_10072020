#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J blastn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node117,node118
#SBATCH --output ID_blastn.out

module purge
module load linuxbrew/colsa

for i in 1 2 3 4 5
do(
cd ~/neonectria_coccinea_genome_reseq_07182024/Neco-$i

formatdb -i contigs.fasta -p F -o T
blastn -db contigs.fasta -query ../Nc_ID.fasta -out ID.blast \
    -evalue 0.00001 -outfmt 6 -culling_limit 1 -num_threads 24
    
)
done
