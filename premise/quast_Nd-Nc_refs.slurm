#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J quast
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24 
#SBATCH --output quast.out
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

srun quast -o Nd_ref_quast -m 0 -e -t 24 Nd_SPANDx_all_seqs/ref.fasta
srun quast -o Nc_ref_quast -m 0 -e -t 24 Nc_SPANDx_all_seqs/ref.fasta
