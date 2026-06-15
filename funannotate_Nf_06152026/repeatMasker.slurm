#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J repMask
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --output repeatMasker.out
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

cd ~/neonectria_pop_gen_analyses/Nd_annotate


srun RepeatMasker -dir . -xsmall -pa 12 -lib $HOME/RepeatMaskerLib/fngrep.ref ref.fasta
