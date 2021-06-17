#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcf2ped"
#SBATCH --output=plink1.9.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate plink1.9

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter

srun plink --vcf test_dat.vcf --mind 0.1 --geno 0.1 --out test_dat --allow-extra-chr --make-bed 

