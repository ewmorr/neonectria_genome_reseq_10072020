#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vcftools"
#SBATCH --output=vcftools.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/


while IFS= read -r line; do
    site1=$(echo $line | cut -f1 -d' ')
    site2=$(echo $line | cut -f2 -d' ')

    srun vcftools --vcf Nf.out.filtered.LD_filtered_0.5_10Kb.vcf \
        --weir-fst-pop ~/neonectria_genome_reseq_10072020/pairwise_site_IDs/${site1}.txt \
        --weir-fst-pop ~/neonectria_genome_reseq_10072020/pairwise_site_IDs/${site2}.txt \
        --out pairwise_fisher/${site1}_vs_${site2}

done < ~/neonectria_genome_reseq_10072020/pairwise_site_IDs/pairwise_site_combs.txt
