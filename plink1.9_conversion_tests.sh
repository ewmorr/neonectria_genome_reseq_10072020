#!/bin/bash
#From local computer
conda activate plink1.9
cd ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/plink_conversion_test

plink --vcf test_dat.vcf --recode12 --mind 0.1 --geno 0.1 --allow-extra-chr --out test_dat.recode12
plink --vcf test_dat.vcf --recode --mind 0.1 --geno 0.1 --allow-extra-chr --out test_dat.recode

plink --vcf test_dat.vcf --make-bed --mind 0.1 --geno 0.1 --allow-extra-chr --out test_dat.make-bed

plink --file test_dat.recode12 --make-bed --mind 0.1 --geno 0.1 --allow-extra-chr --out test_dat.recode12.make-bed
plink --file test_dat.recode --make-bed --mind 0.1 --geno 0.1 --allow-extra-chr --out test_dat.recode.make-bed

plink --bfile test_dat.recode12.make-bed --recode12 --mind 0.1 --geno 0.1 --allow-extra-chr --out test_dat.recode12.make-bed.recode12
plink --bfile test_dat.recode.make-bed --recode --mind 0.1 --geno 0.1 --allow-extra-chr --out test_dat.recode.make-bed.recode

plink --bfile test_dat.make-bed --recode vcf --out test_dat.make-bed --allow-extra-chr
