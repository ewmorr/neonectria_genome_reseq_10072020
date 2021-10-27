require(PopGenome)
require(dplyr)
require(ggplot2)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")


#READING DATA and checking it
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)

scfs = c("tig00000001_pilon", "tig00000008_pilon", "tig00000025_pilon", "tig00000071_pilon", "tig00000091_pilon", "tig00000129_pilon", "tig00000138_pilon", "tig00000233_pilon", "tig00000235_pilon", "tig00000257_pilon", "tig00000314_pilon", "tig00000320_pilon", "tig00000350_pilon", "tig00000778_pilon", "tig00000779_pilon", "tig00000780_pilon", "tig00007946_pilon", "tig00007950_pilon", "tig00007951_pilon", "tig00007953_pilon")


snp <- readVCF("Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf.gz", include.unknown = T, approx = T, frompos = 1, topos = 10, tid = scfs)

#CAN ONLY READ ONE SCF / REGION AT A TIME. SHOULD JUST USE READDATA FUNCTION TO READ IN EACH VCF SEPARATELY

