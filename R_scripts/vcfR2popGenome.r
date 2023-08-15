library(vcfR)


vcf <- read.vcfR("data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf.gz", verbose = FALSE)
#test.vcf = read.vcfR("Nf_post_SPANDx/LD_filter/plink_conversion_test/test_dat.vcf", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', convertNA=TRUE, as.numeric = T)

head(gt)

nrow(gt)
#128629 SNPs
sum(apply(gt, 1, function(x) any(is.na(x))))
#128624 contain NA

vcfR::is.biallelic(vcf) %>% sum
#12689 (all biallelic)


