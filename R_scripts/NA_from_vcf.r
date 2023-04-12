require(vcfR)
require(tidyr)
require(ggplot2)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#filtered VCF
#vcf <- read.vcfR("~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.recode.vcf", verbose = FALSE)
vcf <- read.vcfR("data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.recode.vcf", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
gt.na = is.na(gt)

gt.na.perc = data.frame(sample = colnames(gt.na), missing = apply(gt.na, 2, sum)/nrow(gt.na) )


p1 = ggplot(gt.na.perc, aes(x=sample, y=missing)) +
    geom_col() +
    geom_hline(yintercept = 0.35) +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 90)
    )
p1

quantile(gt.na.perc$missing, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))

quantile(gt.na.perc$missing, probs = 0.9)[1]
             
p1 = ggplot(gt.na.perc, aes(x=sample, y=missing)) +
    geom_col() +
    geom_hline(yintercept = quantile(gt.na.perc$missing, probs = 0.9)[1]) +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 90)
    )
p1

pdf("figures/Nf_missing_by_sample.new_assembly.pdf", width = 24, height = 6)
p1
dev.off()

