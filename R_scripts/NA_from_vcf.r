require(vcfR)
require(tidyr)
require(ggplot2)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#filtered VCF
vcf <- read.vcfR("Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
gt.na = is.na(gt)

gt.na.perc = data.frame(sample = colnames(gt.na), missing = apply(gt.na, 2, sum)/nrow(gt.na) )


p1 = ggplot(gt.na.perc, aes(x=sample, y=missing)) +
    geom_col() +
    geom_hline(yintercept = 0.26) +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 90)
    )

pdf("figures/missing_by_sample.pdf", width = 20, height = 6)
p1
dev.off()

