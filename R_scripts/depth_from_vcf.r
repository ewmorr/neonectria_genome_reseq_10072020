require(vcfR)
require(tidyr)
require(ggplot2)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#filtered VCF
vcf <- read.vcfR("~/SPANDx_Nf_run2/Outputs/Master_vcf/out.filtered.PASS.vcf", verbose = FALSE)
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)

dpf <- pivot_longer(data.frame(dp), names_to = "Sample", values_to = "Depth", cols = everything()) #DAMN YOU HADLEY!!!!
dpf <- dpf[ dpf$Depth > 0,]

p1 <- ggplot(dpf, aes(x=Sample, y=Depth)) +
    geom_violin(fill="#C0C0C0", adjust=1.0, scale = "count", trim=TRUE) +
    stat_summary(fun=median, geom="point", shape=23, size=2) +
    my_gg_theme +
    scale_y_continuous(trans=scales::log2_trans(),
        breaks=c(1, 10, 100, 800),
        minor_breaks=c(1:10, 2:10*10, 2:8*100)
    ) +
    labs(x = "samples") +
    theme(
        axis.text.x = element_blank()
    )

pdf("depth_by_samples.filtered.PASS.pdf", width = 12, height = 6)
p1
dev.off()

#median/mean coverage
median_cov = apply(dp, 2, median, na.rm = T)
mean_cov = apply(dp, 2, mean, na.rm = T)

mean(mean_cov)

cov_mat = matrix(c("median", median(median_cov), "mean", mean(mean_cov) ), nrow = 2)
write.table(data.frame(cov_mat), "VCF.filtered.PASS.coverage.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#unfiltered VCF
vcf <- read.vcfR("~/SPANDx_Nf_run2/Outputs/Master_vcf/out.vcf", verbose = FALSE)
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)

dpf <- pivot_longer(data.frame(dp), names_to = "Sample", values_to = "Depth", cols = everything()) #DAMN YOU HADLEY!!!!
dpf <- dpf[ dpf$Depth > 0,]

p1 <- ggplot(dpf, aes(x=Sample, y=Depth)) +
geom_violin(fill="#C0C0C0", adjust=1.0, scale = "count", trim=TRUE) +
stat_summary(fun=median, geom="point", shape=23, size=2) +
my_gg_theme +
scale_y_continuous(trans=scales::log2_trans(),
breaks=c(1, 10, 100, 800),
minor_breaks=c(1:10, 2:10*10, 2:8*100)
) +
labs(x = "samples") +
theme(
axis.text.x = element_blank()
)

pdf("depth_by_samples.unfiltered.pdf", width = 12, height = 6)
p1
dev.off()

#median/mean coverage
median_cov = apply(dp, 2, median, na.rm = T)
mean_cov = apply(dp, 2, mean, na.rm = T)

mean(mean_cov)

cov_mat = matrix(c("median", median(median_cov), "mean", mean(mean_cov) ), nrow = 2)
write.table(data.frame(cov_mat), "VCF.unfiltered.coverage.txt", quote = F, sep = "\t", row.names = F, col.names = F)
