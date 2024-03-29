library(vcfR)

vcf <- read.vcfR("data/Nd_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.vcf", verbose = FALSE)
#test.vcf = read.vcfR("Nf_post_SPANDx/LD_filter/plink_conversion_test/test_dat.vcf", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', convertNA=TRUE, as.numeric = T)

head(gt)

nrow(gt)

colSums(is.na(gt) )/nrow(gt)
sum(rowSums(is.na(gt) )/ncol(gt) > 0.5)




par(mar=c(8,4,1,1))
#boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", log='y', las=2)
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", las=2)
abline(h=seq(0,1e4, by=100), col="#C0C0C088")
par(mar=c(5,4,4,2))


require(reshape2)
require(ggplot2)
source("~/ggplot_theme.txt")

dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
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

pdf("figures/depth_by_samples.pdf", width = 12, height = 6)
p1
dev.off()

#number of SNPs
nrow(dp)

#median coverage
median_cov = apply(dp, 2, median, na.rm = T)
mean_cov = apply(dp, 2, median, na.rm = T)
