require(pcadapt)
require(ggplot2)
require(dplyr)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")
require(RColorBrewer)
require(gridExtra)

makeBed <- "Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.bed"
makeBed.file <- read.pcadapt(makeBed, type = "bed")
bed_mat.makeBed = pcadapt::bed2matrix(makeBed.file)

sum(bed_mat.makeBed[is.na(bed_mat.makeBed) == F] == 0)
sum(bed_mat.makeBed[is.na(bed_mat.makeBed) == F] == 1)
sum(bed_mat.makeBed[is.na(bed_mat.makeBed) == F] == 2)

x <- pcadapt(input = makeBed.file, K = 20, ploidy = 2) #K must be less than the number of species
#Ploidy is set to 2 because the function throws an error otherwise. This looks like because the conversion from BED is not recognizing haploid data. A look at `bed_mat.makeBed` shows SNPs encoded as 0 or 2 with no outcomes of 1 (i.e., het gentypes).

x_K30 <- pcadapt(input = makeBed.file, K = 30, ploidy = 2) #K must be less than the number of species
x_K50 <- pcadapt(input = makeBed.file, K = 50, ploidy = 2) #K must be less than the number of species

#plot(x, option = "screeplot")
#plot(x_K30, option = "screeplot")
plot(x_K50, option = "screeplot")

plot(x, option = "scores")
plot(x, option = "manhattan") #This is useful for full dataset
plot(x, option = "qqplot") #Another representation of p-value for SNP
plot(x, option = "stat.distribution") #chi.sq distribution

#Extract data for pretty plotting

#Plot scree plot (can use any of the different K val PCAs
#The singular values vector within PCA object contains the square root of porportion variance explained
scree_dat = data.frame(prop.var = x_K50$singular.values^2, PC = 1:length(x_K50$singular.values))

p1 = ggplot(scree_dat, aes(PC, prop.var)) +
geom_point() +
geom_line() +
scale_x_continuous(breaks = c(1,10,20,30,40,50)) +
labs(x = "PC axis", y = "Proportion variance explained") +
my_gg_theme

pdf("figures/Nf.NOT_LD_filtered.PCA_scree_plot.pdf", width = 6, height = 5)
p1
dev.off()

#Plot PCA scores
fam_info = read.table("Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.fam", header = F)
sample_metadata = read.table("sample_metadata/sample_metadata.txt", header = T)

pc_scores.metadata = data.frame(
    PC1 = x$scores[,1],
    PC2 = x$scores[,2],
    PC3 = x$scores[,3],
    PC4 = x$scores[,4],
    PC5 = x$scores[,5],
    PC6 = x$scores[,6],
    PC7 = x$scores[,7],
    PC8 = x$scores[,8],
    sample = fam_info[,1]
)

pc_scores.metadata = left_join(pc_scores.metadata, sample_metadata, by = "sample")

get.color.palette <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
state.colors = get.color.palette("Paired")

p1 = ggplot(pc_scores.metadata, aes(PC1, PC2, color = State)) +
geom_point(size = 3) +
scale_color_manual(values = state.colors, guide = "none") +
my_gg_theme

p2 = ggplot(pc_scores.metadata, aes(PC1, PC3, color = State)) +
geom_point(size = 3) +
scale_color_manual(values = state.colors, guide = "none") +
my_gg_theme

p3 = ggplot(pc_scores.metadata, aes(PC2, PC3, color = State)) +
geom_point(size = 3) +
scale_color_manual(values = state.colors) +
my_gg_theme

p4 = ggplot(pc_scores.metadata, aes(PC1, PC4, color = State)) +
geom_point(size = 3) +
scale_color_manual(values = state.colors, guide = "none") +
my_gg_theme

p5 = ggplot(pc_scores.metadata, aes(PC2, PC4, color = State)) +
geom_point(size = 3) +
scale_color_manual(values = state.colors, guide = "none") +
my_gg_theme

p6 = ggplot(pc_scores.metadata, aes(PC3, PC4, color = State)) +
geom_point(size = 3, alpha = 0.7) +
scale_color_manual(values = state.colors) +
my_gg_theme

pdf("figures/Nf.NOT_LD_filtered.PCA_plot.pdf", width = 16, height = 8)
grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2, widths = c(0.3, 0.3, 0.4))
dev.off()

#For full dataset will want to rerun PCA, determine number of axes of pop structure, and then recalculate with that number K. Then test statistics for outlier SNPs can be assessed

x_K5 <- pcadapt(input = makeBed.file, K = 5, ploidy = 2) #K must be less than the number of species

plot(x_K5, option = "manhattan") #This is useful for full dataset
plot(x_K5, option = "qqplot") #Another representation of p-value for SNP
plot(x_K5, option = "stat.distribution") #chi.sq distribution
