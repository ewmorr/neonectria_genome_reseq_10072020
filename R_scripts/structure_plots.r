require(dplyr)
require(ggplot2)
require(tidyr)
require(gridExtra)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#reading in different levels of K for admixture plotting. Using PC 2-8 based on pcadapt analyses (obviously K = 1 is pointless)

K2 = read.table("Nf_post_SPANDx/LD_filter/structure_threader/str_K2_rep2_f.r_format.txt", header = F)

#each row represents a sample, and each column represents the proportion of the ancestral genotype from which that sample was derived
fam_info = read.table("Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)
state_order = c("WV", "NC", "ME.N", "ME.S", "MI", "NJ", "NY.N", "NY.S", "PA", "NH") #This is based on examining the K2 plot. May want to change order of Me.N and Me.S towards the right
state_order = c("MI", "NC", "WV", "PA", "NJ", "NY.S", "NY.N", "NH", "ME.S", "ME.N") #REset to approximately coincide with site distance


######
#K2 plot
colnames(K2) = c("sample", "V1", "V2")

K2.meta = left_join(K2, sample_metadata)

sample_order.K2 = K2.meta[with(K2.meta, order(State, V1, V2)),] %>% select(sample)

K2.long = pivot_longer(K2.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K2.long$sample = factor(K2.long$sample, levels = as.factor(sample_order.K2$sample))

p.K2 = ggplot(K2.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 2", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
    strip.background=element_rect(colour="black", fill=NA),
#    strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
)

pdf("figures/Nf.LD_filtered.structure.K2.pdf", width = 22, height = 4)
p.K2
dev.off()

