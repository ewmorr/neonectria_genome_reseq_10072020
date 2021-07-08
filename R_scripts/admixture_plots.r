require(dplyr)
require(ggplot2)
require(tidyr)
require(gridExtra)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#reading in different levels of K for admixture plotting. Using PC 2-8 based on pcadapt analyses (obviously K = 1 is pointless)

K2 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.2.Q", header = F)
K3 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.3.Q", header = F)
K4 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.4.Q", header = F)
K5 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.5.Q", header = F)
K6 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.6.Q", header = F)
K7 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.7.Q", header = F)
K8 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.8.Q", header = F)
K9 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.9.Q", header = F)
K10 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.10.Q", header = F)
K11 = read.table("Nf_post_SPANDx/LD_filter/admixture/Nf.out.filtered.LD_filtered_0.5_10Kb.11.Q", header = F)

#each row represents a sample, and each column represents the proportion of the ancestral genotype from which that sample was derived
fam_info = read.table("Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)
state_order = c("WV", "NC", "ME.N", "ME.S", "MI", "NJ", "NY.N", "NY.S", "PA", "NH") #This is based on examining the K2 plot. May want to change order of Me.N and Me.S towards the right

######
#K2 plot
K2.fam = data.frame(K2, sample = fam_info[,1])
K2.meta = left_join(K2.fam, sample_metadata)

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
    strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
)

######
#K3 plot
K3.fam = data.frame(K3, sample = fam_info[,1])
K3.meta = left_join(K3.fam, sample_metadata)

sample_order.K3 = K3.meta[with(K3.meta, order(State, V1, V2, V3)),] %>% select(sample)

K3.long = pivot_longer(K3.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K3.long$sample = factor(K3.long$sample, levels = as.factor(sample_order.K3$sample))

p.K3 = ggplot(K3.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 3", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
strip.text.x = element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank()
)

######
#K4 plot
K4.fam = data.frame(K4, sample = fam_info[,1])
K4.meta = left_join(K4.fam, sample_metadata)

sample_order.K4 = K4.meta[with(K4.meta, order(State, V1, V2, V3)),] %>% select(sample)

K4.long = pivot_longer(K4.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K4.long$sample = factor(K4.long$sample, levels = as.factor(sample_order.K4$sample))

p.K4 = ggplot(K4.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 4", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
strip.text.x = element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank()
)

######
#K5 plot
K5.fam = data.frame(K5, sample = fam_info[,1])
K5.meta = left_join(K5.fam, sample_metadata)

sample_order.K5 = K5.meta[with(K5.meta, order(State, V1, V2, V3)),] %>% select(sample)

K5.long = pivot_longer(K5.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K5.long$sample = factor(K5.long$sample, levels = as.factor(sample_order.K5$sample))

p.K5 = ggplot(K5.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 5", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
strip.text.x = element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank()
)

######
#K6 plot
K6.fam = data.frame(K6, sample = fam_info[,1])
K6.meta = left_join(K6.fam, sample_metadata)

sample_order.K6 = K6.meta[with(K6.meta, order(State, V1, V2, V3)),] %>% select(sample)

K6.long = pivot_longer(K6.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K6.long$sample = factor(K6.long$sample, levels = as.factor(sample_order.K6$sample))

p.K6 = ggplot(K6.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 6", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
strip.text.x = element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank()
)

######
#K7 plot
K7.fam = data.frame(K7, sample = fam_info[,1])
K7.meta = left_join(K7.fam, sample_metadata)

sample_order.K7 = K7.meta[with(K7.meta, order(State, V1, V2, V3)),] %>% select(sample)

K7.long = pivot_longer(K7.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K7.long$sample = factor(K7.long$sample, levels = as.factor(sample_order.K7$sample))

p.K7 = ggplot(K7.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 7", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
strip.text.x = element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank()
)

######
#K8 plot
K8.fam = data.frame(K8, sample = fam_info[,1])
K8.meta = left_join(K8.fam, sample_metadata)

sample_order.K8 = K8.meta[with(K8.meta, order(State, V1, V2, V3)),] %>% select(sample)

K8.long = pivot_longer(K8.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K8.long$sample = factor(K8.long$sample, levels = as.factor(sample_order.K8$sample))

p.K8 = ggplot(K8.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 8", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
axis.text.x = element_blank()
)

######
#K9 plot

state_order = c("WV", "NC", "ME.N", "ME.S", "MI", "NJ", "NY.N", "NY.S", "PA", "NH") #REset to approximately coincide with site distance
state_order = c("MI", "NC", "WV", "PA", "NJ", "NY.S", "NY.N", "NH", "ME.S", "ME.N") #REset to approximately coincide with site distance

K9.fam = data.frame(K9, sample = fam_info[,1])
K9.meta = left_join(K9.fam, sample_metadata)

sample_order.K9 = K9.meta[with(K9.meta, order(State, V1, V2, V3, V5, V6, V7, V8, V9)),] %>% select(sample)

K9.long = pivot_longer(K9.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K9.long$sample = factor(K9.long$sample, levels = as.factor(sample_order.K9$sample))

p.K9 = ggplot(K9.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 9", x = "Genotype") +
scale_fill_manual(values = twelvePaired) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
strip.text.x = element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank()
)

######
#K10 plot
K10.fam = data.frame(K10, sample = fam_info[,1])
K10.meta = left_join(K10.fam, sample_metadata)

sample_order.K10 = K10.meta[with(K10.meta, order(State, V1, V2, V3, V5, V6, V7, V8, V9, V10)),] %>% select(sample)

K10.long = pivot_longer(K10.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K10.long$sample = factor(K10.long$sample, levels = as.factor(sample_order.K10$sample))

p.K10 = ggplot(K10.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 10", x = "Genotype") +
scale_fill_manual(values = twelvePaired) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
strip.text.x = element_blank(),
axis.title.x = element_blank(),
axis.text.x = element_blank()
)

######
#K11 plot
K11.fam = data.frame(K11, sample = fam_info[,1])
K11.meta = left_join(K11.fam, sample_metadata)

sample_order.K11 = K11.meta[with(K11.meta, order(State, V1, V2, V3, V5, V6, V7, V8, V9, V10, V11)),] %>% select(sample)

K11.long = pivot_longer(K11.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K11.long$sample = factor(K11.long$sample, levels = as.factor(sample_order.K11$sample))

p.K11 = ggplot(K11.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(State, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 11", x = "Genotype") +
scale_fill_manual(values = twelvePaired) +
my_gg_theme +
theme(
strip.background=element_rect(colour="black", fill=NA),
axis.text.x = element_blank()
)

pdf("figures/Nf.LD_filtered.admixture.pdf", width = 20, height = 16)
grid.arrange(p.K2,p.K3,p.K4,p.K5,p.K6,p.K7,p.K8,ncol = 1, heights = c(rep(0.14, 6), 0.16) )
dev.off()

pdf("figures/Nf.LD_filtered.admixture.K9-K11.pdf", width = 20, height = 8)
grid.arrange(p.K9,p.K10,p.K11,ncol = 1, heights = c(0.3,0.3,0.4) )
dev.off()

