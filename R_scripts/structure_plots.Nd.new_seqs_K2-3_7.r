require(dplyr)
require(ggplot2)
require(tidyr)
require(gridExtra)
source("R_scripts/ggplot_theme.txt")

#reading in different levels of K for admixture plotting. Using PC 2-8 based on pcadapt analyses (obviously K = 1 is pointless)

K2 = read.table("data/Nd_SPANDx_all_seqs_str_th/str_K2_rep5_f.r_format.txt", header = F)
K3 = read.table("data/Nd_SPANDx_all_seqs_str_th/str_K3_rep5_f.r_format.txt", header = F)
K4 = read.table("data/Nd_SPANDx_all_seqs_str_th/str_K4_rep5_f.r_format.txt", header = F)
K7 = read.table("data/Nd_SPANDx_all_seqs_str_th/str_K7_rep5_f.r_format.txt", header = F)

#each row represents a sample, and each column represents the proportion of the ancestral genotype from which that sample was derived

############
#Read meatdata
fam_info = read.table("data/Nd_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[3:4] = c("DNA_label", "sample")
site_info = read.table("data/sample_metadata/site_coords.txt", header = T)

sample_metadata = left_join(
data.frame(sample = fam_info[,1]),
    sample_ID_map %>% select(sample,Site)
) %>% left_join(.,
    site_info %>% select(Site, state.name)
) %>% select(sample, state.name)

##########



#state_order = c("WV", "NC", "ME.N", "ME.S", "MI", "NJ", "NY.N", "NY.S", "PA", "NH") #This is based on examining the K4 plot. May want to change order of Me.N and Me.S towards the right

state_order = c("ME.N", "NH.CCM", "NH.CW", "PA", "WV", "MI") #REset to approximately coincide with site distance


######
#plot
colnames(K2) = c("sample", "V1", "V2")
colnames(K3) = c("sample", "V1", "V2", "V3")
colnames(K4) = c("sample", "V1", "V2", "V3", "V4")
colnames(K7) = c("sample", "V1", "V2", "V3", "V4", "V5", "V6", "V7")

#join meta
K2.meta = left_join(K2, sample_metadata)
K3.meta = left_join(K3, sample_metadata)
K4.meta = left_join(K4, sample_metadata)
K7.meta = left_join(K7, sample_metadata)

#order df by ancestory sorting
sample_order.K2 = K2.meta[with(K2.meta, order(state.name, V1, V2)),] %>% select(sample)
sample_order.K3 = K3.meta[with(K3.meta, order(state.name, V1, V2, V3)),] %>% select(sample)
sample_order.K4 = K4.meta[with(K4.meta, order(state.name, V1, V2, V3, V4)),] %>% select(sample)
sample_order.K7 = K7.meta[with(K7.meta, order(state.name, V1, V2, V3, V4, V5, V6, V7)),] %>% select(sample)

#long format
K2.long = pivot_longer(K2.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K3.long = pivot_longer(K3.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K4.long = pivot_longer(K4.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K7.long = pivot_longer(K7.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))

#order factor level
K2.long$sample = factor(K2.long$sample, levels = as.factor(sample_order.K2$sample))
K3.long$sample = factor(K3.long$sample, levels = as.factor(sample_order.K3$sample))
K4.long$sample = factor(K4.long$sample, levels = as.factor(sample_order.K4$sample))
K7.long$sample = factor(K7.long$sample, levels = as.factor(sample_order.K7$sample))

#PLOTS
p.K2 = ggplot(K2.long, aes(sample, proportion, fill = ancestor)) +
    geom_bar(stat = "identity") +
    facet_grid(~factor(state.name, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
    guides(fill = "none") +
    labs(y = "K = 2", x = "Genotype") +
    scale_fill_manual(values = rev(cbPalette)) +
    my_gg_theme +
    theme(
        strip.background=element_rect(colour="black", fill=NA),
        strip.text.x = element_text(angle = 0, size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        
    )
p.K2

p.K3 = ggplot(K3.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(state.name, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 3", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
    strip.background=element_rect(colour="black", fill=NA),
    strip.text.x = element_text(angle = 0, size = 15),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 15),

)
p.K3

p.K4 = ggplot(K4.long, aes(sample, proportion, fill = ancestor)) +
    geom_bar(stat = "identity") +
    facet_grid(~factor(state.name, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
    guides(fill = "none") +
    labs(y = "K = 4", x = "Genotype") +
    scale_fill_manual(values = rev(cbPalette)) +
    my_gg_theme +
    theme(
        strip.background=element_rect(colour="black", fill=NA),
        strip.text.x = element_text(angle = 0, size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        
    )
p.K4

p.K7 = ggplot(K7.long, aes(sample, proportion, fill = ancestor)) +
    geom_bar(stat = "identity") +
    facet_grid(~factor(state.name, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
    guides(fill = "none") +
    labs(y = "K = 7", x = "Genotype") +
    scale_fill_manual(values = rev(cbPalette)) +
    my_gg_theme +
    theme(
        strip.background=element_rect(colour="black", fill=NA),
        strip.text.x = element_text(angle = 0, size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        
    )
p.K7


pdf("figures/Nd.LD_filtered.structure.K2-7.pdf", width = 12, height = 8)
grid.arrange(p.K2,p.K3,p.K4,p.K7, nrow = 4)
dev.off()

