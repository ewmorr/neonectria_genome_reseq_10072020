require(dplyr)
require(ggplot2)
require(tidyr)
require(gridExtra)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#reading in different levels of K for admixture plotting. Using PC 2-8 based on pcadapt analyses (obviously K = 1 is pointless)

K4 = read.table("Nf_SPANDx_all_seqs_str_th/str_K4_rep5.r_format.txt", header = F)

#each row represents a sample, and each column represents the proportion of the ancestral genotype from which that sample was derived

############
#Read meatdata
fam_info = read.table("Nf_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
sample_ID_map = read.table("sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[3:4] = c("DNA_label", "sample")
site_info = read.table("sample_metadata/site_coords.txt", header = T)

sample_metadata = left_join(
data.frame(sample = fam_info[,1]),
    sample_ID_map %>% select(sample,Site)
) %>% left_join(.,
    site_info %>% select(Site, state.name)
) %>% select(sample, state.name)

##########



#state_order = c("WV", "NC", "ME.N", "ME.S", "MI", "NJ", "NY.N", "NY.S", "PA", "NH") #This is based on examining the K4 plot. May want to change order of Me.N and Me.S towards the right

state_order = c("MI", "NC", "WV", "VA", "NJ", "PA.W", "PA", "NY.W", "NY.S", "NY.N", "CT", "NH.CW", "NH.CCM", "NH.BART", "NH.SCG", "ME.S", "ME.N") #REset to approximately coincide with site distance


######
#K4 plot
colnames(K4) = c("sample", "V1", "V2", "V3", "V4")

K4.meta = left_join(K4, sample_metadata)

sample_order.K4 = K4.meta[with(K4.meta, order(state.name, V1, V2, V3, V4)),] %>% select(sample)

K4.long = pivot_longer(K4.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K4.long$sample = factor(K4.long$sample, levels = as.factor(sample_order.K4$sample))

p.K4 = ggplot(K4.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(state.name, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 4", x = "Genotype") +
scale_fill_manual(values = rev(cbPalette)) +
my_gg_theme +
theme(
    strip.background=element_rect(colour="black", fill=NA),
    strip.text.x = element_text(angle = 90, size = 30),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 30),
    axis.text.y = element_text(size = 25),

)

pdf("figures/Nf.LD_filtered.structure.K4.pdf", width = 48, height = 4)
p.K4
dev.off()

