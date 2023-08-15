library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
source("R_scripts/ggplot_theme.txt")

K3 = read.table("data/Nf_SPANDx_all_seqs_str_th/str_K3_rep1.r_format.txt", header = F)

#each row represents a sample, and each column represents the proportion of the ancestral genotype from which that sample was derived

############
#Read meatdata
fam_info = read.table("data/Nf_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[3:4] = c("DNA_label", "sample")
site_info = read.table("data/sample_metadata/site_coords.txt", header = T)

sample_metadata = left_join(
        data.frame(sample = fam_info[,1]),
        sample_ID_map %>% select(sample,Site)
    ) %>% left_join(.,
        site_info %>% select(Site, state.name)
    ) %>% 
    select(sample, state.name)

##########



#state_order = c("WV", "NC", "ME.N", "ME.S", "MI", "NJ", "NY.N", "NY.S", "PA", "NH") #This is based on examining the K3 plot. May want to change order of Me.N and Me.S towards the right

state_order = c("MI", "NC", "WV", "VA", "NJ", "PA.W", "PA", "NY.W", "NY.S", "NY.N", "CT", "NH.CW", "NH.CCM", "NH.BART", "NH.SCG", "ME.S", "ME.N") #REset to approximately coincide with site distance


######
#K3 plot
colnames(K3) = c("sample", "V1", "V2", "V3")

K3.meta = left_join(K3, sample_metadata)

sample_order.K3 = K3.meta[with(K3.meta, order(state.name, V1, V2, V3)),] %>% select(sample)

K3.long = pivot_longer(K3.meta, names_to = "ancestor", values_to = "proportion", cols = starts_with("V"))
K3.long$sample = factor(K3.long$sample, levels = as.factor(sample_order.K3$sample))

p.K3 = ggplot(K3.long, aes(sample, proportion, fill = ancestor)) +
geom_bar(stat = "identity") +
facet_grid(~factor(state.name, levels = state_order), scales = "free_x", switch = "x", space = "free_x") + #Note the use of levels here at the facet argument. This is a pain in the butt to get to work (ordering a factor for index level ordering and then also ordering facets
guides(fill = "none") +
labs(y = "K = 3", x = "Genotype") +
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
p.K3

pdf("figures/Nf.LD_filtered.structure.K3.pdf", width = 48, height = 4)
p.K3
dev.off()

