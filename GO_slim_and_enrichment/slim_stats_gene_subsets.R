library(dplyr)
library(ggplot2)
source("R_scripts/ggplot_theme.txt")

#granular terms
granular = read.table("data/blast2GO/Nf_granular_GO.long.uniq.txt", header = T, sep = "\t", quote = "\"")
#yeast GOslim
yeast = read.table("data/blast2GO/Nf_GOslim_Yeast.long.uniq.txt", header = T, sep = "\t", quote = "\"")

#Gene subsets
freezeThaw_unique = read.table("data/Nf_LFMM_tables/freezeThaw.geneIDs.unique.txt", header = F)
hdd4_ppt_shared = read.table("data/Nf_LFMM_tables/hdd4_ppt_shared.geneIDs.shared.txt", header = F)
ppt_unique = read.table("data/Nf_LFMM_tables/ppt.geneIDs.unique.txt", header = F)
hdd4_unique = read.table("data/Nf_LFMM_tables/hdd4.geneIDs.unique.txt", header = F)

colnames(freezeThaw_unique) = "SeqName"
colnames(hdd4_ppt_shared) = "SeqName"
colnames(ppt_unique) = "SeqName"
colnames(hdd4_unique) = "SeqName"

####################################
#Assessing P aspect as this is most common approach
#For us F aspect may be better because there are more annotations

yeast.P = yeast %>% filter(GOaspect == "P")
yeast.F = yeast %>% filter(GOaspect == "F")

granular.P = granular %>% filter(GOaspect == "P")
granular.F = granular %>% filter(GOaspect == "F")

#yeast

#Freeze-thaw
freezeThaw_unique.GO.P = yeast.P %>%
    filter(SeqName %in% freezeThaw_unique$SeqName)
freezeThaw_unique.GO.F = yeast.F %>%
    filter(SeqName %in% freezeThaw_unique$SeqName)

#hdd4
hdd4_unique.GO.P = yeast.P %>%
    filter(SeqName %in% hdd4_unique$SeqName)
hdd4_unique.GO.F = yeast.F %>%
    filter(SeqName %in% hdd4_unique$SeqName)

#ppt
ppt_unique.GO.P = yeast.P %>%
    filter(SeqName %in% ppt_unique$SeqName)
ppt_unique.GO.F = yeast.F %>%
    filter(SeqName %in% ppt_unique$SeqName)

#granular

#Freeze-thaw
freezeThaw_unique.GO_gran.P = granular.P %>%
    filter(SeqName %in% freezeThaw_unique$SeqName)
freezeThaw_unique.GO_gran.F = granular.F %>%
    filter(SeqName %in% freezeThaw_unique$SeqName)

#hdd4
hdd4_unique.GO_gran.P = granular.P %>%
    filter(SeqName %in% hdd4_unique$SeqName)
hdd4_unique.GO_gran.F = granular.F %>%
    filter(SeqName %in% hdd4_unique$SeqName)

#ppt
ppt_unique.GO_gran.P = granular.P %>%
    filter(SeqName %in% ppt_unique$SeqName)
ppt_unique.GO_gran.F = granular.F %>%
    filter(SeqName %in% ppt_unique$SeqName)

###################
#Count summaries

freezeThaw_unique.GO.P.count = freezeThaw_unique.GO.P %>%
    group_by(GOname) %>%
    summarize(n = n())

freezeThaw_unique.GO.F.count = freezeThaw_unique.GO.F %>%
    group_by(GOname) %>%
    summarize(n = n())

hdd4_unique.GO.P.count = hdd4_unique.GO.P %>%
    group_by(GOname) %>%
    summarize(n = n())

hdd4_unique.GO.F.count = hdd4_unique.GO.F %>%
    group_by(GOname) %>%
    summarize(n = n())

ppt_unique.GO.P.count = ppt_unique.GO.P %>%
    group_by(GOname) %>%
    summarize(n = n())

ppt_unique.GO.F.count = ppt_unique.GO.F %>%
    group_by(GOname) %>%
    summarize(n = n())

###################
#Plots
#yeast C,F,P
p1 = ggplot(freezeThaw_unique.GO.P.count,
            aes(x = reorder(GOname, -n), y = n)
) +
    geom_col() +
    labs(y = "No. mappings", x = "P: Biological process") +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        axis.title.x = element_blank()
    )
p1

p2 = ggplot(hdd4_unique.GO.P.count,
    aes(x = reorder(GOname, -n), y = n)
) +
    geom_col() +
    labs(y = "No. mappings", x = "P: Biological process") +
    scale_y_continuous(breaks = c(0,1)) +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        axis.title.x = element_blank()
    )
p2

p3 = ggplot(ppt_unique.GO.P.count,
            aes(x = reorder(GOname, -n), y = n)
) +
    geom_col() +
    labs(y = "No. mappings", x = "P: Biological process") +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 65, 
                                   hjust = 1,
                                   size = 14
        )
    )
p3


pdf("figures/LFMM_gene_subsets.freezeThaw.pdf", width = 13, height = 6)
p1
dev.off()

pdf("figures/LFMM_gene_subsets.hdd4.pdf", width = 8, height = 6)
p2
dev.off()
