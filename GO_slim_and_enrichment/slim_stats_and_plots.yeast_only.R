require(dplyr)
require(ggplot2)
source("R_scripts/ggplot_theme.txt")

#granular terms
granular = read.table("data/blast2GO/Nf_granular_GO.long.uniq.txt", header = T, sep = "\t", quote = "\"")
#yeast GOslim
yeast = read.table("data/blast2GO/Nf_GOslim_Yeast.long.uniq.txt", header = T, sep = "\t", quote = "\"")

########
#Basic counts
#number genes
granular$SeqName %>% unique %>% length
yeast$SeqName %>% unique %>% length
#14289 (number of genes)

#count number of aspect hits and NA
#aspect counts are the total number of hits to the aspect, and one gene can hit more than once. NA is just the number of NA genes.
granular %>%
    group_by(GOaspect) %>%
    summarize(n = n())
#NA = 4661
14289-4661
#9628, corresponds to the number of high-quality B2GO mappings (as it should)
yeast %>%
    group_by(GOaspect) %>%
    summarize(n = n())
#NA = 4666
#This difference of 5 NA from granular to yeast corresponds to the number of "obsolete" mappings (i.e., marked as ANNOTATED in the slim instead of GO-SLIM flag)

#######
#Number of genes mapped within different aspects. Genes can hit more than one aspect

#granular
granular %>% filter(GOaspect == "C") %>% pull(SeqName) %>% unique() %>% length()
#5865
granular %>% filter(GOaspect == "F") %>% pull(SeqName) %>% unique() %>% length()
#7841
granular %>% filter(GOaspect == "P") %>% pull(SeqName) %>% unique() %>% length()
#6077

#yeast
yeast %>% filter(GOaspect == "C") %>% pull(SeqName) %>% unique() %>% length()
#4096
yeast %>% filter(GOaspect == "F") %>% pull(SeqName) %>% unique() %>% length()
#7588
yeast %>% filter(GOaspect == "P") %>% pull(SeqName) %>% unique() %>% length()
#5078

##############
#Number of categories and counts within different aspects
#granular
granular.C = granular %>% filter(GOaspect == "C") %>% group_by(GOname) %>% summarize(n = n())
granular.F = granular %>% filter(GOaspect == "F") %>% group_by(GOname) %>% summarize(n = n())
granular.P = granular %>% filter(GOaspect == "P") %>% group_by(GOname) %>% summarize(n = n())
#yeast
yeast.C = yeast %>% filter(GOaspect == "C") %>% group_by(GOname) %>% summarize(n = n())
yeast.F = yeast %>% filter(GOaspect == "F") %>% group_by(GOname) %>% summarize(n = n())
yeast.P = yeast %>% filter(GOaspect == "P") %>% group_by(GOname) %>% summarize(n = n())

#number cats
#granular
granular.C %>% nrow()
#648
granular.F %>% nrow()
#1576
granular.P %>% nrow()
#1991

#yeast
yeast.C %>% nrow()
#33
yeast.F %>% nrow()
#57
yeast.P %>% nrow()
#122

###########
#In general F aspect has more genes mapped (about 2k more, or close to half the total genes) but lower number of categories (like half) than P

###################
#Plots
#yeast C,F,P
p1 = ggplot(yeast.C,
            aes(x = reorder(GOname, -n), y = n)
) +
    geom_col() +
    labs(y = "No. mappings", x = "C: Cellular component") +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p1

p2 = ggplot(yeast.F,
    aes(x = reorder(GOname, -n), y = n)
) +
    geom_col() +
    labs(y = "No. mappings", x = "F: Molecular function") +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p2

p3 = ggplot(yeast.P,
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


pdf("figures/go_slim_overview.yeast.pdf", width = 24, height = 8)
p1
p2
p3
dev.off()

