require(dplyr)
require(ggplot2)
source("R_scripts/ggplot_theme.txt")

#granular terms
granular = read.table("data/blast2GO/makerFINAL.all.maker.proteins-1.long.uniq.txt", header = T, sep = "\t", quote = "\"")
#yeast GOslim
yeast = read.table("data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.long.uniq.txt", header = T, sep = "\t", quote = "\"")
#generic GOslim
generic = read.table("data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.long.uniq.txt", header = T, sep = "\t", quote = "\"")

########
#Basic counts
#number genes
granular$SeqName %>% unique %>% length
yeast$SeqName %>% unique %>% length
generic$SeqName %>% unique %>% length
#14064

#count number of aspect hits and NA
#aspect counts are the total number of hits to the aspect, and one gene can hit more than once. NA is just the number of NA genes.
granular %>%
    group_by(GOaspect) %>%
    summarize(n = n())
#NA = 4433
14064-4433
#9631, corresponds to the number of high-quality B2GO mappings (as it should)
yeast %>%
    group_by(GOaspect) %>%
    summarize(n = n())
#NA = 4439
#This difference of six NA from granular to yeast corresponds to the number of "obsolete" mappings (i.e., marked as ANNOTATED in the slim)
generic %>%
    group_by(GOaspect) %>%
    summarize(n = n())
#NA = 5902
#difference corresponds to 1469 with ANNOTATED instead of GO-SLIM flag

#######
#Number of genes mapped within different aspects. Genes can hit more than one aspect

#granular
granular %>% filter(GOaspect == "C") %>% pull(SeqName) %>% unique() %>% length()
#5875
granular %>% filter(GOaspect == "F") %>% pull(SeqName) %>% unique() %>% length()
#7841
granular %>% filter(GOaspect == "P") %>% pull(SeqName) %>% unique() %>% length()
#6045

#yeast
yeast %>% filter(GOaspect == "C") %>% pull(SeqName) %>% unique() %>% length()
#3919
yeast %>% filter(GOaspect == "F") %>% pull(SeqName) %>% unique() %>% length()
#7573
yeast %>% filter(GOaspect == "P") %>% pull(SeqName) %>% unique() %>% length()
#5041

#generic
generic %>% filter(GOaspect == "C") %>% pull(SeqName) %>% unique() %>% length()
#2902
generic %>% filter(GOaspect == "F") %>% pull(SeqName) %>% unique() %>% length()
#7120
generic %>% filter(GOaspect == "P") %>% pull(SeqName) %>% unique() %>% length()
#4927

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
#generic
generic.C = generic %>% filter(GOaspect == "C") %>% group_by(GOname) %>% summarize(n = n())
generic.F = generic %>% filter(GOaspect == "F") %>% group_by(GOname) %>% summarize(n = n())
generic.P = generic %>% filter(GOaspect == "P") %>% group_by(GOname) %>% summarize(n = n())

#number cats
#granular
granular.C %>% nrow()
#643
granular.F %>% nrow()
#1571
granular.P %>% nrow()
#1971

#yeast
yeast.C %>% nrow()
#28
yeast.F %>% nrow()
#38
yeast.P %>% nrow()
#71

#generic
generic.C %>% nrow()
#25
generic.F %>% nrow()
#34
generic.P %>% nrow()
#60

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

#generic C,F,P
p4 = ggplot(generic.C,
            aes(x = reorder(GOname, -n), y = n)
) +
    geom_col() +
    labs(y = "No. mappings", x = "C: Cellular component") +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, 
                                   hjust = 1,
                                   size = 14
        )
    )
p4

p5 = ggplot(generic.F,
            aes(x = reorder(GOname, -n), y = n)
) +
    geom_col() +
    labs(y = "No. mappings", x = "F: Molecular function") +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, 
                                   hjust = 1,
                                   size = 14
        )
    )
p5

p6 = ggplot(generic.P,
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
p6

pdf("figures/go_slim_overview.yeast.pdf", width = 24, height = 8)
p1
p2
p3
dev.off()

pdf("figures/go_slim_overview.generic.pdf", width = 24, height = 8)
p4
p5
p6
dev.off()
