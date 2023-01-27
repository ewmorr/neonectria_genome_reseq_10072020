require(dplyr)
require(ggplot)
source("R_scripts/ggplot_theme.txt")

yeast = read.table("data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.long.uniq.txt", header = T, sep = "\t", quote = "\"")

generic = read.table("data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.long.uniq.txt", header = T, sep = "\t", quote = "\"")

