library(vcfR)
library(tidyr)
library(ggplot2)
source("R_scripts/ggplot_theme.txt")

#filtered VCF
vcf <- read.vcfR("data/Nf_SPANDx_all_seqs/out.vcf", verbose = FALSE)

#just pulling in dp to get sample IDs

sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()
is.matrix(dp)

#sample_ids
first_set.ids = paste0("NG", 1:99) 
second_set.ids = paste0("NG", 101:163)
third_set.ids = paste0("NG", 170:196)

lib_one = first_set.ids[first_set.ids %in% sample_ids]
lib_two = second_set.ids[second_set.ids %in% sample_ids]
lib_three = third_set.ids[third_set.ids %in% sample_ids]

write.table(lib_one, "data/Nf/lib_one_ids.txt")
write.table(lib_two, "data/Nf/lib_two_ids.txt")
write.table(lib_three, "data/Nf/lib_three_ids.txt")