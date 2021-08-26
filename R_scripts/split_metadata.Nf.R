require(dplyr)

metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)

#Will need to adjust this if site/geographical coding becomes more complex, e.g., multiple independent sites in a single state. For now we are coding independent sites within a state with unique modifiers (e.g., NY.N vs NY.S)

levels_site = levels(metadata$State)

pairwise_comps = combn(levels_site, 2)

setwd("sample_metadata")
system("mkdir pairwise")
setwd("pairwise")

write.table(t(pairwise_comps), file = "pairwise_site_combs.txt", quote = F, sep =  "\t", row.names = F, col.names = F)


for(i in 1:length(levels_site)){
    temp.tab = metadata %>% filter(State == levels_site[i]) %>% select(sample)
    write.table(temp.tab, file = paste(levels_site[i], ".txt", sep = ""), col.names = F, row.names = F, quote = F)
}


