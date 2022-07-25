source("~/repo/neonectria_genome_reseq_10072020/R_scripts/make_site_metadata.r")

print(sample_metadata %>% group_by(state.name) %>% summarize(n()))
