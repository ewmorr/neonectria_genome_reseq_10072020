source("~/repo/neonectria_genome_reseq_10072020/R_scripts/make_site_metadata.r")
low_n_states = ( (sample_metadata %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 4) )$state.name

retain_samps = (sample_metadata %>% filter(!state.name %in% low_n_states) )$sample

write.table(retain_samps, file = "Nf_SPANDx_all_seqs/retain_samples.txt", quote = F, col.names = F, row.names = F)
