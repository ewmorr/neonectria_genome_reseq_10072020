require(tidyr)

fam_info = read.table("Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)

ind_file = left_join(data.frame(sample = fam_info[,1]), sample_metadata %>% select(sample,State) )

write.table(ind_file, "Nf_post_SPANDx/LD_filter/ind_file.structure", quote = F, sep = "\t", row.names = F, col.names = F)
