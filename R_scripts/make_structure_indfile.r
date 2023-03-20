require(tidyr)
require(dplyr)

fam_info = read.table("data/Nd_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[3:4] = c("DNA_label", "sample")
site_info = read.table("data/sample_metadata/site_coords.txt", header = T)

ind_file = left_join(
    data.frame(sample = fam_info[,1]),
        sample_ID_map %>% select(sample,Site)
    ) %>% left_join(.,
        site_info %>% select(Site, state.name)
    ) %>% select(sample, state.name)

write.table(ind_file, "data/Nd_SPANDx_all_seqs/ind_file.structure", quote = F, sep = "\t", row.names = F, col.names = F)
