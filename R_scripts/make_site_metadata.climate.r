library(dplyr)

fam_info.Nf = read.table("data/Nf_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
fam_info.Nf$sp = "Nf"
fam_info.Nd = read.table("data/Nd_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
fam_info.Nd$sp = "Nd"
fam_info = rbind(fam_info.Nf, fam_info.Nd)

sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
colnames(sample_ID_map)[3:4] = c("DNA_label", "sample")
site_info = read.table("data/sample_metadata/site_coords.txt", header = T)


sample_metadata = left_join(
data.frame(sample = fam_info[,1], sp = fam_info$sp),
    sample_ID_map %>% select(sample,Site)
) %>% left_join(.,
    site_info %>% select(Site, state.name)
) %>% select(sample, state.name, sp)

site.info = read.csv("data/sample_metadata/site_info.csv")
site.GDD = read.table("data/sample_metadata/site_climate.GDD.txt", header = T)
site.climate = read.table("data/sample_metadata/sites_climate.txt", header = T)
site.GDD$freezeThaw.annual_mean = site.GDD$freezeThaw.mean_growing + site.GDD$freezeThaw.mean_nongrowing

site_metadata = left_join(site.GDD, site.info %>% select(Site, lat, lon, duration_infection), by = "Site") %>%
    left_join(., site.climate %>% select(Site, tmin, tmax, ppt, MAT, lat, lon, state.name), by = c("Site", "lat", "lon") )

sample_metadata.site_info = left_join(sample_metadata, site_metadata, by = "state.name")
nrow(sample_metadata.site_info)

write.csv(sample_metadata.site_info, file = "data/sample_metadata/sample_metadata.filtered_samples.csv", quote = F, row.names = F)
