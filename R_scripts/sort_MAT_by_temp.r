require(dplyr)

setwd("GARNAS_neonectria_genome_reseq_10072020")

#Read the data
sites_climate = read.table("sample_metadata/sites_climate.txt", header = T)
seqd_samples = read.table("sample_metadata/first_set_MAT_metadata.txt", header = T, sep = "\t")

#Filter sequenced samples by species
seqd_samples.Nf = seqd_samples %>% filter(Spp_based_on_ITS_map == "Nf")
seqd_samples.Nd = seqd_samples %>% filter(Spp_based_on_ITS_map == "Nd")

#Filter site data to include sites with sequenced samples of either species
sites.seqd_Nf = sites_climate %>% filter(Site %in% seqd_samples.Nf$Site)
sites.seqd_Nd = sites_climate %>% filter(Site %in% seqd_samples.Nd$Site)

#Calculate 10th and 90th percentiles of site MAT for Nf sites and Nd sites
sites.seqd_Nf$MAT %>% quantile(probs = c(0.1, 0.9))
sites.seqd_Nd$MAT %>% quantile(probs = c(0.1, 0.9))

Nd_sites_MAT_quantiles = sites.seqd_Nd$MAT %>% quantile(probs = c(0.1, 0.9))


#Filter Nd sites by quantiles (we use Nd because there are less sites available)
sites.seqd_Nd %>% filter(MAT < Nd_sites_MAT_quantiles[1])
sites.seqd_Nd %>% filter(MAT > Nd_sites_MAT_quantiles[2])

#filter available samples by the two sites
seqd_samples.Nf %>% filter(Site == "TSP1" | Site == "MEN1")
seqd_samples.Nd %>% filter(Site == "TSP1" | Site == "MEN1")

