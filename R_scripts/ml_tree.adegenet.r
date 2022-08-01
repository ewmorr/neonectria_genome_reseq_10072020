require(vcfR)
require(tidyr)
require(dplyr)
require(ggplot2)
require(pegas)
require(adegenet)
require(ape)
#require(dartR) #no longer needed bc not using gl2gi
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#metadata
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/make_site_metadata.r")
site_coords = read.table("sample_metadata/site_coords_for_map.txt", header = T)
ind.metrics = left_join(data.frame(sample = fam_info[,1]), sample_metadata %>% select(sample,state.name) ) %>%
left_join(., site_coords %>% select(state.name, lat, lon))
#########


#filtered VCF
vcf <- read.vcfR("Nf_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.vcf", verbose = FALSE)
gl = vcfR2genlight(vcf)

#Set metadata in genLight
gl@other$ind.metrics = ind.metrics
gl@other$latlong = ind.metrics[,3:4]
gl@other$latlong
gl@pop = as.factor(ind.metrics$state.name) #need to set pop for the ibd test to work
gl@ploidy = rep(as.integer(1), nInd(gl))
nPop(gl)
pop(gl)

# https://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf p. 42
