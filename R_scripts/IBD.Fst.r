require(vcfR)
require(tidyr)
require(dplyr)
require(ggplot2)
require(pegas)
require(dartR)
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
gl@ploidy = rep(as.integer(1), nInd(gl)
nPop(gl)
pop(gl)

#SOME SITES HAVE SMALL SAMPLE SIZE
#Set min sample size to either 3 or 4
low_n_states = ( (sample_metadata %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 4) )$state.name

rm.ind.list = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(pop.gl %in% low_n_states)
gl.subset = gl[!gl@ind.names %in% rm.ind.list$gl.ind.names]

gl.subset
#calculate Mercator prjected distance
xy <- dismo::Mercator(gl.subset@other$latlong[,c("lon","lat")])

#calculate Fst
gl.subset.Fst = StAMPP::stamppFst(gl.subset, nboots=1, nclusters = 4) #nclusters is the number of threads. Set to four for local
Dgen <- as.dist(gl.subset.Fst)


#Most values are negative and should be treated as zero (negative Fst doesn't make much sense)
#positive values are between NH-NY.S, NH-WV, NH-NC

x.Fst.pval = StAMPP::stamppFst(x, nboots=100, nclusters = 4, percent = 95) #nclusters is the number of threads. Set to four for local

#Fst test produces the same as from dartR
