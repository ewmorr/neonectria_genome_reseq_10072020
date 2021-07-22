require(vcfR)
require(tidyr)
require(dplyr)
require(ggplot2)
require(pegas)
require(dartR)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#filtered VCF
vcf <- read.vcfR("Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf", verbose = FALSE)
gl = vcfR2genlight(vcf)
#Warning message:
#In vcfR2genlight(vcf) : Found 3027 loci with more than two alleles.
#Objects of class genlight only support loci with two alleles.
#3027 loci will be omitted from the genlight object.

#This might be why PDSpider is some dropping alleles
#26354-3027 = 23327
#so it's not all of them (~21K left after PDSpider)

gi = gl2gi(gl)

fam_info = read.table("Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.fam", header = F)
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)
site_coords = read.table("sample_metadata/site_coords.txt", header = T)

ind.metrics = left_join(data.frame(sample = fam_info[,1]), sample_metadata %>% select(sample,State,Site) ) %>%
left_join(., site_coords %>% select(Site, lat, lon)) %>% select(-Site)

gl@other$ind.metrics = ind.metrics
gl@other$latlong = ind.metrics[,3:4]
gl@other$latlong
gl@pop = ind.metrics$State #need to set pop for the ibd test to work
gl@ploidy = rep(as.integer(1), 65)
nPop(gl)
pop(gl)

#remove the ME.N and NJ samples because of small sample size
rm.ind.list = data.frame(gl@ind.names, pop(gl)) %>% filter(pop.gl. == "ME.N" | pop.gl. == "NJ")
gl.subset = gl[!gl@ind.names %in% rm.ind.list$gl.ind.names]

x = gl.subset
#calculate Mercator prjected distance
xy <- dismo::Mercator(x@other$latlong[,c("lon","lat")])

x.Fst = StAMPP::stamppFst(x, nboots=1, nclusters = 4) #nclusters is the number of threads. Set to four for local
Dgen <- as.dist(x.Fst)


#Most values are negative and should be treated as zero (negative Fst doesn't make much sense)
#positive values are between NH-NY.S, NH-WV, NH-NC

x.Fst.pval = StAMPP::stamppFst(x, nboots=100, nclusters = 4, percent = 95) #nclusters is the number of threads. Set to four for local

#Fst test produces the same as from dartR
