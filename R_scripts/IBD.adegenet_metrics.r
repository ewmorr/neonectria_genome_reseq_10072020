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
rm.ind.list = data.frame(gl@ind.names, pop(gl)) %>% filter(pop.gl. == "ME.N" | pop.gl. == "NJ" | pop.gl. == "MI")
gl.subset = gl[!gl@ind.names %in% rm.ind.list$gl.ind.names]

x = gl.subset
#calculate Mercator prjected distance
xy <- dismo::Mercator(x@other$latlong[,c("lon","lat")])

#Need to convert to genpop for adegenet dist.genpop
gi.subset = gl2gi(x)
gi.subset@ploidy = rep(as.integer(1), 59)

gp = genind2genpop(gi.subset)

Dgen.1 <- dist.genpop(gp,method=1)
Dgen.2 <- dist.genpop(gp,method=2)
Dgen.3 <- dist.genpop(gp,method=3)
Dgen.4 <- dist.genpop(gp,method=4)
Dgen.5 <- dist.genpop(gp,method=5)

#Get pariwise site distances
site_coords.subset = site_coords %>% filter(state.name %in% c("ME.S", "NC", "NH.CW", "NY.N", "NY.S", "PA", "WV"))
site_order = c("ME.S", "NC", "NH.CW", "NY.N", "NY.S", "PA", "WV")
site_coords.subset.order = site_coords.subset[match(site_order, site_coords.subset$state.name), ]
rownames(site_coords.subset.order) = site_coords.subset.order$state.name

Dgeo <- dist(dismo::Mercator(site_coords.subset.order[,c("lon", "lat")]))

ibd.1 <- mantel.randtest(Dgen.1,Dgeo)
ibd.2 <- mantel.randtest(Dgen.2,Dgeo)
ibd.3 <- mantel.randtest(Dgen.3,Dgeo)
ibd.4 <- mantel.randtest(Dgen.4,Dgeo)
ibd.5 <- mantel.randtest(Dgen.5,Dgeo)

ibd.1.log <- mantel.randtest(Dgen.1,log(Dgeo))
ibd.2.log <- mantel.randtest(Dgen.2,log(Dgeo))
ibd.3.log <- mantel.randtest(Dgen.3,log(Dgeo))
ibd.4.log <- mantel.randtest(Dgen.4,log(Dgeo))
ibd.5.log <- mantel.randtest(Dgen.5,log(Dgeo))

require(reshape2)
Dgeo.long = subset(melt(Dgeo %>% as.matrix), value != 0)
Dgen.long = subset(melt(Dgen.2 %>% as.matrix), value != 0)

Dgen.Dgeo = data.frame(gen = Dgen.long$value, geo = Dgeo.long$value)

p1 = ggplot(Dgen.Dgeo, aes(x = geo/1000, y = gen)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "Pairwise site geographic distance (km)", y = expression(paste("Pairwise site genetic distance (D"["CSE"], ")"))) +
my_gg_theme

p2 = ggplot(Dgen.Dgeo, aes(x = geo/1000, y = gen)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "Pairwise site geographic distance (km)", y = expression(paste("Pairwise site genetic distance (D"["CSE"], ")"))) +
scale_x_continuous(trans = "log", breaks = c(200, 500, 1000, 2000)) +
my_gg_theme

pdf("figures/IBD.Dcse_test.pdf", width = 9, height = 5)
p1
dev.off()

