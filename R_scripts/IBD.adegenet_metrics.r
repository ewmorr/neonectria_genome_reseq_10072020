require(vcfR)
require(tidyr)
require(dplyr)
require(ggplot2)
require(pegas)
require(reshape2)
require(adegenet)
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

#SOME SITES HAVE SMALL SAMPLE SIZE
#Set min sample size to either 3 or 4
#low_n_states = ( (sample_metadata %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 4) )$state.name
#rm.ind.list = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(pop.gl %in% c(low_n_states))

#rveerse filtering of above, but then randomly select four samples per site
low_n_states = ( (sample_metadata %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 4) )$state.name
keep.ind.list = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(!pop.gl %in% c(low_n_states))

keep.ind.rand = list()
pops_incl = keep.ind.list$pop.gl %>% as.character %>% unique

for(i in 1:length(pops_incl)){
    temp = keep.ind.list %>% filter(pop.gl == pops_incl[i])
    keep.ind.rand[[pops_incl[i]]] = temp[sample(1:nrow(temp), 4), ] # random sample of four rows
}

keep.ind.rand.df = bind_rows(keep.ind.rand)
##########################################

gl.subset = gl[gl@ind.names %in% keep.ind.rand.df$gl.ind.names]

gl.subset
#calculate Mercator prjected distance
xy <- dismo::Mercator(gl.subset@other$latlong[,c("lon","lat")])

#NOTE THAT gl2gi FUNCTION IS A HACK FROM DARTR AND MAY BE MESSING UP DOWNSTREAM ANALYSES BECAUSE IT RECODES AS DIPLOID
#Need to convert to genpop for adegenet dist.genpop
#gi.subset = gl2gi(gl.subset)
#gi.subset@ploidy = rep(as.integer(1), nInd(gl.subset))

#INSTEAD OF DARTR METHOD TRY THE DATAFRAME CONVERSION AS RECOMMENDED BY ADEGENET AUTHORS (https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2014-May/000840.html)
#First convert to a data.frame whih will give a table of 0, 1, NA, and then add 1 to values to have correct conversion of NA (0 is default)

#x=glSim(10,10)
#y = as.data.frame(x)
#z=df2genind(y + 1, ploidy=1) #IT APPEARS THAT THIS METHOD STILL WORKS WITHOUT THE +1 BUT WE WILL RETAIN FOR SAFETY'S SAKE
#To remove noninformative SNPs

#retrieve the colnames of sites with only one allele
to_remove = names(z@loc.n.all[z@loc.n.all == 1 ])
#get the col index
rm_indx = which(colnames(y) %in% to_remove)

z.rm = z[loc=-rm_indx]

y = as.data.frame(gl.subset)
gi.subset = df2genind(y + 1, ploidy=1)
#reset pop
gi.subset@pop = gl.subset@pop


gp = genind2genpop(gi.subset)


#There are several distance metrics available
#Method 2 is "Angular distance or Edward's distance" D[CSE]
#Dcse performs well in IBD tests "Sere et al 2017 Heredity (2017) 119, 55–63"
#Dgen.1 <- dist.genpop(gp,method=1)
Dgen.2 <- dist.genpop(gp,method=2)
#Dgen.3 <- dist.genpop(gp,method=3)
#Dgen.4 <- dist.genpop(gp,method=4)
#Dgen.5 <- dist.genpop(gp,method=5)

#Get pariwise site distances

#filter to included sites
site_coords.subset = site_coords %>% filter(state.name %in% levels(pop(gl.subset)))
#need to order the df to compare distances correctly
site_order = levels(pop(gl.subset))
site_coords.subset.order = site_coords.subset[match(site_order, site_coords.subset$state.name), ]
rownames(site_coords.subset.order) = site_coords.subset.order$state.name

Dgeo <- dist(dismo::Mercator(site_coords.subset.order[,c("lon", "lat")]))

#ibd.1 <- mantel.randtest(Dgen.1,Dgeo)
ibd.2 <- mantel.randtest(Dgen.2,Dgeo) # r =  0.1889 , P = 0.143
#ibd.3 <- mantel.randtest(Dgen.3,Dgeo)
#ibd.4 <- mantel.randtest(Dgen.4,Dgeo)
#ibd.5 <- mantel.randtest(Dgen.5,Dgeo)

#ibd.1.log <- mantel.randtest(Dgen.1,log(Dgeo))
ibd.2.log <- mantel.randtest(Dgen.2,log(Dgeo)) # r = 0.24 , P = 0.11
#ibd.3.log <- mantel.randtest(Dgen.3,log(Dgeo))
#ibd.4.log <- mantel.randtest(Dgen.4,log(Dgeo))
#ibd.5.log <- mantel.randtest(Dgen.5,log(Dgeo))


Dgeo.long = subset(reshape2::melt(Dgeo %>% as.matrix), value != 0)
Dgen.long = subset(reshape2::melt(Dgen.2 %>% as.matrix), value != 0)

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
p2
dev.off()

######################################
#Running without VA (very young site)#
######################################

#SOME SITES HAVE SMALL SAMPLE SIZE
#Set min sample size to either 3 or 4
low_n_states = ( (sample_metadata %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 4) )$state.name
rm.ind.list = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(pop.gl %in% c(low_n_states))

#rveerse filtering of above, but then randomly select four samples per site
low_n_states = ( (sample_metadata %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 4) )$state.name
keep.ind.list = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(!pop.gl %in% c(low_n_states, "VA"))

keep.ind.rand = list()
pops_incl = keep.ind.list$pop.gl %>% as.character %>% unique

for(i in 1:length(pops_incl)){
    temp = keep.ind.list %>% filter(pop.gl == pops_incl[i])
    keep.ind.rand[[pops_incl[i]]] = temp[sample(1:nrow(temp), 4), ] # random sample of four rows
}

keep.ind.rand.df = bind_rows(keep.ind.rand)

#keep.ind.rand.df = structure(list(gl.ind.names = c("NG34", "NG55", "NG21", "NG16",
#"NG135", "NG142", "NG140", "NG126", "NG105", "NG115", "NG133",
#"NG104", "NG12", "NG9", "NG11", "NG28", "NG128", "NG114", "NG159",
#"NG134", "NG137", "NG121", "NG136", "NG119", "NG154", "NG35",
#"NG44", "NG40", "NG74", "NG8", "NG30", "NG29", "NG68", "NG72",
#"NG75", "NG71", "NG56", "NG38", "NG48", "NG37", "NG63", "NG66",
#"NG59", "NG62"), pop.gl = structure(c(12L, 12L, 12L, 12L, 6L,
#6L, 6L, 6L, 13L, 13L, 13L, 13L, 17L, 17L, 17L, 17L, 15L, 15L,
#15L, 15L, 1L, 1L, 1L, 1L, 8L, 8L, 8L, 8L, 3L, 3L, 3L, 3L, 11L,
#11L, 11L, 11L, 5L, 5L, 5L, 5L, 14L, 14L, 14L, 14L), .Label = c("CT",
#"ME.N", "ME.S", "MI", "NC", "NH.BART", "NH.CCM", "NH.CW", "NH.SCG",
#"NJ", "NY.N", "NY.S", "NY.W", "PA", "PA.W", "VA", "WV"), class = "factor")), row.names = c(NA,
#-44L), class = "data.frame")
##########################################

gl.subset = gl[gl@ind.names %in% keep.ind.rand.df$gl.ind.names]


gl.subset
#calculate Mercator prjected distance
xy <- dismo::Mercator(gl.subset@other$latlong[,c("lon","lat")])

#Need to convert to genpop for adegenet dist.genpop
y = as.data.frame(gl.subset)
gi.subset = df2genind(y + 1, ploidy=1)
#reset pop
gi.subset@pop = gl.subset@pop

gp = genind2genpop(gi.subset)

#There are several distance metrics available
#Method 2 is "Angular distance or Edward's distance" D[CSE]
#Dcse performs well in IBD tests "Sere et al 2017 Heredity (2017) 119, 55–63"
#Dgen.1 <- dist.genpop(gp,method=1)
Dgen.2 <- dist.genpop(gp,method=2)
#Dgen.3 <- dist.genpop(gp,method=3)
#Dgen.4 <- dist.genpop(gp,method=4)
#Dgen.5 <- dist.genpop(gp,method=5)

#Get pariwise site distances

#filter to included sites
site_coords.subset = site_coords %>% filter(state.name %in% levels(pop(gl.subset)))
#need to order the df to compare distances correctly
site_order = levels(pop(gl.subset))
site_coords.subset.order = site_coords.subset[match(site_order, site_coords.subset$state.name), ]
rownames(site_coords.subset.order) = site_coords.subset.order$state.name

Dgeo <- dist(dismo::Mercator(site_coords.subset.order[,c("lon", "lat")]))

#ibd.1 <- mantel.randtest(Dgen.1,Dgeo)
ibd.2 <- mantel.randtest(Dgen.2,Dgeo) # r =  0.55 , P = 0.005
#ibd.3 <- mantel.randtest(Dgen.3,Dgeo)
#ibd.4 <- mantel.randtest(Dgen.4,Dgeo)
#ibd.5 <- mantel.randtest(Dgen.5,Dgeo)

#ibd.1.log <- mantel.randtest(Dgen.1,log(Dgeo))
ibd.2.log <- mantel.randtest(Dgen.2,log(Dgeo)) # r = 0.59 , P = 0.002
#ibd.3.log <- mantel.randtest(Dgen.3,log(Dgeo))
#ibd.4.log <- mantel.randtest(Dgen.4,log(Dgeo))
#ibd.5.log <- mantel.randtest(Dgen.5,log(Dgeo))


Dgeo.long = subset(reshape2::melt(Dgeo %>% as.matrix), value != 0)
Dgen.long = subset(reshape2::melt(Dgen.2 %>% as.matrix), value != 0)

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


pdf("figures/IBD.Dcse_test.no_VA.pdf", width = 9, height = 5)
p1
p2
dev.off()
