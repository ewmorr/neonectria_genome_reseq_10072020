library(PopGenome)
library(dplyr)
library(ggplot2)
source("R_scripts/ggplot_theme.txt")


#READING DATA and checking it
source("R_scripts/make_site_metadata.r")

#Concatenate chromosomes to perform calculations across whole genome

#snp.concat = readData("data/Nf_SPANDx_all_seqs/scaffolds_split_rm_low_n", format="VCF", include.unknown = T)
snp.concat = readData("data/Nf_SPANDx_all_seqs/noINDEL_fasta", format="FASTA", include.unknown = T, FAST = T, big.data = T)
sum(snp.concat@n.biallelic.sites)
sum(snp.concat@n.biallelic.sites) + sum(snp.concat@n.polyallelic.sites)
#116513 SNPs biallelic, no poly

rm(snp.concat)

snp.concat = readData("data/Nf_SPANDx_all_seqs/noINDEL_fasta_rm_low_n", format="FASTA", include.unknown = T, FAST = T, big.data = T)
#takes about 5 minutes. holding 6.2 G RAM

sum(snp.concat@n.biallelic.sites)
#114618
sum(snp.concat@n.biallelic.sites) + sum(snp.concat@n.polyallelic.sites)
#114618
sum(snp.concat@n.sites)
#42805617 (this a little short... total len by quast is 42948211 = 142594 dif)
sum(snp.concat@n.gaps)
#0
sum(snp.concat@n.valid.sites)
#NaN
sum(snp.concat@n.unknowns)
#0

snp.concat = concatenate.regions(snp.concat)

get.sum.data(snp.concat)
show.slots(snp.concat)

sum(snp.concat@n.biallelic.sites)
#114618
sum(snp.concat@n.biallelic.sites) + sum(snp.concat@n.polyallelic.sites)
#114618


#########
#ADD POPULATION INFO to genome object
pops <- get.individuals(snp.concat)[[1]]
sample_metadata.sorted = left_join(data.frame(sample = pops), sample_metadata)
pop.levels = as.character(unique(sample_metadata.sorted$state.name))

pop_levels.sampleID.list = list()

for(i in 1:length(pop.levels)){
  pop_levels.sampleID.list[[ pop.levels[i] ]] = as.character((filter(sample_metadata.sorted, state.name == pop.levels[i]))$sample)
}

site_sample_nums = lapply(pop_levels.sampleID.list, length) %>% unlist

#Join pops to site data
site.info = read.csv("data/sample_metadata/site_info.dur_inf.csv")

site.info.all = left_join(data.frame(state.name = pop.levels), site.info, by = "state.name") %>% 
  left_join(., data.frame(state.name = names(site_sample_nums), site_sample_nums))

snp.concat  <- set.populations(snp.concat, pop_levels.sampleID.list)
snp.concat@populations

############################################
#calculate diversity stats on whole genome
snp.concat <- F_ST.stats(snp.concat) # this does the calculations and
# adds the results to the appropriate slots
snp.concat <- neutrality.stats(snp.concat) # this does the calculations and
# adds the results to the appropriate slots
#snp.concat <- linkage.stats(snp.concat) # this does the calculations and
# adds the results to the appropriate slots
snp.concat <- sweeps.stats(snp.concat) # this does the calculations and #THIS DOESN'T REALLY MAKE SENSE at chromosome level
# adds the results to the appropriate slots
#snp.concat <- diversity.stats(snp.concat, pi = T) #This is now working after filtering out the low sample sites
#This last command overwrites the Pi cal from F_ST.stats andwrites all zeroes...


##############
#IT APPEARS THIS SUCKS UP A TON OF MEMORY, PARTICULARLY THE LINKAGE STATS (MEM DUMPED)
#MAY NEED TO PROCEED WITH IND CHROMOSOME LEVEL OR LEAVE OUT LINKAGE AT whole genome level)
#The other metrics run fine

saveRDS(snp.concat, "data/intermediate_RDS/pops.GENOME.concat.rm_low_n.rds")
snp.concat = readRDS("data/intermediate_RDS/pops.GENOME.concat.rm_low_n.rds")

# Print FST
get.F_ST(snp.concat) 
snp.concat@nucleotide.F_ST

# Print diversities
get.diversity(snp.concat)
get.diversity(snp.concat)[[1]] # pop1 (B)
get.diversity(snp.concat)[[2]] # pop2 (b)
snp.concat@nuc.diversity.within
#From te docs:
#Note: the @nuc.divesrity.within have to be normalized/divided by the total number of nucleotides in a given window/region!

show.slots(snp.concat)

snp.concat@Pi #according to the docs this is not valid with NA values
#snp@Pi

##################
#For each divestity metric it is now added to the genome object and can be accessed by SLOTS and is the calculated by pairwise pop comparisons by row
#Columns are the chromosomes (if not concatenated)

#combining with site names and plotting

site_div = left_join(
    data.frame(
        state.name = pop.levels,
        pop.name = paste("pop", seq(1:12)),
        Pi = snp.concat@nuc.diversity.within[1,]/snp.concat@n.sites,
        Tajima.D = snp.concat@Tajima.D[1,],
        nuc_div_within = snp.concat@nuc.diversity.within[1,],
        n_segregating_sites = snp.concat@n.segregating.sites[1,],
        n_sites = snp.concat@n.sites
    ),
    site.info.all,
    by = "state.name"
)

harmonic_number = function(x){
  harms = vector(mode = "numeric", length = length(x))
  for(u in 1:length(x)){
    harm = 0
    for(i in 1:x[u]){
      harm = harm + (1/i)
    }
    harms[u] = harm
  }
  return(harms)
}


site_div$theta = site_div$n_segregating_sites/site_div$n_sites/harmonic_number(site_div$site_sample_nums)

write.table(site_div, "data/intermediate_RDS/population_level_diversity.rm_low_n.txt", col.names = T, row.names = F, quote = F, sep = "\t")

saveRDS(snp.concat, "data/intermediate_RDS/pops.GENOME.concat.rm_low_n.rds")

plot(Pi ~ duration_infection, data = site_div)
summary(lm(Pi ~ duration_infection, data = site_div))
#no relationship Pi ~ dur_inf

plot(Tajima.D ~ duration_infection, data = site_div )
summary(lm(Tajima.D ~ duration_infection, data = site_div))

plot((nuc_div_within/n_sites)-theta ~ duration_infection, data = site_div )
plot(theta ~ duration_infection, data = site_div )
plot(Pi-theta ~ duration_infection, data = site_div )

plot(Tajima.D ~ duration_infection, data = site_div %>% filter(state.name != "VA") )
summary(lm(Tajima.D ~ duration_infection, data = site_div %>% filter(state.name != "VA") ))
# est = -0.0196, R2 = 0.375, p = 0.0452

plot(nuc_div_within ~ duration_infection, data = site_div) #accrding to manual "have to be divided by the slot @n.sites to obtain diversity per site
summary(lm(nuc_div_within ~ duration_infection, data = site_div))

plot(nuc_div_within/n_sites ~ duration_infection, data = site_div)
summary(lm(nuc_div_within/n_sites ~ duration_infection, data = site_div))

plot(n_segregating_sites ~ duration_infection, data = site_div)
summary(lm(n_segregating_sites ~ duration_infection, data = site_div))

#segregating sites needs to be corrected for sample number
plot(n_segregating_sites/site_sample_nums ~ duration_infection, data = site_div)
summary(lm(n_segregating_sites/site_sample_nums ~ duration_infection, data = site_div))

#PLOTS
p1 = ggplot(site_div, aes(x = duration_infection, y = nuc_div_within/n_sites)) +
geom_point() +
#geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "BBD infection duration (yrs)", y = expression(pi)) +
my_gg_theme

p2 = ggplot(site_div, aes(x = duration_infection, y = Tajima.D)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "BBD infection duration (yrs)", y = "Tajima's D") +
my_gg_theme

p3 = ggplot(site_div %>% filter(state.name != "VA"), aes(x = duration_infection, y = Tajima.D)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "BBD infection duration (yrs)", y = "Tajima's D") +
my_gg_theme



pdf("figures/Pi_v_dur_inf.pdf", width = 6, height = 4)
p1
dev.off()


pdf("figures/TajD_v_dur_inf.pdf", width = 6, height = 4)
p2
dev.off()

pdf("figures/TajD_v_dur_inf.no_VA.pdf", width = 6, height = 4)
p3
dev.off()

