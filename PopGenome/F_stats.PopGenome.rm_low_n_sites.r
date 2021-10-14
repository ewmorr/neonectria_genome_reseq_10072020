require(PopGenome)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")
require(dplyr)

#READING DATA and checking it
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)


snp <- readData("Nf_post_SPANDx/scaffolds_split_rm_low_n", format="VCF", include.unknown = T)

get.sum.data(snp)
show.slots(snp)

sum(snp@n.biallelic.sites)
sum(snp@n.biallelic.sites) + sum(snp@n.polyallelic.sites)

#ADD POPULATION INFO to genome object
pops <- get.individuals(snp)[[1]]
sample_metadata.sorted = left_join(data.frame(sample = pops), sample_metadata)
pop.levels = as.character(unique(sample_metadata.sorted$State))

pop_levels.sampleID.list = list()

for(i in 1:length(pop.levels)){
    pop_levels.sampleID.list[[ pop.levels[i] ]] = as.character((filter(sample_metadata.sorted, State == pop.levels[i]))$sample)
}


snp  <- set.populations(snp, pop_levels.sampleID.list)
snp@populations

#Join pops to site data
site.info = read.csv("sample_metadata/site_info.csv")
site.GDD = read.table("sample_metadata/site_info.GDD.txt", header = T)
site.climate = read.table("sample_metadata/sites_climate.txt", header = T)
site.coords = read.table("sample_metadata/site_coords.txt", header = T)

site.info.all = left_join(data.frame(state = pop.levels), site.info, by = "state") %>%
    left_join(., site.GDD, by = "Site") %>%
    left_join(., site.climate, by = "Site") #%>%
#    left_join(., site.coords, by = "Site")



#DIVERSITY calcs
# Diversities and FST (by scaffold)
snp <- F_ST.stats(snp) # this does the calculations and
# adds the results to the appropriate slots
snp <- neutrality.stats(snp) # this does the calculations and
# adds the results to the appropriate slots
snp <- linkage.stats(snp) # this does the calculations and
# adds the results to the appropriate slots
snp <- sweeps.stats(snp) # this does the calculations and #THIS DOESN'T REALLY MAKE SENSE at chromosome level
# adds the results to the appropriate slots
#snp <- diversity.stats(snp, pi = T) #This is now working after filtering out the low sample sites

# Print FST
get.F_ST(snp) # each line is a scaffold
snp@nucleotide.F_ST

# Print diversities
get.diversity(snp)
get.diversity(snp)[[1]] # pop1 (B)
get.diversity(snp)[[2]] # pop2 (b)
snp@nuc.diversity.within

#At this point it would make sense to filter out chromosomes where there is no or very low within pop nucleotide diversity
snp@n.sites
#E.G., scaffold 7 (len 278338) , 14 (len 7), 16 (24493), 18 (27481)

#or can use conactenate command
snp.concat = readData("Nf_post_SPANDx/scaffolds_split_rm_low_n", format="VCF", include.unknown = T)
snp.concat = concatenate.regions(snp.concat)

get.sum.data(snp.concat)
show.slots(snp.concat)

sum(snp@n.biallelic.sites)
sum(snp@n.biallelic.sites) + sum(snp@n.polyallelic.sites)

sum(snp.concat@n.biallelic.sites)
sum(snp.concat@n.biallelic.sites) + sum(snp.concat@n.polyallelic.sites)

#The concatenate.regions functoin throws a warning but seems to have worked

#########
# set pop
snp.concat  <- set.populations(snp.concat, pop_levels.sampleID.list)
snp.concat@populations

############################################
#Recalculate diversity stats on whole genome
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

# Print FST
get.F_ST(snp.concat) # each line is a scaffold
snp.concat@nucleotide.F_ST

# Print diversities
get.diversity(snp.concat)
get.diversity(snp.concat)[[1]] # pop1 (B)
get.diversity(snp.concat)[[2]] # pop2 (b)
snp.concat@nuc.diversity.within

show.slots(snp.concat)

snp.concat@Pi
snp@Pi

##################
#For each divestity metric it is now added to the genome object and can be accessed by SLOTS and is the calculated by pairwise pop comparisons by row
#Columns are the chromosomes (if not concatenated)

#E.g.
show.slots(snp)
snp@nuc.F_ST.pairwise
snp@Pi #Pi is already nomalized by sample within the calculation

snp@haplotype.F_ST
snp@hap.F_ST.pairwise
snp@Tajima.D

#combining with site names and plotting

site_div = left_join(
    data.frame(
        state = pop.levels,
        pop.name = paste("pop", seq(1:7)),
        Pi = snp.concat@Pi[1,],
        Tajima.D = snp.concat@Tajima.D[1,],
        nuc_div_within = snp.concat@nuc.diversity.within[1,]
    ),
    site.info.all,
    by = "state"
)

plot(Pi ~ duration_infection, data = site_div)
summary(lm(Pi ~ duration_infection, data = site_div))
plot(Tajima.D ~ duration_infection, data = site_div)
summary(lm(Tajima.D ~ duration_infection, data = site_div))
plot(nuc_div_within ~ duration_infection, data = site_div)
summary(lm(nuc_div_within ~ duration_infection, data = site_div))

#For sliding window calcs may want to NOT concatenate so that regions are not arbitrarily added together
#Note that it looks like theta estimators are only calculated in sliding window

