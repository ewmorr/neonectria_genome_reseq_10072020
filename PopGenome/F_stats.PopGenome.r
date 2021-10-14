require(PopGenome)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")
require(dplyr)

#READING DATA and checking it
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)

snp <- readData("Nf_post_SPANDx/scaffolds_split", format="VCF", include.unknown = T)

get.sum.data(snp)
show.slots(snp)

sum(snp@n.biallelic.sites)
sum(snp@n.biallelic.sites) + sum(snp@n.polyallelic.sites)
#ADD POPULATION INFO to genome object
pops <- get.individuals(snp)[[1]]
sample_metadata.sorted = left_join(data.frame(sample = pops), sample_metadata)
pop.levels = as.character(levels(sample_metadata.sorted$State))

pop_levels.sampleID.list = list()

for(i in 1:length(pop.levels)){
    pop_levels.sampleID.list[[ pop.levels[i] ]] = as.character((filter(sample_metadata.sorted, State == pop.levels[i]))$sample)
}


snp  <- set.populations(snp, pop_levels.sampleID.list)
snp@populations


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
snp <- diversity.stats(snp, pi = T)

# Print FST
get.F_ST(snp) # each line is a scaffold
snp@nucleotide.F_ST

# Print diversities
get.diversity(snp)
get.diversity(snp)[[1]] # pop1 (B)
get.diversity(snp)[[2]] # pop2 (b)
snp@nuc.diversity.within


##################
#For each divestity metric it is now added to the genome object and can be accessed by SLOTS and is the calculated by pairwise pop comparisons by row
#Columns are the chromosomes
#can then avg across columns but CAREFUL not all columns are populated
#This is likely beacuase of MAT related genes/scaffolds as it looks like column 16 and 18 are only with NA and they have v high propotions of NA (in pairwise fisher) and are relatively small (20kbp each)

#E.g.
show.slots(snp)
snp@nuc.F_ST.pairwise
snp@Pi #FOR THIS NEED TO TAKE INTO ACCOUNT NUMBER OF INDS PER POP (AND FOR ALL OTHER METRICS????)
#EG POP1 (i.e., ME.N with 1 sample, shows up as 0 diversity (duh)
snp@haplotype.F_ST
snp@hap.F_ST.pairwise
snp@Tajima.D
#The pops 1 (ME.N), 3 (MI), and 6 (NJ) are not be calculated, likely bc of low sample numbers (1, 2, and 3 samples, respectively). Should chuck these for now.
#scaffolds 14, 16, 18 throwing lots of NAs or zeroes so should exclude for all calcs probably

pop_size = unlist(lapply(pop_levels.sampleID.list, length))

apply(snp@Pi, 2, ave)[1,]/pop_size
#May not want to normalize by pop size because this is giving some pretty funky numbers
#then regress this with infection duration


#Because of issues with low sample number, we will filter the VCF with BCFTOOLS and then rerun the PopGenome analyses using a new script `F_stats.PopGenome.exclude_low_n_pops.R`

write.table(
    c(pop_levels.sampleID.list[["ME.S"]], pop_levels.sampleID.list[["NC"]], pop_levels.sampleID.list[["NH"]],  pop_levels.sampleID.list[["NY.N"]], pop_levels.sampleID.list[["NY.S"]], pop_levels.sampleID.list[["PA"]], pop_levels.sampleID.list[["WV"]]),
    file = "Nf_post_SPANDx/low_n_samples.txt",
    row.names = F,
    quote = F,
    col.names = F
)


