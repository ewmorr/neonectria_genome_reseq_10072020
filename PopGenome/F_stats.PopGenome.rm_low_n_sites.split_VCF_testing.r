require(PopGenome)


snp <- readData("Nf_post_SPANDx/scaffolds_split_rm_low_n", format="VCF", include.unknown = T)

snp.scf_lens = data.frame(scf_ids = names(snp@n.sites), lens = snp@n.sites, stringsAsFactors = F)

chr.1 = readVCF(
    "Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.rm_low_n_sites.vcf.gz",
    include.unknown = T,
    frompos = 1,
    topos = snp.scf_lens$lens[1],
    tid = snp.scf_lens$scf_ids[1],
    approx = T,
    numcols = snp.scf_lens$lens[1]
)

#readVCF is throwing an error on genotype field

#Can make this a list and lapply

chr.1 = readData(
"Nf_post_SPANDx/scf_1",
include.unknown = T,
format = "VCF"
)

chr.2 = readData(
"Nf_post_SPANDx/scf_2",
include.unknown = T,
format = "VCF"
)

#This still throws an error
chr.12 = concatenate.classes(list(chr.1, chr.2))




chr.list = list()

for(i in 1:length(snp.scf_lens$scf_ids)){
    
    chr.list
    
}

##SLIDING WINDOW ON INDIVIDUAL CHR WITH READ DATA WORKS
#Try sliding window before loading other packages and setting populations
snp_sw = sliding.window.transform(chr.1, 10000, 10000, type = 2, whole.data = T)

sum(chr.1@n.biallelic.sites)
sum(chr.1@n.biallelic.sites) + sum(chr.1@n.polyallelic.sites)




#Trying set pop
require(dplyr)
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)

pops <- get.individuals(chr.1)[[1]]
sample_metadata.sorted = left_join(data.frame(sample = pops), sample_metadata)
pop.levels = as.character(unique(sample_metadata.sorted$State))

pop_levels.sampleID.list = list()

for(i in 1:length(pop.levels)){
    pop_levels.sampleID.list[[ pop.levels[i] ]] = as.character((filter(sample_metadata.sorted, State == pop.levels[i]))$sample)
}


chr.1  <- set.populations(chr.1, pop_levels.sampleID.list)
chr.1@populations


snp_sw = sliding.window.transform(chr.1, 10000, 10000, type = 2, whole.data = T)
#Still works after pops set. Try stats

snp_sw <- F_ST.stats(snp_sw, mode = "nucleotide")
snp_sw@nucleotide.F_ST

