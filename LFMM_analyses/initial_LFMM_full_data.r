require(lfmm)
require(dplyr)

#The genotype data can simply be read in as a matrix (according the docs)
#OR can try loading LEA and using readLfmm()
Y = as.matrix(read.table("Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode01missing9.lfmm", header = F))
require(LEA)
y = read.lfmm("Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode01missing9.lfmm")
#Both seem to give the same output

pc = prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(7,pc$sdev[7]^2, type = "h", lwd = 3, col = "blue")


#read env data. Needs to be sorted by PED .sampleIDs info
sampleIDs = read.table("Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode01missing9.sampleIDs", header = F)
colnames(sampleIDs) = "sample"
sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)

sample_metadata.filtered = left_join(sampleIDs, sample_metadata)
colnames(sample_metadata.filtered)[ncol(sample_metadata.filtered)] = "state"

#Join pops to site data
site.info = read.csv("sample_metadata/site_info.csv")
site.GDD = read.table("sample_metadata/site_info.GDD.txt", header = T)
site.climate = read.table("sample_metadata/sites_climate.txt", header = T)
site.coords = read.table("sample_metadata/site_coords.txt", header = T)

sample_metadata.site_info = left_join(sample_metadata.filtered, site.info, by = "state") %>%
    left_join(., site.GDD, by = "Site") %>%
    left_join(., site.climate, by = "Site") %>%
    left_join(., site.coords, by = "Site")

