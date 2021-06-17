require(pcadapt)
path_to_file <- "Nf_post_SPANDx/LD_filter/test_dat.bed"
filename <- read.pcadapt(path_to_file, type = "bed") #This looks like prefers conversion to BED format first

#To check bed file

bed_mat = pcadapt::bed2matrix(filename)

x <- pcadapt(input = filename, K = 20, ploidy = 1)
pca(filename, scale = T)
