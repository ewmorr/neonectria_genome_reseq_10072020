require(pcadapt)
path_to_file <- "Nf_post_SPANDx/LD_filter/test_dat.vcf"
filename <- read.pcadapt(path_to_file, type = "vcf") #This looks like prefers conversion to BED format first


