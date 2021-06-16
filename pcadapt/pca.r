require(pcadapt)
path_to_file <- "Nf_post_SPANDx/LD_filter/test_dat.ped"
filename <- read.pcadapt(path_to_file, type = "ped") #This looks like prefers conversion to BED format first

x <- pcadapt(input = filename, K = 20, ploidy = 1)
