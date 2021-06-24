require(pcadapt)
makeBed <- "Nf_post_SPANDx/LD_filter/plink_conversion_test/test_dat.make-bed.bed"
makeBed.file <- read.pcadapt(makeBed, type = "bed")
bed_mat.makeBed = pcadapt::bed2matrix(makeBed.file)


x <- pcadapt(input = makeBed.file, K = 15, ploidy = 2) #K must be less than the number of species

bed_mat.makeBed.mat = pcadapt::bed2matrix(makeBed.file)

plot(x, option = "screeplot")

plot(x, option = "scores")
