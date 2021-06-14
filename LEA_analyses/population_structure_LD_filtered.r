require(LEA)

#This is the read command
#The object takes the path to the geno file
test.geno = vcf2geno("Nf_post_SPANDx/LD_filter/test_dat.vcf")

#reading the full data set causes segfault
#test.geno = vcf2geno("Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf")

Nf.K1-10 = snmf(test.geno, K = 1:10, ploidy = 1, entropy = T,
CPU = 4, project = "new") #four CPUs local. Change for server

#This is throwing a sefault, but probably because of the error "Error: Unknown element '2' in the data file"
Nf.K1-10 = snmf(vcf2geno("~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/test_dat.vcf"), K = 1:10, ploidy = 1, entropy = T,
CPU = 4, project = "new") #four CPUs local. Change for server
