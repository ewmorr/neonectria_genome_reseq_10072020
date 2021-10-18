require(LEA)

#This is the read command
#The object takes the path to the geno file
#Note that the package does not appear to accept VCF v4.2 so we need to first convert to PED format. This is done using plink v1.9 (outside of R)
test.geno = ped2geno("Nf_post_SPANDx/LD_filter/test_dat.ped")
pca(test.geno,scale = T)

#reading the full data set causes segfault
#test.geno = vcf2geno("Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf")

#This is still throwing a segfault on the error of an unrecognized character. The pca function in this package throws a missing data error. Trying tess3r
#Nf.K1-10 = snmf(test.geno, K = 1:10, ploidy = 1, entropy = T, CPU = 4, project = "new") #four CPUs local. Change for server


