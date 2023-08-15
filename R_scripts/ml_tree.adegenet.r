require(adegenet)
require(ape)
require(phangorn)

# read SNP multi seq aligment fasta
seqs = fasta2DNAbin(file = "data/Nf_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.noINDEL.fasta")

# https://adegenet.r-forge.r-project.org/files/PRstats/practical-introphylo.1.0.pdf
# convert to phangorn format
seqs.phyDat = as.phyDat(seqs)

#initialize tree with NJ (using original DNAbinobject, not phangorn format)
#
tre.ini = nj(dist.dna(seqs, pairwise.deletion = T)) #must use pairwise deletion because of high number of gaps
plot(tre.ini)

#optimization
fit.ini = pml(tre.ini, seqs.phyDat, k = 4)

#this takes quite a while to run (like 30 minutes or more)
fit = optim.pml(fit.ini, optNni=T, optBf = T, optQ = T, optGamma = T)

# computed
fit

fit$tree
plot(fit$tree)
anova(fit.ini,fit)

plot(tre)

#write the ml tree before plotting so don't have to rerun

saveRDS(fit, "data/intermediate_RDS/Nf_ML_tree.rds")
