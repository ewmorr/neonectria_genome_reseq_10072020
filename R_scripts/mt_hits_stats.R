library(dplyr)

#Nd
blast_res = read.table("data/N_ditissima_ref_genome/Nd_fonseca_mt.blast", header = F)
scf_len = read.table("data/N_ditissima_ref_genome/scaffold_lengths.txt", header = F)


colnames(blast_res) = c(
    "contig",
    "dbID",
    "perc_sim",
    "hit_length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "e-val",
    "bitscore"
)

colnames(scf_len) = c("contig", "contig_len")

mt_hits = blast_res %>%
    group_by(contig) %>%
    summarize(tot_hit_len = sum(hit_length), n_hit = n()) %>%
    left_join(., scf_len, by = "contig") %>%
    mutate(perc_contig = tot_hit_len/contig_len*100)

range(mt_hits$perc_contig)
#0.006783116 69.375145383
plot(sort(mt_hits$perc_contig))

write.table(mt_hits, "data/N_ditissima_ref_genome/mt_hits_summary.txt", col.names = T, row.names = F, quote = F)    
# considering anything over 10% as mt 
# 

#Nc
blast_res = read.table("data/N_coccinea_ref_genome/Nc_fonseca_mt.blast", header = F)
scf_len = read.table("data/N_coccinea_ref_genome/scaffold_lengths.txt", header = F)


colnames(blast_res) = c(
    "contig",
    "dbID",
    "perc_sim",
    "hit_length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "e-val",
    "bitscore"
)

colnames(scf_len) = c("contig", "contig_len")

mt_hits = blast_res %>%
    group_by(contig) %>%
    summarize(tot_hit_len = sum(hit_length), n_hit = n()) %>%
    left_join(., scf_len, by = "contig") %>%
    mutate(perc_contig = tot_hit_len/contig_len*100)

range(mt_hits$perc_contig)
#0.0249575 71.8243040
plot(sort(mt_hits$perc_contig))
plot(mt_hits$perc_contig ~ mt_hits$contig_len)

write.table(mt_hits, "data/N_coccinea_ref_genome/mt_hits_summary.txt", col.names = T, row.names = F, quote = F)    
