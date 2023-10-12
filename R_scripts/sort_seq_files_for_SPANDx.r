library(dplyr)

#setwd("~/GARNAS_neonectria_genome_reseq_10072020")
seq_files_dat = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")

seq_files_dat.Nf = seq_files_dat |> filter(Spp_based_on_ITS_map == "Nf")
seq_files_dat.Nd = seq_files_dat |> filter(Spp_based_on_ITS_map == "Nd")

Nf_samples_to_exclude = c("NG106", "NG108", "NG57", "NG111", "NG160") #NG106 and NG108 were sequenced on first run and had duplicates on the second; we retain the samples with higher output (second run). See README.SNP_calling.md for more info. Samples NG57, NG111, NG160 were resequenced on the third run and we will retain these samples due to poor SNP calling success on first rounds.
n_of_samples.Nf = seq_files_dat.Nf |> group_by(sample) |> summarize(n = n())
n_of_samples.Nd = seq_files_dat.Nd |> filter(Sequence_label != "NG20") |> group_by(sample) |> summarize(n = n()) #filtering out NG20 because it was resequenced on the second run

Nf_dups = ( n_of_samples.Nf |> filter(n > 1) )$sample
Nd_dups = ( n_of_samples.Nd |> filter(n > 1) )$sample

Nf_singles = ( n_of_samples.Nf |> filter(n == 1) )$sample
Nd_singles = ( n_of_samples.Nd |> filter(n == 1) )$sample

Nf_first_set = ( seq_files_dat %>% filter(sample %in% Nf_singles & set == "first") )$Sequence_label
Nf_second_set = ( seq_files_dat %>% filter(sample %in% Nf_singles & set == "second") )$Sequence_label
Nf_second_plus = ( seq_files_dat %>% filter(sample %in% Nf_singles & set == "second_plus") )$Sequence_label

Nd_first_set = ( seq_files_dat %>% filter(sample %in% Nd_singles & set == "first") )$Sequence_label
Nd_second_set = ( seq_files_dat %>% filter(sample %in% Nd_singles & set == "second") )$Sequence_label
Nd_second_plus = ( seq_files_dat %>% filter(sample %in% Nd_singles & set == "second_plus") )$Sequence_label

write.table(Nf_first_set, "sample_metadata/sample_IDs.Nf.10072020.txt", quote = F, col.names = F, row.names = F)
write.table(Nf_second_set, "sample_metadata/sample_IDs.Nf.03312022.txt", quote = F, col.names = F, row.names = F)
write.table(Nf_second_plus, "sample_metadata/sample_IDs.Nf.03312022_adtl_reads.txt", quote = F, col.names = F, row.names = F)
write.table(Nd_first_set, "sample_metadata/sample_IDs.Nd.10072020.txt", quote = F, col.names = F, row.names = F)
write.table(Nd_second_set, "sample_metadata/sample_IDs.Nd.03312022.txt", quote = F, col.names = F, row.names = F)
write.table(Nd_second_plus, "sample_metadata/sample_IDs.Nd.03312022_adtl_reads.txt", quote = F, col.names = F, row.names = F)


seq_files_dat %>% filter(sample %in% Nf_dups)

write.table(seq_files_dat %>% filter(sample %in% Nf_dups), "sample_metadata/sample_IDs.Nf.duplicate_samples_to_cat.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#the ony Nd dup is MI1 1.1.1 which was sequenced on first round and second round inluding additional reads. We'll just use the reads from the second round as the first round sequenced poorly. Filtering out NG20 above
seq_files_dat %>% filter(sample %in% Nd_dups)

