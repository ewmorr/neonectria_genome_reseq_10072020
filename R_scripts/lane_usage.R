require(dplyr)

cov = read.table("data/coverage/Nf_all_seqs/coverage_by_sample.dedup_bam.txt")
colnames(cov) = c("sample", "avg_cov")

cov$sample = as.numeric(gsub("NG", "", cov$sample))

cov.2nd_set = cov %>% filter(sample > 100)

sum(cov.2nd_set$avg_cov * 42749761)/10^9

counts = read.table("data/bbduk_trimmed_counts/bbduk_filtered_read_count_03312022.txt")
sum(counts[1:63,])
(sum(counts[1:63,])/2)/(400*10^6) #800*10^6 reads per flow cell, half that per lane. The number of reads is total so divide by 2 for pairs

