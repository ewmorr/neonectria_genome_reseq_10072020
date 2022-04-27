require(ggplot2)
require(dplyr)
source("~/ggplot_theme.txt")

sample_IDs.new = read.table("bbduk_trimmed_counts/sample_IDs_03312022.txt")
sample_IDs.old = read.table("bbduk_trimmed_counts/sample_IDs.txt")
read_counts.new = read.table("bbduk_trimmed_counts/bbduk_filtered_read_count_03312022.txt")
read_counts.old = read.table("bbduk_trimmed_counts/bbduk_filtered_read_count.txt")


reads.new = cbind(sample_IDs.new, read_counts.new)
reads.new = reads.new[-64,]
reads.old = cbind(sample_IDs.old, read_counts.old)

reads.new$run = "03312022"
reads.old$run = "11012020"

colnames(reads.new) = c("sample", "reads", "run")
colnames(reads.old) = c("sample", "reads", "run")


reads_all = rbind(reads.new, reads.old)
quantile(reads.new$reads, probs = c(0.1, 0.5, 0.9))
quantile(reads.old$reads, probs = c(0.1, 0.5, 0.9))
quantile(reads.old$reads, probs = c(0.25, 0.5, 0.75))

reads.new %>% filter(reads <= 3590578) #median of old
reads.new %>% filter(reads <= 3312613) #25th percentile of old
reads.new[order(reads.new$reads),]

pdf("bbduk_trimmed_counts/read_count.pdf", height = 6, width = 10)
ggplot(reads_all, aes(x = reorder(sample, reads), y = reads, color = as.factor(run))) +
geom_point(size = 1) +
scale_y_log10() +
labs(x = "sample", color = "run") +
scale_color_manual(values = cbPalette) +
my_gg_theme +
theme(
axis.text.x = element_blank()
)
dev.off()

ggplot(reads_all %>% filter(run == "03312022"), aes(x = reads)) +
geom_histogram(bins = 10, color = "white") +
scale_x_log10()

ggplot(reads_all %>% filter(run == "11012020"), aes(x = reads)) +
geom_histogram(bins = 10, color = "white") +
scale_x_log10()

ggplot(reads_all, aes(x = reads, fill = as.factor(run))) +
geom_histogram(bins = 10, color = "white") +
scale_x_log10()

ggplot(reads_all, aes(x = reads, fill = as.factor(run))) +
geom_histogram(bins = 5, color = "white") +
scale_x_log10()

ggplot(reads_all, aes(x = reads, fill = as.factor(run))) +
geom_histogram(bins = 20, color = "white") +
scale_x_log10()

ggplot(reads_all, aes(x = reads, fill = as.factor(run))) +
geom_histogram(bins = 50, color = "white") +
scale_x_log10()

ggplot(reads_all %>% filter(run == "03312022"), aes(x = reads)) +
geom_histogram(bins = 5, color = "white")
