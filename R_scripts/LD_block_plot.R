require(ggplot2)
require(dplyr)

LD_blocks = read.table("data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_block.blocks.det", header = T)

p <- ggplot(LD_blocks,aes(x=KB))+
  geom_density(size=0.5,colour="grey40")+
  labs(x="LD block length (Kb)",y="Density")+
  theme_bw()
p

median(LD_blocks$KB)
# 0.1535 KB
mean(LD_blocks$KB)
# 0.6506

p <- ggplot(LD_blocks %>% filter(KB < 10),aes(x=KB))+
  geom_density(size=0.5,colour="grey40")+
  labs(x="LD block length (Kb)",y="Density")+
  geom_vline(xintercept = median(LD_blocks$KB), color = "red") +
  geom_vline(xintercept = mean(LD_blocks$KB), color = "blue") +
  theme_bw()
p

p <- ggplot(LD_blocks,aes(x=KB))+
  geom_histogram(binwidth = 0.25)+
  labs(x="LD block length (Kb)",y="Count (250 bp bins)")+
  geom_vline(xintercept = median(LD_blocks$KB), color = "red", alpha = 0.5) +
  geom_vline(xintercept = mean(LD_blocks$KB), color = "blue", alpha = 0.5) +
  theme_bw()
p

pdf("figures/LD_block_distribution.pdf")
p
dev.off()
