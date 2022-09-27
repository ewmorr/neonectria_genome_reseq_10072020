require(dplyr)
require(ggplot2)
require(stringr)

source("R_scripts/ggplot_theme.txt")

r_tab = read.table("data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_decay.ld.test", header = T)

#calculate SNP dist
r_tab$distance = r_tab$BP_B - r_tab$BP_A

#calculate means
r_tab$distc <- cut(r_tab$distance,breaks=seq(from=0,to=ceiling(max(r_tab$distance)/1000)*1000,by=1000) )


r_tab_avg <- r_tab %>% group_by(distc, CHR_A) %>% summarise(mean=mean(R2),median=median(R2))

#calc mean distance
r_tab_avg <- r_tab_avg %>% mutate(
    start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
    end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
    mid=start+((end-start)/2)
    )

head(r_tab)

p1 = ggplot(r_tab, aes(x = distance/10^3, y = R2)) +
  geom_point() +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

p2 = ggplot(r_tab, aes(x = distance/10^3, y = R2, color = CHR_A)) +
  geom_point() +
  scale_color_manual(values = c25) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

p3 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean)) +
  geom_point() +
  geom_line() +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

p4 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean, color = CHR_A)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c25) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

p5 = ggplot() +
  geom_point(data = r_tab, aes(x = distance/10^3, y = R2)) +
  geom_line(data = r_tab_avg, aes(x = mid/10^3, y = mean) ) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

pdf("figures/LD_decay.test.pdf", width = 16, height = 6)
p1
p2
p3
p4
p5
dev.off()

