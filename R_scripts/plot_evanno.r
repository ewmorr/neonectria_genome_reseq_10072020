require(ggplot2)
require(gridExtra)
require(gtable)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#evanno = read.table("data/Nf_SPANDx_all_seqs_str_th/evanno_out_K1-10.r_format.txt", header = T, sep = "\t")
evanno = read.table("data/Nd_SPANDx_all_seqs_str_th/evanno_out.K1-8.r_format.txt", header = T, sep = "\t")

p1 = ggplot(evanno, aes(K, Mean.LnP.K)) +
geom_point() +
#geom_line() +
geom_errorbar(data = evanno, aes(x = K, ymin = Mean.LnP.K - Stdev.LnP.K, ymax = Mean.LnP.K + Stdev.LnP.K)) +
scale_x_continuous(limits = c(1,10), breaks = c(2,5,10)) +
labs(y = "L(K)") +
my_gg_theme

p2 = ggplot(evanno, aes(K, Ln.prime.K)) +
geom_point() +
#geom_line() +
geom_errorbar(data = evanno, aes(x = K, ymin = Ln.prime.K - Stdev.LnP.K, ymax = Ln.prime.K + Stdev.LnP.K)) +
scale_x_continuous(limits = c(1,10), breaks = c(2,5,10)) +
labs(y = "L'(K)") +
my_gg_theme

p3 = ggplot(evanno, aes(K, Ln.primer.primer.K)) +
geom_point() +
#geom_line() +
geom_errorbar(data = evanno, aes(x = K, ymin = Ln.primer.primer.K - Stdev.LnP.K, ymax = Ln.primer.primer.K + Stdev.LnP.K)) +
scale_x_continuous(limits = c(1,10), breaks = c(2,5,10)) +
labs(y = "L''(K)") +
my_gg_theme

p4 = ggplot(evanno, aes(K, Delta.K)) +
geom_point() +
geom_line() +
#geom_errorbar(data = evanno, aes(x = K, ymin = L.primer.primer.K - Stdev.LnP.K., ymax = L.primer.primer.K + Stdev.LnP.K.)) +
scale_x_continuous(limits = c(1,10), breaks = c(2,5,10)) +
labs(y = expression(paste(Delta,"K"))) +
my_gg_theme

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
gp3<-ggplotGrob(p3)
gp4<-ggplotGrob(p4)

gp13 = rbind(gp1, gp3)
gp24 = rbind(gp2, gp4)

pdf("figures/Nd.evanno.new_seqs.pdf", width = 8, height = 6)
grid.arrange(gp13, gp24, ncol = 2)
dev.off()
