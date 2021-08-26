require(ggplot)
require(gridExtra)
require(gtable)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

evanno = read.table("Nf_post_SPANDx/LD_filter/structure_threader/bestK/evanno_r_format.txt", header = T)

p1 = ggplot(evanno, aes(K, Mean.LnP.K.)) +
geom_point() +
#geom_line() +
geom_errorbar(data = evanno, aes(x = K, ymin = Mean.LnP.K. - Stdev.LnP.K., ymax = Mean.LnP.K. + Stdev.LnP.K.)) +
scale_x_continuous(limits = c(1,15), breaks = c(2,5,10,15)) +
labs(y = "L(K)") +
my_gg_theme

p2 = ggplot(evanno, aes(K, Ln.prime.K.)) +
geom_point() +
#geom_line() +
geom_errorbar(data = evanno, aes(x = K, ymin = Ln.prime.K. - Stdev.LnP.K., ymax = Ln.prime.K. + Stdev.LnP.K.)) +
scale_x_continuous(limits = c(1,15), breaks = c(2,5,10,15)) +
labs(y = "L'(K)") +
my_gg_theme

p3 = ggplot(evanno, aes(K, L.primer.primer.K)) +
geom_point() +
#geom_line() +
geom_errorbar(data = evanno, aes(x = K, ymin = L.primer.primer.K - Stdev.LnP.K., ymax = L.primer.primer.K + Stdev.LnP.K.)) +
scale_x_continuous(limits = c(1,15), breaks = c(2,5,10,15)) +
labs(y = "L''(K)") +
my_gg_theme

p4 = ggplot(evanno, aes(K, Delta.K)) +
geom_point() +
geom_line() +
#geom_errorbar(data = evanno, aes(x = K, ymin = L.primer.primer.K - Stdev.LnP.K., ymax = L.primer.primer.K + Stdev.LnP.K.)) +
scale_x_continuous(limits = c(1,15), breaks = c(2,5,10,15)) +
labs(y = expression(paste(Delta,"K"))) +
my_gg_theme

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
gp3<-ggplotGrob(p3)
gp4<-ggplotGrob(p4)

gp13 = rbind(gp1, gp3)
gp24 = rbind(gp2, gp4)

pdf("figures/Nf.evanno.pdf", width = 8, height = 6)
grid.arrange(gp13, gp24, ncol = 2)
dev.off()
