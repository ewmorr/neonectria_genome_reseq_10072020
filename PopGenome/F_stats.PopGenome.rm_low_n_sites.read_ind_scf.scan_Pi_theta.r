require(PopGenome)
require(dplyr)
require(tidyr)
require(ggplot2)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")


#Read scf into list
scf_list = list()
scf_dirs = paste("scf", seq(1,20,1), sep = "_")

for(i in 1:length(scf_dirs)){
    scf_list[[scf_dirs[i]]] = readData(
        paste("Nf_post_SPANDx/", scf_dirs[i], sep = ""),
        include.unknown = T,
        format = "VCF"
    )
}

#Trying set pop

sample_metadata = read.table("sample_metadata/sample_metadata.Nf.txt", header = T)

pops <- get.individuals(scf_list[["scf_1"]])[[1]]
sample_metadata.sorted = left_join(data.frame(sample = pops), sample_metadata)
pop.levels = as.character(unique(sample_metadata.sorted$State))

pop_levels.sampleID.list = list()

for(i in 1:length(pop.levels)){
    pop_levels.sampleID.list[[ pop.levels[i] ]] = as.character((filter(sample_metadata.sorted, State == pop.levels[i]))$sample)
}

scf_list.pops = lapply(scf_list, set.populations, pop_levels.sampleID.list)

for(i in 1:length(scf_dirs)){
    print(scf_dirs[i])
    print(scf_list.pops[[scf_dirs[i]]]@populations)
}

#Read additional metadata
#Join pops to site data
site.info = read.csv("sample_metadata/site_info.csv")
site.GDD = read.table("sample_metadata/site_info.GDD.txt", header = T)
site.climate = read.table("sample_metadata/sites_climate.txt", header = T)
site.coords = read.table("sample_metadata/site_coords.txt", header = T)

site.info.all = left_join(data.frame(state = pop.levels), site.info, by = "state") %>%
left_join(., site.GDD, by = "Site") %>%
left_join(., site.climate, by = "Site")


scf_stats_list = list()
# get n.sites
for(i in 1:length(scf_dirs)){
    scf_stats_list[[scf_dirs[i]]] = scf_list.pops[[scf_dirs[i]]]@n.sites
}

scf_names = do.call(rbind.data.frame, lapply(scf_stats_list, names))
colnames(scf_names) = "scf_name"
scf_names$scf = scf_dirs

scf_stats.nSites = do.call(rbind.data.frame, scf_stats_list)
#n.sites by scf

    scf_list.pops[[scf_dirs[i]]]@theta_Watterson


#combining with site names and plotting
site_info_pops = left_join(
data.frame(
state = pop.levels,
pop.name = paste("pop", seq(1:7))
),
site.info.all[site.info.all$Site != "ADN2" & site.info.all$Site != "ADS2", ],
by = "state"
)

######################
#Trying sliding window

#first run sliding window with NO pops to get dataset wide estimates
#Rm element 14
scf_list.no_pops.rm = scf_list
scf_list.no_pops.rm[["scf_14"]] = NULL
scf_dirs.rm = scf_dirs[! scf_dirs == "scf_14"]
scf_names.rm = scf_names[-14,]

scf_list.no_pops.sw = lapply(scf_list.no_pops.rm, sliding.window.transform, 10000, 10000, type = 2, whole.data = T)

#Rerun stats
scf_list.no_pops.sw = lapply(scf_list.no_pops.sw, F_ST.stats)
scf_list.no_pops.sw = lapply(scf_list.no_pops.sw, neutrality.stats)
scf_list.no_pops.sw = lapply(scf_list.no_pops.sw, linkage.stats)
scf_list.no_pops.sw = lapply(scf_list.no_pops.sw, sweeps.stats)
scf_list.no_pops.sw = lapply(scf_list.no_pops.sw, diversity.stats)

scf_list.no_pops.sw$scf_1@region.names

#mean position extract
scf.pos <- sapply(scf_list.no_pops.sw$scf_1@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    val <- mean(as.numeric(split))
    return(val)
})

#scf list of positions
genome_pos.list = list()

for(i in 1:length(scf_dirs.rm)){
    genome_pos.list[[scf_dirs.rm[i]]] = sapply(
    scf_list.no_pops.sw[[scf_dirs.rm[i]]]@region.names,
    function(x){
        split <- strsplit(x," ")[[1]][c(1,3)]
        val <- mean(as.numeric(split))
        return(val)
    }
    )
}

scf_list.no_pops.sw$scf_1@nucleotide.F_ST
scf_list.no_pops.sw$scf_2@Tajima.D %>% head
scf_list.no_pops.sw$scf_1@theta_Watterson %>% head
scf_list.no_pops.sw$scf_1@Pi

#Some dplyr BS to format tables

#loop list
######
#Theta
theta_sw_list = list()

for(i in 1:length(scf_dirs.rm)){
    
    theta_sw_list[[ scf_dirs.rm[i] ]] = data.frame(
    scf_list.no_pops.sw[[i]]@theta_Watterson,
    stringsAsFactors = F
    )
    theta_sw_list[[ scf_dirs.rm[i] ]]$pos = unname(genome_pos.list[[i]])
    theta_sw_list[[ scf_dirs.rm[i] ]]$scf = as.character(scf_names.rm[i,1])
    #print(as.character(scf_names.rm[i,1]))
}

theta_sw.df = do.call(rbind.data.frame, theta_sw_list)

#######
#Pi
pi_sw_list = list()

for(i in 1:length(scf_dirs.rm)){
    
    pi_sw_list[[ scf_dirs.rm[i] ]] = data.frame(
    scf_list.no_pops.sw[[i]]@nuc.diversity.within/10000,
    stringsAsFactors = F
    )
    pi_sw_list[[ scf_dirs.rm[i] ]]$pos = unname(genome_pos.list[[i]])
    pi_sw_list[[ scf_dirs.rm[i] ]]$scf = as.character(scf_names.rm[i,1])
    #print(as.character(scf_names.rm[i,1]))
}

pi_sw.df = do.call(rbind.data.frame, pi_sw_list)

#######
#Tajima.D
tajD_sw_list = list()

for(i in 1:length(scf_dirs.rm)){
    
    tajD_sw_list[[ scf_dirs.rm[i] ]] = data.frame(
    scf_list.no_pops.sw[[i]]@Tajima.D,
    stringsAsFactors = F
    )
    tajD_sw_list[[ scf_dirs.rm[i] ]]$pos = unname(genome_pos.list[[i]])
    tajD_sw_list[[ scf_dirs.rm[i] ]]$scf = as.character(scf_names.rm[i,1])
    #print(as.character(scf_names.rm[i,1]))
}

tajD_sw.df = do.call(rbind.data.frame, tajD_sw_list)

####
#Combine allto list

colnames(theta_sw.df) = c("theta", "pos", "scf")
colnames(pi_sw.df) = c("pi", "pos", "scf")
colnames(tajD_sw.df) = c("tajD", "pos", "scf")

stats_sw.no_pop.df = left_join(theta_sw.df, pi_sw.df) %>% left_join(., tajD_sw.df)

#Add lengths to filter scf by len
scf_lens = read.table("Nf_post_SPANDx/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scf", "length")
stats_sw.no_pop.df = left_join(stats_sw.no_pop.df, scf_lens)
stats_sw.no_pop.df$scf %>% unique

stats_sw.no_pop.df.gt100kb = stats_sw.no_pop.df %>% filter(length > 100000)
#Outlier tests
# identify the 95% percentile
#theta
my_threshold <- quantile(stats_sw.no_pop.df.gt100kb$theta, 0.95, na.rm = T)
# make an outlier column in the data.frame
stats_sw.no_pop.df.gt100kb <- stats_sw.no_pop.df.gt100kb %>% mutate(outlier.theta = ifelse(theta > my_threshold, "outlier", "background"))
stats_sw.no_pop.df.gt100kb %>% group_by(outlier.theta) %>% tally()

#Number of outliers

#pi
my_threshold <- quantile(stats_sw.no_pop.df.gt100kb$pi, 0.95, na.rm = T)
# make an outlier column in the data.frame
stats_sw.no_pop.df.gt100kb <- stats_sw.no_pop.df.gt100kb %>% mutate(outlier.pi = ifelse(pi > my_threshold, "outlier", "background"))
#Number of outliers
stats_sw.no_pop.df.gt100kb %>% group_by(outlier.pi) %>% tally()

#tajD
stats_sw.no_pop.df.gt100kb %>% group_by(outlier.tajD) %>% tally()
my_threshold <- quantile(stats_sw.no_pop.df.gt100kb$tajD, c(0.025, 0.975), na.rm = T)
# make an outlier column in the data.frame
stats_sw.no_pop.df.gt100kb <- stats_sw.no_pop.df.gt100kb %>% mutate(outlier.tajD = ifelse(tajD > my_threshold[2], "outlier", ifelse(tajD < my_threshold[1], "outlier", "background")))

#Number of outliers
stats_sw.no_pop.df.gt100kb %>% group_by(outlier.tajD) %>% tally()

ggplot(stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = theta, color = outlier.theta)) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
#geom_line(color = "black", alpha = 0.5) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(theta)) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

ggplot(stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = pi, color = outlier.pi)) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
#geom_line(color = "black", alpha = 0.5) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(Pi)) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)


ggplot(stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = tajD, color = outlier.tajD)) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
#geom_line(color = "black", alpha = 0.5) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "Tajima's D") +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)


pv.hdd.with_pos = read.table("Nf_LFMM_tables/hdd4_lfmm.txt", header = T)
pv.ft.with_pos = read.table("Nf_LFMM_tables/freezeThaw_lfmm.txt", header = T)
pv.ppt.with_pos = read.table("Nf_LFMM_tables/ppt_lfmm.txt", header = T)
pv.dur_inf.with_pos = read.table("Nf_LFMM_tables/dur_inf_lfmm.txt", header = T)

#need same column names on facet
colnames(pv.hdd.with_pos) = c("calibrated.p", "effect_size", "scf", "pos", "man.adj.p",
"length", "outlier", "FDR.p", "FDR.sig", "FDR.p.man", "FDR.sig.man")
colnames(pv.ft.with_pos) = c("calibrated.p", "effect_size", "scf", "pos", "man.adj.p",
"length", "outlier", "FDR.p", "FDR.sig", "FDR.p.man", "FDR.sig.man")
colnames(pv.ppt.with_pos) = c("calibrated.p", "effect_size", "scf", "pos", "man.adj.p",
"length", "outlier", "FDR.p", "FDR.sig", "FDR.p.man", "FDR.sig.man")
colnames(pv.dur_inf.with_pos) = c("calibrated.p", "effect_size", "scf", "pos", "man.adj.p",
"length", "outlier", "FDR.p", "FDR.sig", "FDR.p.man", "FDR.sig.man")

pi_min_max = range(stats_sw.no_pop.df.gt100kb$pi, na.rm = T)
theta_min_max = range(stats_sw.no_pop.df.gt100kb$theta, na.rm = T)
tajD_min_max = range(stats_sw.no_pop.df.gt100kb$tajD, na.rm = T)

p1 = ggplot() +
geom_point(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = pi, color = outlier.pi), alpha = 0.5, size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = pi), color = "black", alpha = 0.5) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_color_manual(values = c("white", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(Pi)) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
#axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8),
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p2 = ggplot() +
geom_linerange(
data = pv.hdd.with_pos %>% filter(FDR.sig == "sig" & length > 100000),
aes(x = pos/10^6, ymin = pi_min_max[1], ymax = pi_min_max[2]),
color = "blue",
alpha = 0.5
) +
geom_point(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = pi, alpha = outlier.pi), color = "black", size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = pi), color = "black", alpha = 0.5) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_alpha_manual(values = c(0, 0.5), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(Pi)) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
#axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8),
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p3 = ggplot() +
geom_point(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = theta, color = outlier.theta), alpha = 0.5, size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = theta), color = "black", alpha = 0.5) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_color_manual(values = c("white", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(theta)) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
#axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8),
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p4 = ggplot() +
geom_linerange(
    data = pv.hdd.with_pos %>% filter(FDR.sig == "sig" & length > 100000),
    aes(x = pos/10^6, ymin = theta_min_max[1], ymax = theta_min_max[2]),
    color = "blue",
    alpha = 0.5
) +
geom_point(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = theta, alpha = outlier.theta), color = "black", size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = theta), color = "black", alpha = 0.5) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_alpha_manual(values = c(0, 0.5), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(theta)) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
#axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8),
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p5 = ggplot() +
geom_point(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = tajD, color = outlier.tajD), alpha = 0.5, size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = tajD), color = "black", alpha = 0.5) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_color_manual(values = c("white", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "Tajima's D") +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

p6 = ggplot() +
geom_linerange(
data = pv.hdd.with_pos %>% filter(FDR.sig == "sig" & length > 100000),
aes(x = pos/10^6, ymin = tajD_min_max[1], ymax = tajD_min_max[2]),
color = "blue",
alpha = 0.5
) +
geom_point(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = tajD, alpha = outlier.tajD), color = "black", size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = tajD), color = "black", alpha = 0.5) +
facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_alpha_manual(values = c(0, 0.5), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "Tajima's D") +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

require(gtable)
require(gridExtra)
require(grid)

plots = list(p1, p3, p5)
grobs <- lapply(plots, ggplotGrob)
g <- do.call(gridExtra::gtable_rbind, grobs)
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1,1,1), "null")

pdf("figures/pi_theta_tajD.sw.pdf", width = 16, height = 8)
grid.newpage()
grid.draw(g)
dev.off()

plots = list(p2, p4, p6)
grobs <- lapply(plots, ggplotGrob)
g <- do.call(gridExtra::gtable_rbind, grobs)
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1,1,1), "null")

pdf("figures/pi_theta_tajD.sw.LFMM.pdf", width = 16, height = 8)
grid.newpage()
grid.draw(g)
dev.off()


p7 = ggplot() +
geom_point(data = stats_sw.no_pop.df.gt100kb %>% filter(scf == "tig00000025_pilon"), aes(x = pos/10^6, y = pi, color = outlier.pi), alpha = 1, size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = pi), color = "black", alpha = 0.5) +
#facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_color_manual(values = c("white", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(Pi)) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
#axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 12),
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p8 = ggplot() +
geom_point(data = stats_sw.no_pop.df.gt100kb %>% filter(scf == "tig00000025_pilon"), aes(x = pos/10^6, y = theta, color = outlier.theta), alpha = 1, size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = theta), color = "black", alpha = 0.5) +
#facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_color_manual(values = c("white", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(theta)) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
#axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 12),
axis.text.x = element_blank(),
axis.title.x = element_blank()
)


p9 = ggplot() +
geom_point(data = stats_sw.no_pop.df.gt100kb %>% filter(scf == "tig00000025_pilon"), aes(x = pos/10^6, y = tajD, color = outlier.tajD), alpha = 1, size = 1) +
geom_line(data = stats_sw.no_pop.df.gt100kb, aes(x = pos/10^6, y = tajD), color = "black", alpha = 0.5) +
#facet_grid(. ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_color_manual(values = c("white", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "Tajima's D") +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 12),
axis.text.y = element_text(size = 12)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

plots = list(p7, p8, p9)
grobs <- lapply(plots, ggplotGrob)
g <- do.call(gridExtra::gtable_rbind, grobs)
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1,1,1), "null")

pdf("figures/pi_theta_tajD.sw.one_tig.pdf", width = 16, height = 8)
grid.newpage()
grid.draw(g)
dev.off()


#########################
#########################
#sliding window with pops
#########################
#########################

scf_list.sw = lapply(scf_list.pops, sliding.window.transform, 10000, 10000, type = 2, whole.data = T)
#This throws an error about attepmting to select one element. Probably the scf with only SNP at pos 7

sw_try = sliding.window.transform(scf_list.pops[["scf_14"]], 10000, 10000, type = 2, whole.data = T)
sw_try = sliding.window.transform(scf_list.pops[["scf_11"]], 10000, 10000, type = 2, whole.data = T)
#That seems to be the problem

#Rm element 14
scf_list.pops.rm = scf_list.pops
scf_list.pops.rm[["scf_14"]] = NULL

#Rerun sliding window
scf_list.sw = lapply(scf_list.pops.rm, sliding.window.transform, 10000, 10000, type = 2, whole.data = T)

#Rerun stats
scf_list.sw = lapply(scf_list.sw, F_ST.stats)
scf_list.sw = lapply(scf_list.sw, neutrality.stats, theta = T)
scf_list.sw = lapply(scf_list.sw, linkage.stats)
scf_list.sw = lapply(scf_list.sw, sweeps.stats)

scf_list.sw$scf_1@region.names

#mean position extract
scf.pos <- sapply(scf_list.sw$scf_1@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    val <- mean(as.numeric(split))
    return(val)
})

#scf list of positions
genome_pos.list = list()

for(i in 1:length(scf_dirs.rm)){
    genome_pos.list[[scf_dirs.rm[i]]] = sapply(
        scf_list.sw[[scf_dirs.rm[i]]]@region.names,
        function(x){
            split <- strsplit(x," ")[[1]][c(1,3)]
            val <- mean(as.numeric(split))
            return(val)
        }
    )
}

scf_list.sw$scf_1@nucleotide.F_ST
scf_list.sw$scf_1@Tajima.D
scf_list.sw$scf_1@theta_Watterson
scf_list.sw$scf_1@Pi

#Some dplyr BS to format tables

#loop list
######
#Theta
theta_sw_list = list()

for(i in 1:length(scf_dirs.rm)){
    
    theta_sw_list[[ scf_dirs.rm[i] ]] = data.frame(
        scf_list.sw[[i]]@theta_Watterson,
        stringsAsFactors = F
        )
    theta_sw_list[[ scf_dirs.rm[i] ]]$pos = unname(genome_pos.list[[i]])
    theta_sw_list[[ scf_dirs.rm[i] ]]$scf = as.character(scf_names.rm[i,1])
    #print(as.character(scf_names.rm[i,1]))
}

theta_sw.df = do.call(rbind.data.frame, theta_sw_list)
theta_sw.df.ong = pivot_longer(cols = c(-pos, -scf), theta_sw.df, names_to = "pop", values_to = "theta")

#######
#Pi
pi_sw_list = list()

for(i in 1:length(scf_dirs.rm)){
    
    pi_sw_list[[ scf_dirs.rm[i] ]] = data.frame(
    scf_list.sw[[i]]@nuc.diversity.within/10000,
    stringsAsFactors = F
    )
    pi_sw_list[[ scf_dirs.rm[i] ]]$pos = unname(genome_pos.list[[i]])
    pi_sw_list[[ scf_dirs.rm[i] ]]$scf = as.character(scf_names.rm[i,1])
    #print(as.character(scf_names.rm[i,1]))
}

pi_sw.df = do.call(rbind.data.frame, pi_sw_list)
pi_sw.df.long = pivot_longer(cols = c(-pos, -scf), pi_sw.df, names_to = "pop", values_to = "pi")

#######
#Tajima.D
tajD_sw_list = list()

for(i in 1:length(scf_dirs.rm)){
    
    tajD_sw_list[[ scf_dirs.rm[i] ]] = data.frame(
    scf_list.sw[[i]]@Tajima.D,
    stringsAsFactors = F
    )
    tajD_sw_list[[ scf_dirs.rm[i] ]]$pos = unname(genome_pos.list[[i]])
    tajD_sw_list[[ scf_dirs.rm[i] ]]$scf = as.character(scf_names.rm[i,1])
    #print(as.character(scf_names.rm[i,1]))
}

tajD_sw.df = do.call(rbind.data.frame, tajD_sw_list)
tajD_sw.df.long = pivot_longer(cols = c(-pos, -scf), tajD_sw.df, names_to = "pop", values_to = "tajD")

####
#Combine allto list

stats_sw.list = left_join(theta_sw.df, pi_sw.df) %>% left_join(., tajD_sw.df)
stats_sw.list.long = pivot_longer(stats_sw.list, cols = c(pi, theta, tajD), names_to = "stat", values_to = "value")
