require(PopGenome)
require(dplyr)
require(tidyr)
require(ggplot2)
source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

#Can make this a list and lapply

chr.1 = readData(
"Nf_post_SPANDx/scf_1",
include.unknown = T,
format = "VCF"
)

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


#Calculate stats on individual scf to see if these (avgs across scf) reflect the previous calcs from concatenated version

scf_list.pops = lapply(scf_list.pops, F_ST.stats)
scf_list.pops = lapply(scf_list.pops, neutrality.stats)
scf_list.pops = lapply(scf_list.pops, linkage.stats)
scf_list.pops = lapply(scf_list.pops, diversity.stats)

#Calculate mean stats per pop
scf_stats_list = list()

#Tajima's D
for(i in 1:length(scf_dirs)){
    #print(scf_dirs[i])
    #print(scf_list.pops[[scf_dirs[i]]]@Tajima.D)
    scf_stats_list[[scf_dirs[i]]] = scf_list.pops[[scf_dirs[i]]]@Tajima.D
}
scf_stats.tajD = do.call(rbind.data.frame, scf_stats_list)
pop_tajD = apply(scf_stats.tajD, 2, mean, na.rm = T)

#Pi
for(i in 1:length(scf_dirs)){
    #print(scf_dirs[i])
    #print(scf_list.pops[[scf_dirs[i]]]@Tajima.D)
    scf_stats_list[[scf_dirs[i]]] = scf_list.pops[[scf_dirs[i]]]@nuc.diversity.within
}
scf_stats.Pi = do.call(rbind.data.frame, scf_stats_list)
pop_Pi = apply(scf_stats.Pi, 2, mean, na.rm = T)

#Looks like Pi is not normalized by se len so get n.sites
for(i in 1:length(scf_dirs)){
    #print(scf_dirs[i])
    #print(scf_list.pops[[scf_dirs[i]]]@Tajima.D)
    scf_stats_list[[scf_dirs[i]]] = scf_list.pops[[scf_dirs[i]]]@n.sites
}

scf_names = do.call(rbind.data.frame, lapply(scf_stats_list, names))
colnames(scf_names) = "scf_name"
scf_names$scf = scf_dirs

scf_stats.nSites = do.call(rbind.data.frame, scf_stats_list)
#n.sites by scf

#norm Pi by seq len
scf_stats.Pi_norm = scf_stats.Pi/scf_stats.nSites[,1]
pop_Pi.norm = apply(scf_stats.Pi_norm, 2, mean, na.rm = T)

#Try theta (it's not listed in slots but found it below)
for(i in 1:length(scf_dirs)){
    #print(scf_dirs[i])
    #print(scf_list.pops[[scf_dirs[i]]]@Tajima.D)
    scf_stats_list[[scf_dirs[i]]] = scf_list.pops[[scf_dirs[i]]]@theta_Watterson
}
scf_stats.theta = do.call(rbind.data.frame, scf_stats_list)
pop_theta = apply(scf_stats.theta, 2, mean, na.rm = T)

#norm theta by seq len
scf_stats.theta_norm = scf_stats.theta/scf_stats.nSites[,1]
pop_theta.norm = apply(scf_stats.theta_norm, 2, mean, na.rm = T)


#combining with site names and plotting
site_div = left_join(
data.frame(
state = pop.levels,
pop.name = paste("pop", seq(1:7)),
Pi = pop_Pi,
Pi.norm = pop_Pi.norm,
Tajima.D = pop_tajD,
theta = pop_theta
),
site.info.all[site.info.all$Site != "ADN2" & site.info.all$Site != "ADS2", ],
by = "state"
)

plot(Pi.norm ~ duration_infection, data = site_div)
summary(lm(Pi.norm ~ duration_infection, data = site_div))
plot(Tajima.D ~ duration_infection, data = site_div)
summary(lm(Tajima.D ~ duration_infection, data = site_div)) #R2 = 0.64, P = 0.03
plot(theta ~ duration_infection, data = site_div)
summary(lm(theta ~ duration_infection, data = site_div))

p1 = ggplot(site_div, aes(x = duration_infection, y = Pi.norm)) +
geom_point() +
#geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "BBD infection duration (yrs)", y = expression(pi)) +
my_gg_theme

p2 = ggplot(site_div, aes(x = duration_infection, y = Tajima.D)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "BBD infection duration (yrs)", y = "Tajima's D") +
my_gg_theme

pdf("figures/Pi_v_dur_inf.pdf", width = 6, height = 4)
p1
dev.off()

pdf("figures/TajD_v_dur_inf.pdf", width = 6, height = 4)
p2
dev.off()


######################
#Trying sliding window

scf_list.sw = lapply(scf_list.pops, sliding.window.transform, 10000, 10000, type = 2, whole.data = T)
#This throws an error about attepmting to select one element. Probably the scf with only SNP at pos 7

sw_try = sliding.window.transform(scf_list.pops[["scf_14"]], 10000, 10000, type = 2, whole.data = T)
sw_try = sliding.window.transform(scf_list.pops[["scf_11"]], 10000, 10000, type = 2, whole.data = T)
#That seems to be the problem

#Rm element 14
scf_list.pops.rm = scf_list.pops
scf_list.pops.rm[["scf_14"]] = NULL
scf_dirs.rm = scf_dirs[! scf_dirs == "scf_14"]
scf_names.rm = scf_names[-14,]

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

#Some dplyr BS to format tables

#First trying on one scf
scf_1.pw_nuc_FST = data.frame(t(scf_list.sw[[1]]@nuc.F_ST.pairwise), stringsAsFactors = F)
scf_1.pw_nuc_FST$pos = unname(genome_pos.list[[1]])
scf_1.pw_nuc_FST$scf = as.character(scf_names.rm[1,1])

scf_2.pw_nuc_FST = data.frame(t(scf_list.sw[[2]]@nuc.F_ST.pairwise), stringsAsFactors = F)
scf_2.pw_nuc_FST$pos = unname(genome_pos.list[[2]])
scf_2.pw_nuc_FST$scf = as.character(scf_names.rm[2,1])

#loop list
FST_sw_list = list()

for(i in 1:length(scf_dirs.rm)){
    
    FST_sw_list[[ scf_dirs.rm[i] ]] = data.frame(
        t(scf_list.sw[[i]]@nuc.F_ST.pairwise),
        stringsAsFactors = F
        )
    FST_sw_list[[ scf_dirs.rm[i] ]]$pos = unname(genome_pos.list[[i]])
    FST_sw_list[[ scf_dirs.rm[i] ]]$scf = as.character(scf_names.rm[i,1])
    #print(as.character(scf_names.rm[i,1]))
}

FST_sw.df = do.call(rbind.data.frame, FST_sw_list)
FST_sw.df.long = pivot_longer(cols = c(-pos, -scf), FST_sw.df, names_to = c("pop1", "pop2"), names_sep = "\\.", values_to = "Fst")


#List to add state names instead of pop#

pops_and_states1 = data.frame(
state1 = pop.levels,
pop1 = paste("pop", seq(1:7), sep = ""),
stringsAsFactors = F
)
pops_and_states2 = data.frame(
state2 = pop.levels,
pop2 = paste("pop", seq(1:7), sep = ""),
stringsAsFactors = F
)

FST_sw.df.long = left_join(FST_sw.df.long, pops_and_states1) %>% left_join(., pops_and_states2)
FST_sw.df.long$state1 %>% unique
FST_sw.df.long$scf %>% unique

#Add lengths to filter scf by len
scf_lens = read.table("Nf_post_SPANDx/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scf", "length")
FST_sw.df.long = left_join(FST_sw.df.long, scf_lens)
FST_sw.df.long$scf %>% unique

FST_sw.df.long %>% nrow
FST_sw.df.long %>% filter(state1 == "ME.S" | state2 == "ME.S") %>% nrow
( FST_sw.df.long %>% filter(state1 == "ME.S" | state2 == "ME.S") )$state1 %>% unique
( FST_sw.df.long %>% filter(state1 == "ME.S" | state2 == "ME.S") )$state2 %>% unique

#set negative Fst to zero
FST_sw.df.long$Fst[FST_sw.df.long$Fst < 0] = 0

#Add column of state comparison
FST_sw.df.long$comp = paste(FST_sw.df.long$state1, FST_sw.df.long$state2, sep = "-")

site_div %>% select(state, duration_infection, HDD4.mean_nongrowing, MAT)
#Identify outliers on data filtered by state comps

#ME.S COMPS
FST_sw.df.long.MES = FST_sw.df.long %>% filter(state1 == "ME.S" | state2 == "ME.S") %>% filter(length > 100000)
FST_sw.df.long.MES$scf %>% unique

#Outlier tests
# identify the 95% percentile
my_threshold <- quantile(FST_sw.df.long.MES$Fst, 0.95, na.rm = T)
# make an outlier column in the data.frame
FST_sw.df.long.MES <- FST_sw.df.long.MES %>% mutate(outlier = ifelse(Fst > my_threshold, "outlier", "background"))
#Number of outliers
FST_sw.df.long.MES %>% group_by(outlier) %>% tally()

#arrange comps by coldest to warmest
site_div %>% select(state, duration_infection, HDD4.mean_nongrowing, MAT)
FST_sw.df.long.MES$comp %>% unique
FST_sw.df.long.MES$comp = factor(FST_sw.df.long.MES$comp, levels = c("ME.S-NY.N", "ME.S-NY.S", "ME.S-PA", "ME.S-NH", "ME.S-WV", "ME.S-NC"))

#Read LFMM results
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

Fst_min_max = c(0,1)


#Plot with ggplot
p1 = ggplot(FST_sw.df.long.MES, aes(x = pos/10^6, y = Fst, color = outlier)) +
facet_grid(comp ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(breaks = c(0,0.5,1)) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("F"[ST]))) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

p2 = ggplot(FST_sw.df.long.MES, aes(x = pos/10^6, y = Fst, alpha = outlier)) +
facet_grid(comp ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
geom_point(size = 1, color = "black") +
geom_line(color = "black", alpha = 0.5) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(breaks = c(0,0.5,1)) +
scale_alpha_manual(values = c(0, 0.5), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("F"[ST]))) +
theme(
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

p3 = ggplot() +
geom_linerange(
    data = pv.hdd.with_pos %>% filter(FDR.sig == "sig" & length > 100000),
    aes(x = pos/10^6, ymin = pi_min_max[1], ymax = Fst_min_max[2]),
    color = "blue",
    alpha = 0.5
) +
geom_point(data = FST_sw.df.long.MES, aes(x = pos/10^6, y = Fst, alpha = outlier), color = "black", size = 1) +
geom_line(data = FST_sw.df.long.MES, aes(x = pos/10^6, y = Fst), color = "black", alpha = 0.5) +
facet_grid(comp ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(breaks = c(0,0.5,1)) +
scale_alpha_manual(values = c(0, 0.5), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("F"[ST]))) +
theme(
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

pdf("figures/Fst.sw_scan.MES_comps.LFMM.pdf", width = 18, height = 6)
p2
p3
dev.off()


#NH COMPS
FST_sw.df.long.NH = FST_sw.df.long %>% filter(state1 == "NH" | state2 == "NH") %>% filter(length > 100000)
FST_sw.df.long.NH$scf %>% unique

#Outlier tests
# identify the 95% percentile
my_threshold <- quantile(FST_sw.df.long.NH$Fst, 0.95, na.rm = T)
# make an outlier column in the data.frame
FST_sw.df.long.NH <- FST_sw.df.long.NH %>% mutate(outlier = ifelse(Fst < my_threshold, "outlier", "background"))
#Number of outliers
FST_sw.df.long.NH %>% group_by(outlier) %>% tally()


#Plot with ggplot
ggplot(FST_sw.df.long.NH, aes(x = pos/10^6, y = Fst, color = outlier)) +
#facet_wrap(~scaffold) +
facet_grid(comp ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
#scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("black", "grey"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("F"[ST]))) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)


#NC COMPS
FST_sw.df.long.NC = FST_sw.df.long %>% filter(state1 == "NC" | state2 == "NC") %>% filter(length > 100000)
FST_sw.df.long.NC$scf %>% unique

#Outlier tests
# identify the 95% percentile
my_threshold <- quantile(FST_sw.df.long.NC$Fst, 0.95, na.rm = T)
# make an outlier column in the data.frame
FST_sw.df.long.NC <- FST_sw.df.long.NC %>% mutate(outlier = ifelse(Fst < my_threshold, "outlier", "background"))
#Number of outliers
FST_sw.df.long.NC %>% group_by(outlier) %>% tally()


#Plot with ggplot
ggplot(FST_sw.df.long.NC, aes(x = pos/10^6, y = Fst, color = outlier)) +
#facet_wrap(~scaffold) +
facet_grid(comp ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
#scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("black", "grey"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("F"[ST]))) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)


#NYN COMPS
FST_sw.df.long.NY.N = FST_sw.df.long %>% filter(state1 == "NY.N" | state2 == "NY.N") %>% filter(length > 100000)
FST_sw.df.long.NY.N$scf %>% unique

#Outlier tests
# identify the 95% percentile
my_threshold <- quantile(FST_sw.df.long.NY.N$Fst, 0.95, na.rm = T)
# make an outlier column in the data.frame
FST_sw.df.long.NY.N <- FST_sw.df.long.NY.N %>% mutate(outlier = ifelse(Fst < my_threshold, "outlier", "background"))
#Number of outliers
FST_sw.df.long.NY.N %>% group_by(outlier) %>% tally()


#Plot with ggplot
ggplot(FST_sw.df.long.NY.N, aes(x = pos/10^6, y = Fst, color = outlier)) +
#facet_wrap(~scaffold) +
facet_grid(comp ~ scf, scales = "free_x", space='free_x') + #grid by state2 bc state1 is all ME.S
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = faNY.Ny_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
#scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("black", "grey"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = expression(paste("F"[ST]))) +
theme(
strip.text.x = element_blank(),
strip.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

