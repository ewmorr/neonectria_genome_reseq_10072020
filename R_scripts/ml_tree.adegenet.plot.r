require(dplyr)
require(ggplot2)
require(phylogram)
require(ggdendro)
require(dendextend)
require(zoo)

#############
#Read tree
source("R_scripts/ggplot_theme.txt")

fit = readRDS("data/intermediate_RDS/Nf_ML_tree.rds")
fit$tree
Nf.dendro = as.dendrogram.phylo(fit$tree)

#Read metadata
#metadata
source("R_scripts/make_site_metadata.r")
sample_metadata
site_coords = read.table("data/sample_metadata/site_coords_for_map.txt", header = T)
ind.metrics = 
    left_join(
      data.frame(sample = fam_info[,1]), 
      sample_metadata %>% select(sample,state.name) 
    ) %>% 
    left_join(., site_coords %>% select(state.name, lat, lon))
#########

colnames(sample_metadata) = c("label", "state.name")

###################################
#pretty plot of tree with ggdendro and using dendextend functions to extract leaves

leaves = Nf.dendro %>% dendextend::get_nodes_attr("leaf")
leaf_height = Nf.dendro %>% dendextend::get_leaves_attr("height")


#######################################
ddata <- ggdendro::dendro_data(Nf.dendro, type = "rectangle")

ddata$labels$leaf_height = leaf_height
ddata$labels = full_join(ddata$labels, sample_metadata)

p1 = ggplot() +
  geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = ddata$labels, aes(x = x, y = leaf_height, color = state.name), size = 3) +
  scale_color_manual(values = c25) +
  coord_flip() +
 # coord_polar() +
  scale_y_reverse() +
  labs(color = "State/region") +
  theme_dendro() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
p1

pdf("figures/Nf_ML_tree.pdf", height = 12, width = 6)
p1
dev.off()


p2 = ggplot() +
  geom_segment(data = ddata$segments, aes(x = x, y = y/10, xend = xend, yend = yend/10)) +
  geom_point(data = ddata$labels, aes(x = x, y = leaf_height/10, color = state.name), size = 3) +
  scale_color_manual(values = c25) +
  coord_flip() +
  coord_polar() +
  scale_y_reverse() +
  labs(color = "State/region") +
  theme_dendro() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

pdf("figures/Nf_ML_tree_polar.pdf", height = 10, width = 10)
p2
dev.off()


################################
#ggdebdro with coloring of branches

cut <- 3   # Number of clusters
dendr <- ggdendro::dendro_data(Nf.dendro, type = "rectangle") 
clust <- cutree(Nf.dendro, k = cut)               # find 'cut' clusters
clust.df <- data.frame(label = names(clust), cluster = clust)

# Split dendrogram into upper grey section and lower coloured section
height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
cut.height <- mean(c(height[cut], height[cut-1]))
dendr$segments$line <- ifelse(dendr$segments$y == dendr$segments$yend &
                                dendr$segments$y > cut.height, 1, 2)
dendr$segments$line <- ifelse(dendr$segments$yend  > cut.height, 1, dendr$segments$line)

# Number the clusters
dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
change <- which(dendr$segments$cluster == 1)
for (i in 1:cut) dendr$segments$cluster[change[i]] = i + 1
dendr$segments$cluster <-  ifelse(dendr$segments$line == 1, 1, 
                                  ifelse(dendr$segments$cluster == 0, NA, dendr$segments$cluster))
dendr$segments$cluster <- zoo::na.locf(dendr$segments$cluster) 


#####################
#Add heights and cols to leaves for plotting points
leaves = Nf.dendro %>% dendextend::get_nodes_attr("leaf")
leaf_height = Nf.dendro %>% dendextend::get_leaves_attr("height")
dendr$labels$leaf_height = leaf_height
dendr$labels = full_join(dendr$labels, sample_metadata)




#with points

p3 = ggplot() + 
  geom_segment(data = segment(dendr), 
               aes(x=x, y=y, xend=xend, yend=yend, size=factor(line), colour=factor(cluster)), 
               lineend = "square", show.legend = FALSE) + 
  scale_colour_manual(values = c("grey60", rainbow(cut))) +
  scale_size_manual(values = c(.1, 1)) +
  geom_point(data = dendr$labels, aes(x = x, y = leaf_height, fill = state.name), size = 3, shape = 21) +
  scale_fill_manual(values = c25) +
  scale_y_reverse(expand = c(0.2, 0)) + 
  labs(x = NULL, y = NULL, fill = "State/region") +
  coord_flip() +
  theme_dendro() +
  theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)
        )
p3

pdf("figures/Nf_ML_tree.k_4.pdf", height = 12, width = 6)
p3
dev.off()




p4 = ggplot() + 
  geom_segment(data = segment(dendr), 
               aes(x=x, y=y, xend=xend, yend=yend, size=factor(line), colour=factor(cluster)), 
               lineend = "square", show.legend = FALSE) + 
  scale_colour_manual(values = c("grey60", rainbow(cut))) +
  scale_size_manual(values = c(.1, 1)) +
  geom_point(data = dendr$labels, aes(x = x, y = leaf_height, fill = state.name), size = 3, shape = 21) +
  scale_fill_manual(values = c25) +
  scale_y_reverse(expand = c(0.2, 0)) + 
  labs(x = NULL, y = NULL, fill = "State/region") +
  coord_flip() +
  coord_polar() +
  theme_dendro() +
  theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)
  )
p4

pdf("figures/Nf_ML_tree_polar.k_4.pdf", height = 10, width = 10)
p4
dev.off()
