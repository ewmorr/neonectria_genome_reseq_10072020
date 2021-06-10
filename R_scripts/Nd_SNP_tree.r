#require(ggtree)
require(ape)
require(ggdendro)
require(dendextend)
require(phytools)
require(RColorBrewer)
require(dplyr)
#require(tidytree)
require(ggtree)
require(ggplot2)
get.color.palette <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)

metadata = read.table("~/GARNAS_neonectria_genome_reseq_10072020/sample_metadata/sample_metadata.txt", header = T)
colnames(metadata) = c("Site", "Tree", "Plug", "Peri", "label", "species", "State")


states = metadata$State %>% unique
state.colors = get.color.palette("Paired")

state.colors = data.frame(State = states, color = state.colors[1:length(states)], stringsAsFactors = F)
metadata.state_colors = full_join(metadata, state.colors)
metadata.state_colors.Nd = metadata.state_colors %>% filter(species == "Nd")
metadata.state_colors.Nd$State = as.character(metadata.state_colors.Nd$State)
metadata.state_colors.Nd$color = as.character(metadata.state_colors.Nd$color)

metadata.state_colors.Nd[metadata.state_colors.Nd$label == "REF","State"] = as.character("REF")
metadata.state_colors.Nd[metadata.state_colors.Nd$abel == "REF","color"] = "black"

Nd.ML_tree = read.tree("~/GARNAS_neonectria_genome_reseq_10072020/SPANDx_Nd/ML_phylogeny.tre")
Nd.MP_tree = read.tree("~/GARNAS_neonectria_genome_reseq_10072020/SPANDx_Nd/MP_phylogeny.tre")

#tidytree way . The plotting step sucks up A TON of CPU and doesn't produce output after several minutes
#x = as_tibble(Nd.ML_tree)

#x.metadata = full_join(x, metadata, by = "label")
#x.treedat = as.treedata(x.metadata)
#ggtree(x.treedat)


#ggtree way
#Tip nodes are arranges as 1:N (i.e., the first nodes in the list) so that calling the number of samples to subset will plot only tip nodes

ggtree(Nd.ML_tree, branch.length='branch length', layout="circular") +
geom_point2(size = 3,
    aes(
        subset=(
            node %in% seq(1:length(Nd.ML_tree$tip.label))
        ),
        color = c(metadata.state_colors.Nd$State, rep("white", length(Nd.ML_tree$node.label)))
    )
) +
scale_color_brewer(palette = "Paired") +
labs(color = "State/region")


ggtree(Nd.ML_tree, branch.length='none', layout="circular") +
geom_point2(size = 3,
aes(
subset=(
node %in% seq(1:length(Nd.ML_tree$tip.label))
),
color = c(metadata.state_colors.Nd$State, rep("white", length(Nd.ML_tree$node.label)))
)
) +
scale_color_brewer(palette = "Paired") +
labs(color = "State/region")


ggtree(Nd.ML_tree, branch.length='branch length', layout="rectangular") +
geom_point2(size = 3,
aes(
subset=(
node %in% seq(1:length(Nd.ML_tree$tip.label))
),
color = c(metadata.state_colors.Nd$State, rep("white", length(Nd.ML_tree$node.label)))
)
) +
scale_color_brewer(palette = "Paired") +
labs(color = "State/region")

ggtree(Nd.MP_tree, branch.length='branch length', layout="rectangular") +
geom_point2(size = 3,
aes(
subset=(
node %in% seq(1:length(Nd.ML_tree$tip.label))
),
color = c(metadata.state_colors.Nd$State, rep("white", length(Nd.ML_tree$node.label)))
)
) +
scale_color_brewer(palette = "Paired") +
labs(color = "State/region")

ggtree(Nd.ML_tree, branch.length='none', layout="rectangular") +
geom_point2(size = 3,
aes(
subset=(
node %in% seq(1:length(Nd.ML_tree$tip.label))
),
color = c(metadata.state_colors.Nd$State, rep("white", length(Nd.ML_tree$node.label)))
)
) +
scale_color_brewer(palette = "Paired") +
labs(color = "State/region")

#convert to derogram using ape and phytools
Nd.ML_tree.dendro = as.dendrogram(force.ultrametric(ape::root(Nd.ML_tree, outgroup = "REF")))
#ggdendrogram(Nd.ML_tree.dendro, hang = T)
ggdendrogram(Nd.ML_tree.dendro)
ddata <- dendro_data(Nd.ML_tree.dendro, type = "rectangle")

ddata$labels = full_join(ddata$labels, metadata.state_colors.Nd)

ggplot() +
geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
geom_point(data = ddata$labels, aes(x = x, y = y, color = State)) +
scale_color_brewer(palette = "Paired")




Nd_tree.labels = Nd.ML_tree.dendro %>% labels

Nd_tree.metadata = full_join(
    data.frame(label = Nd_tree.labels),
    metadata %>% select(label, State)
) %>% full_join(
    .,
    state.colors
)
Nd_tree.metadata$State = Nd_tree.metadata$State %>% as.character
Nd_tree.metadata[Nd_tree.metadata$Sequence_label == "REF","State"] = "REF"
Nd_tree.metadata[Nd_tree.metadata$Sequence_label == "REF","color"] = "black"


Nd.ML_tree.dendro %>%
    set("leaves_col", Nd_tree.metadata.Nd$color) %>%
    set("leaves_pch", 19) %>%
    set("leaves_cex", 2)
as.ggdend(Nd.ML_tree.dendro)
Nd.ML_tree.dendro %>% get_nodes_attr("node_par")

plot(Nd.ML_tree.dendro)

ggplot() +
geom_segments(ddata$segmetns)







