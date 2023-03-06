# Tree for Timberlake MAGs

# Libs
library(ape)
library(picante)
library(dendextend)

# Tree from KBase made with 49 COGs
TL_tree <- read.tree("~/Desktop/Timberlake/dRep_Tree_20-labels.newick")

# Prune the one non-MAG
TL_tree <- drop.tip(TL_tree, tip = "Frateuriasp.Soil773")

# Trim the tip labels
TL_tree$tip.label <- substr(TL_tree$tip.label, start = 1, stop = nchar(TL_tree$tip.label)-11)

# Trim the node labels
TL_tree$node.label <- as.character(round(as.numeric(TL_tree$node.label), digits = 2))
TL_tree$node.label <- replace_na(TL_tree$node.label, replace = "")

# Plot with ape
png("~/Desktop/Timberlake/Tree.png", width = 6.5, height = 6, units = "in", res = 300)
par(oma = c(1,0,0,1))
plot.phylo(TL_tree,
           align.tip.label = T,
           no.margin = T,
           font = 1,
           cex = 0.6,
           edge.width = 2,
           show.node.label = T,
           node.pos = 1,
           label.offset = 0.005,
           adj = 0)
add.scale.bar(x = 0.5, y = 0.5)
title("MAG Phylogeny\nfrom 49 COGs", adj = 0.1, line = -3)
text(x = 0.15, y = 1.5, label = "Archaea")
text(x = 0.15, y = 6.6, label = "Bacteria")
dev.off()
