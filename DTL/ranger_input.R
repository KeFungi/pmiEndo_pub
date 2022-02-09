library(treeio)
library(ape)
library(tidyverse)
source("phylotree/read_raxml_tree.R")
GH88_tree <- read.newick("DTL/GH88.newick")
GH88_tree$tip.label <-
  GH88_tree$tip.label %>%
  str_replace("^[^|]+\\|([^|]+\\|[^|]+).+", "\\1") %>%
  str_replace("_[0-9]", "") %>%
  str_replace("\\|", "_")

raxml_phylo_fungi$tip.label <-
  raxml_phylo_fungi$tip.label %>%
  str_replace("_[0-9]", "") %>%
  str_replace("\\|", "_")

sp_list <-
  GH88_tree$tip.label %>%
  str_replace("_.+", "") %>%
  unique()

raxml_phylo_fungi_GH88 <-
  raxml_phylo_fungi %>%
  keep.tip(sp_list)

raxml_phylo_fungi$edge.length <- NULL
raxml_phylo_fungi$node.label <- NULL
GH88_tree$edge.length <- NULL
GH88_tree$node.label <- NULL
write.tree(raxml_phylo_fungi, "DTL/sp_tree.newick")
write.tree(GH88_tree, "DTL/gene_tree.newick")

l1 <- read_lines("DTL/sp_tree.newick")
l2 <- read_lines("DTL/gene_tree.newick")

write_lines(c(l1, l2), "DTL/rager_input.newick")
write_lines(c(l1, l2), "plots/Supp17_rager_input.newick")
