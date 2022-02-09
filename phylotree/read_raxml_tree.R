require(treeio)
outg <- c("Stano2", "Ampqui1", "PhaPMI808", "Talpro1", "Talst1_2")

raxml_phylo_fungi <- read.newick("phylotree/RAxML_busco_65.newick")
raxml_phylo_fungi <- phytools::reroot(raxml_phylo_fungi,70, 0.267216461/2)

raxml_phylo_fungi_bt <- read.newick("phylotree/RAxML_busco_65_bootstrap.newick")
raxml_phylo_fungi_bt_values <- prop.clades(raxml_phylo_fungi, raxml_phylo_fungi_bt)

raxml_phylo_fungi$node.label <- raxml_phylo_fungi_bt_values
