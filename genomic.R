source("shared_setup.R")
## genome stats
genome_size_tb <- read_tsv("meta_tables/genome_size.tsv")
gene_count_tb <- read_tsv("meta_tables/genome_gene_count.tsv")
gene_length_tb <- read_tsv("meta_tables/genome_gene_length.tsv")
SMC_tb <- read_csv("SMC/JGI_SMC.csv") %>%
  rename(SMC=Total, spcode=portal_id)

genome_stats <-
  full_join(genome_size_tb, gene_count_tb) %>%
  full_join(gene_length_tb) %>%
  mutate(CDS_length=gene_length*gene_count)

## completeness
completeness_tb <- read_tsv("meta_tables/genome_completeness.tsv")

## cazyme table
cazyme_rawtb <- read_csv("cazyme/Cazyme_table.csv") #full cazyme table
cazyme_rawtb <- cazyme_rawtb %>%
  select(-`Annotation Description`, -`Total`) #filter non-species columns
cazyme_rawtb <- select_if(cazyme_rawtb, ~!all(is.na(.))) #filter non-information species (all NA)
cazyme_rawtb <- cazyme_rawtb %>%
  rename(cazyme=1) #rename cazyme column
cazyme_rawtb[is.na(cazyme_rawtb)] <- 0 #fill NA with 0

cazyme_tb <- transpose_tibble(cazyme_rawtb, "cazyme", "spcode") #transpose table -> cazyme as columns

## Secreted Protein
secretome_tb <-
  read_tsv("secretome/secretome_count.tsv")

## Repeat elements
repeat_tb <-
  read_tsv("results/repeat_length.tsv", col_names = c("spcode", "TE_len"))

## Resistomes
rgi_tb <- 
  read_tsv("rgi/rgi_tb.tsv")

## protease
protease_tb <-
  read_tsv("results/merops.tsv", col_names = FALSE) %>%
  rename(spcode=X1, protease=X2)

## lipase
lipase_tb <-
  read_tsv("results/led.tsv", col_names = FALSE) %>%
  rename(spcode=X1, lipase=X2)

## all stats
genome_allstat_tb <-
  select(all_genome_tb, spcode) %>%
  left_join(select(taxa_tb, spcode, species, phylum, class),
            by="spcode") %>%
  left_join(select(cazyme_tb, spcode, CAZy)) %>%
  left_join(completeness_tb) %>%
  left_join(genome_stats) %>%
  left_join(select(SMC_tb, spcode, SMC)) %>%
  left_join(secretome_tb) %>%
  left_join(select(rgi_tb, spcode, resistome=total)) %>%
  left_join(repeat_tb) %>%
  mutate(`duplicated BUSCOs (%)`= 100*`Complete and duplicated BUSCOs (D)`/`Total BUSCO groups searched`) %>%
  mutate(genome_size=genome_size/1e6) %>%
  mutate(TE_len=TE_len/1e6) %>%
  left_join(lipase_tb) %>%
  left_join(protease_tb)

write_csv(genome_allstat_tb, "results/genome_stats.csv")

# Phylogeny
all_phylo_dt <-
  eco_tb %>%
  left_join(taxa_tb) %>%
  left_join(select(all_genome_tb, spcode, Name)) %>%
  sub_treeXdata("spcode", raxml_phylo_fungi)

all_phylo_input <-
  all_phylo_dt$tree %>%
  as_tibble() %>%
  left_join(all_phylo_dt$data, by=c("label"="spcode")) %>%
  mutate(anno=NA) %>%
  mutate(anno=replace(anno, node==68, "Mucoromycotina")) %>%
  mutate(anno=replace(anno, node==70, "Mortierellomycotina")) %>%
  mutate(anno=replace(anno, node==72, "Basidiomycota")) %>%
  mutate(anno=replace(anno, node==80, "Ascomycota")) %>%
  mutate(ecology=factor(ecology, levels = ecology_order)) %>%
  as.treedata()

all_phylo_spcode_plot <-
  ggtree(all_phylo_input) +
  geom_tiplab(aes(label=label, color=ecology), size=2.5) +
  geom_nodelab(aes(label=anno), vjust=1, size=3) +
  xlim(0, 2.5) +
  ylim(-0.5, 65)

ggsave("plots/phylogeny.pdf",
       all_phylo_spcode_plot,
       width = 10, height = 6)

# genomic stats contrast and phylogenetic comparison
input_vars <- c("CAZy", "genome_size", "gene_count", "SMC",
                "n_secreted", "TE_len", "resistome",
                "protease", "lipase")

constrast_tree_dt<-
  contrast_tb %>%
  left_join(genome_allstat_tb) %>%
  filter(spcode %in% pull(contrast_tb, spcode)) %>%
  sub_treeXdata("spcode", raxml_phylo_fungi)

constrast_tree_dt$tree$node.label <- NULL

constrast_tree_input <-
  constrast_tree_dt$tree %>%
  as_tibble() %>%
  left_join(constrast_tree_dt$data, by=c("label"="spcode")) %>%
  mutate(label=replace(label, node==39, "Mucoromycota")) %>%
  mutate(label=replace(label, node==41, "Mortierellomycota")) %>%
  mutate(label=replace(label, node==43, "Basidiomycota")) %>%
  mutate(label=replace(label, node==48, "Ascomycota")) %>%
  as.treedata()

# Wilcoxon test
contrast_w_tb <-
  map_dfr(input_vars, ~wilcox_wrapper(.data = constrast_tree_dt$data, pre_col = "ecology", .x)) %>%
  arrange(p)

write_csv(contrast_w_tb, "results/genome_Wilcoxon.csv")

contrast_phylolm_fits <-
  map(input_vars, ~phylolm_wrapper(.data = constrast_tree_dt$data, id_col = "spcode", pre_cols = "ecology", res_cols = .x, phy = constrast_tree_dt$tree, model="OUfixedRoot"))

contrast_phylolm_tb <-
  map_dfr(contrast_phylolm_fits, extract_phylolm_fit) %>%
  filter(term!="(Intercept)") %>%
  select(-term)

write_csv(contrast_phylolm_tb, "results/genome_PLR.csv")

## plot
sig1 <- geom_signif(comparisons=list(c("saprobe", "poplar endophyte")), annotations="*", textsize=6, vjust = 0.5)
sig2 <- geom_signif(comparisons=list(c("saprobe", "poplar endophyte")), annotations="**", textsize=6, vjust = 0.5)
sig3 <- geom_signif(comparisons=list(c("saprobe", "poplar endophyte")), annotations="***", textsize=6, vjust = 0.5)

box_theme <-
  theme(legend.position="none",
        plot.margin = unit(c(0.2,0,0.2,0.1), "inch"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, angle=90),
        axis.title.y = element_blank()
  )

univar_plotlist <- list()

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "genome_size")
univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      sig2 +
      ylab("genome size (Mbp)") +
      xlab(NULL) +
      box_theme)
  )

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "gene_count")
univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      sig2 +
      ylab("gene count") +
      xlab(NULL) +
      box_theme)
  )

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "TE_len")
univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      ylab("repeated element length (Mbp)") +
      xlab(NULL) +
      box_theme)
  )

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "CAZy")

univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      sig2 +
      ylab("cazyme count") +
      xlab(NULL) +
      box_theme)
  )

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "n_secreted")
univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      sig2 +
      ylab("small secreted proteins") +
      xlab(NULL) +
      box_theme)
  )

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "protease")
univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      sig2 +
      ylab("protease") +
      xlab(NULL) +
      box_theme)
  )

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "lipase")
univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      sig2 +
      ylab("lipase") +
      xlab(NULL) +
      box_theme)
  )

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "resistome")
univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      ylab("antibiotic resistance genes") +
      xlab(NULL) +
      box_theme)
  )

p <- univariate_plots(constrast_tree_dt$data, id_col = "spcode", group_col = "ecology", var_col = "SMC")

univar_plotlist <-
  c(univar_plotlist, list(
    p[[2]] +
      ylab("secondary metabolite gene clusters") +
      xlab(NULL) +
      box_theme)
  )

univar_plotlist <- c(univar_plotlist)

constrast_boxplots <- plot_grid(plotlist=univar_plotlist, greedy=FALSE, nrow=1)

ggsave("plots/genome_boxplot.pdf",
       constrast_boxplots,
       width=17,
       height=4,
       units="cm")

fea_var_order <-
  c("genome\nsize\n(Mbp)**",
    "gene\n  count(k)**",
    "repeated\nelement length(Mbp)",
    "cazyme\n  count**",
    "small secreted\n  proteins**",
    "protease**",
    "lipase**",
    "antibiotic\nresistance\ngenes",
    "secondary\nmetabolite\nclusters"
  )

constrast_tree <- ggtree(constrast_tree_input, branch.length="none") +
  geom_tiplab(aes(label=species, color=ecology), size=2) +
  xlim(0, 20) +
  theme(legend.position = "none")

colh_plots <-
  constrast_tree$data %>%
  select(spcode=label, ecology, all_of(input_vars)) %>%
  pivot_longer(input_vars) %>%
  mutate(
    name=replace(name, name=="genome_size", "genome\nsize\n(Mbp)**"),
    name=replace(name, name=="gene_count", "gene\n  count(k)**"),
    name=replace(name, name=="n_secreted", "small secreted\n  proteins**"),
    name=replace(name, name=="SMC", "secondary\nmetabolite\nclusters"),
    name=replace(name, name=="CAZy", "cazyme\n  count**"),
    name=replace(name, name=="protease", "protease**"),
    name=replace(name, name=="lipase", "lipase**"),
    name=replace(name, name=="TE_len", "repeated\nelement length(Mbp)"),
    name=replace(name, name=="resistome", "antibiotic\nresistance\ngenes")) %>%
  mutate(name=factor(name,levels=fea_var_order)) %>%
  ggplot(aes(x=value, y=spcode, fill=ecology)) +
  geom_colh() +
  facet_wrap("name", scales="free", nrow=1) +
  theme_classic() +
  ylab(NULL) +
  scale_y_discrete(breaks=NULL) +
  xlab(NULL) +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(size=4)
  )

treecol_plot <-
  insert_left(colh_plots, constrast_tree, 2/8)

ggsave("plots/genome_treecol.pdf",
       treecol_plot,
       width=16,
       height=14,
       units="cm")
