# Libraries
library(tidyverse)
library(magrittr)
library(pheatmap)
library(caret)
library(phylolm)
library(treeio)
library(ggtree)
library(ape)
library(aplot)
library(ggVennDiagram)
library(ggsignif)
library(ggstance)
library(ranger)
library(cowplot)
library(exactRankTests)

# Source
source("scripts/analysis_wrapper.R")
source("phylotree/read_raxml_tree.R")
source("phylotree/tree_funs.R")
source("phylotree/phyloglmstep_fixed.R")
source("phylotree/comparative_fixed.data.R")
set.seed(1)

# Read data
## taxa
taxa_tb <- read_csv("meta_tables/genome_relatives.csv") %>%
  select(
    spcode=portal_ID,
    strain,
    phylum,
    class,
    order,
    family,
    genus,
    species
  ) %>%
  column_to_rownames("spcode") %>%
  as_tibble(rownames = "spcode")

## genome metadata and ecology
all_genome_tb <-
  read_csv("meta_tables/genome_relatives.csv") %>%
  rename(spcode=portal_ID)

eco_tb <-
  all_genome_tb %>%
  select(spcode, ecology) %>%
  mutate(ecology=replace(ecology, ecology=="poplar endo", "poplar endophyte")) %>%
  mutate(ecology=as.factor(ecology))

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
  read_tsv("repeatmasker/repeat_tb.tsv")

## Resistomes
rgi_tb <- 
  read_tsv("rgi/rgi_tb.tsv")

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
  left_join(select(repeat_tb, spcode, TE=total_percent)) %>%
  mutate(`duplicated BUSCOs (%)`= 100*`Complete and duplicated BUSCOs (D)`/`Total BUSCO groups searched`) %>%
  mutate(TE_len=TE*genome_size/100) %>%
  mutate(genome_size=genome_size/1e6) %>%
  mutate(TE_len=TE_len/1e6)

write_csv(genome_allstat_tb, "genome stats.csv")

# define color codes
color_scale <-
  c("poplar endophyte"="#F8766D",
    "endophyte"="red",
    "saprobe"="#00BFC4",
    "canker"="purple",
    "ERMF"="gold",
    "fungal parasite"="black",
    "pathogen"="green",
    "ECM"="brown"
  )

ecology_order <-
  c("poplar endophyte",
    "endophyte",
    "ERMF",
    "ECM",
    "pathogen",
    "canker",
    "saprobe",
    "fungal parasite"
  )

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

ggsave("phylogeny.pdf",
       all_phylo_spcode_plot,
       width = 10, height = 6)

# genomic stats contrast and phylogenetic comparison
input_vars <- c("CAZy", "genome_size", "gene_count", "SMC", "n_secreted", "TE_len", "resistome")

contrast_tb <-
  eco_tb %>%
  left_join(select(taxa_tb, spcode, phylum, class),
            by="spcode") %>%
  select(-phylum, -class) %>%
  arrange(ecology) %>%
  mutate(ecology=fct_collapse(
    ecology,
    `poplar endophyte`=c("poplar endophyte"),
    saprobe="saprobe",
    other_level = "Other"
  )
  ) %>%
  filter(ecology!="Other") %>%
  mutate(ecology=fct_drop(ecology))

contrast_tb <-
  contrast_tb %>%
  filter(!(spcode %in% c("Amore1", "Hesve2finisherSC"))) %>% #extra saprobe
  filter(!(spcode %in% c("Gloci1"))) %>% #dulicate in one lineage
  filter(!(spcode %in% c("Ilysp1","MarPMI226"))) %>% #dulicate in one lineage
  filter(!(spcode %in% c("Hyafin1"))) #dulicate in one lineage

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

## Wilcoxon test
wilcox_wrapper <-
  function(.data, pre_col, res_col){
    pre_col <- rlang::sym(pre_col)
    res_col <- rlang::sym(res_col)
    
    .data <- filter(.data, !is.na(!!pre_col), !is.na(!!res_col))
    x <- pull(.data, !!pre_col) %>% as.factor()
    y <- pull(.data, !!res_col)
    
    w_results <- wilcox.exact(y ~ x)
    
    tibble(var=as.character(res_col), W=w_results$statistic, p=w_results$p.value)
  }

contrast_w_tb <-
  map_dfr(input_vars, ~wilcox_wrapper(.data = constrast_tree_dt$data, pre_col = "ecology", .x)) %>%
  arrange(p)

write_csv(contrast_w_tb, "genome Wilcoxon rank.csv")

## phylogenetic least square regression
phylolm_wrapper <-
  function(.data, id_col, pre_cols, res_cols, phy, ...) {
    .data <- column_to_rownames(.data, id_col)
    in_formula <-
      formulate_vars(res_vars=res_cols, pre_vars=pre_cols)
    phylolm(in_formula, .data, phy, ...)
  }

extract_phylolm_fit <-
  function(fit) {
    var_name <- fit$formula %>% terms() %>% attr("variables") %>% .[[2]] %>% as.character()
    summary(fit) %>%
      .$coefficients %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      mutate(var=var_name) %>%
      select(var, term, everything())
  }

contrast_phylolm_fits <-
  map(input_vars, ~phylolm_wrapper(.data = constrast_tree_dt$data, id_col = "spcode", pre_cols = "ecology", res_cols = .x, phy = constrast_tree_dt$tree, model="OUfixedRoot"))

contrast_phylolm_fits <-
  map(input_vars,
      ~phylolm_wrapper(.data = constrast_tree_dt$data, id_col = "spcode", pre_cols = "ecology", res_cols = .x, phy = constrast_tree_dt$tree))

contrast_phylolm_tb <-
  map_dfr(contrast_phylolm_fits, extract_phylolm_fit) %>%
  filter(term!="(Intercept)") %>%
  select(-term)

write_csv(contrast_phylolm_tb, "genome phylogenetic linear regression.csv")

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


constrast_tree_dt$data

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

ggsave("genome contrast boxplot plot.pdf",
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

ggsave("genome contrast treecol plot.pdf",
       treecol_plot,
       width=16,
       height=14,
       units="cm")

# PFAM
## read data table
pfam_rawtb <-
  read_csv("pfam/pfam.csv")

pfam_tb <-
  pfam_rawtb %>%
  select(-`Annotation Description`) %>%
  mutate_all(~ifelse(is.na(.x), 0, .x)) %>%
  rename(pfam=1) %>%
  transpose_tibble(index_col = "pfam", new_var ="spcode" ) %>%
  standarlize(col_excluded = "spcode", scale = F) %>%
  column_to_rownames("spcode")

pfam_tb <-
  pfam_tb %>%
  rownames_to_column("spcode") %>%
  as_tibble()

contrast_asco_tb <-
  contrast_tb %>%
  left_join(select(taxa_tb, spcode, phylum)) %>%
  filter(phylum=="Ascomycota") %>%
  select(-phylum)

contrast_allecoasco_tb <-
  eco_tb %>%
  left_join(select(taxa_tb, spcode, phylum, class),
            by="spcode") %>%
  filter(phylum %in% c("Ascomycota")) %>%
  select(-phylum, -class) %>%
  arrange(ecology, spcode)

contrast_allecoasco_pfam_tb<-
  contrast_allecoasco_tb %>%
  inner_join(pfam_tb) %>%
  column_to_rownames("spcode")

pfam_raw_indata <-
  pfam_tb %>%
  right_join(contrast_asco_tb) %>%
  arrange(ecology, spcode) %>%
  column_to_rownames("spcode")

pfam_indata <-
  pfam_raw_indata %>%
  standarlize("ecology", scale=F)

pfam_vars <- colnames(select(pfam_indata, -ecology))

## wilcoxon test
w_test_tb <-
  map_dfr(pfam_vars, ~wilcox_wrapper(.data = pfam_indata, pre_col = "ecology", .x)) %>%
  arrange(p)

w_summary <-
  w_test_tb %>%
  arrange(p) %>%
  filter(p<0.0001) %>%
  mutate_at(-1, ~round(.x, digits = 5)) %>%
  left_join(select(pfam_rawtb, var=`Annotations/Genomes`, annotation=`Annotation Description`))

w_vars <- w_test_tb %>% arrange(p) %>% filter(p<0.0001) %>% pull(var)

write_csv(w_summary, "pfam wilcoxon.csv")

## PLR
pfam_treedata <-
  pfam_indata %>%
  rownames_to_column("spcode") %>%
  sub_treeXdata("spcode", raxml_phylo_fungi)

pfam_treedata$tree$node.label <- NULL

plr_fits <-
  map(pfam_vars, ~phylolm_wrapper(.data = pfam_treedata$data, id_col = "spcode", pre_cols = "ecology", res_cols = .x, phy = pfam_treedata$tree, model="OUfixedRoot"))

plr_tb <-
  map_dfr(plr_fits, extract_phylolm_fit) %>%
  filter(term!="(Intercept)") %>%
  select(-term) %>%
  arrange(`p.value`)

plr_summary_tb <-
  plr_tb %>%
  filter(`p.value`<0.0001) %>%
  select(var, t.value, p.value) %>%
  mutate_at(-1, ~round(.x, digits = 5)) %>%
  left_join(select(pfam_rawtb, var=`Annotations/Genomes`, annotation=`Annotation Description`))

plr_vars <-
  plr_summary_tb %>%pull(var)

write_csv(plr_summary_tb, "pfam PLR.csv")

## random forest
set.seed(1)
fitControl <- trainControl(method = "LOOCV")

rf_loocv_vita <- train(ecology ~ .,
                       data=pfam_indata,
                       method="ranger",
                       trControl=fitControl,
                       importance="impurity_corrected"
)

rf_summary_tb <-
  rf_loocv_vita$finalModel %>%
  importance_pvalues() %>%
  as_tibble(rownames = "var") %>%
  arrange(pvalue) %>%
  filter(pvalue<0.0001)

rf_vars <-
  rf_summary_tb %>% pull(var)

write_csv(rf_summary_tb, "pfam random forest.csv")

## summary
pfam_summary <-
  tibble(var=w_vars, w=TRUE) %>%
  full_join(tibble(var=rf_vars, rf=TRUE)) %>%
  full_join(tibble(var=plr_vars, PLR=TRUE)) %>%
  mutate_all(.funs = ~replace_na(., FALSE)) %>% 
  rowwise() %>%
  mutate(n_method=sum(c_across(-1))) %>%
  arrange(desc(n_method)) %>%
  rename(
    Wilcoxon=w,
    `Random Forest`=rf,
    pfam=var
  ) %>%
  left_join(select(pfam_rawtb, pfam=`Annotations/Genomes`, annotation=`Annotation Description`))

write_csv(pfam_summary, "pfam summary.csv")

ven_core_plot <-
  ggVennDiagram(
    list(
      Wilcoxon=w_vars,
      `Random forest`=rf_vars,
      PLR=plr_vars),
    label="count",
    edge_lty=0) +
  scale_fill_gradient(low="white", high="white") +
  xlim(-8, 8) +
  theme(legend.position = "none")

ggsave("pfam Venn Diagram.pdf", ven_core_plot)

core_pfam <-
  pfam_summary %>%
  filter(n_method>=2) %>%
  pull(pfam)

## plots
heat_annotation_color <-
  list(ecology=c(`poplar endophyte`="#00BFC4",
                 `saprobe`="#F8766D"),
       Wilcoxon=c(`TRUE`="black", `FALSE`="white"),
       `Random Forest`=c(`TRUE`="black", `FALSE`="white"),
       `PLR`=c(`TRUE`="black", `FALSE`="white")
  )

pfam_methods <-
  pfam_summary %>%
  select(-n_method, -annotation) %>%
  mutate_at(vars(-pfam), ~as.character(.x)) %>%
  column_to_rownames("pfam")

asco_core_heat <-
  pfam_indata %>%
  select(all_of(core_pfam)) %>%
  pheatmap(cluster_rows = F, cluster_cols = T, annotation_row = select(pfam_indata, ecology), annotation_col = pfam_methods, scale = "column",clustering_distance_rows = "manhattan", clustering_distance_cols = "manhattan", annotation_colors = heat_annotation_color, method="ward.D2")

ggsave("Asco_core_pfam.png", asco_core_heat, width = 12, height = 8)

## pfam clustering
### co-phylogeny
allecoasco_indata <-
  pfam_tb %>%
  inner_join(contrast_allecoasco_tb) %>%
  mutate(ecology=factor(ecology, levels=ecology_order)) %>%
  column_to_rownames("spcode") %>%
  standarlize("ecology", scale=F) %>%
  select(ecology, all_of(core_pfam))

genome_hclust <-
  allecoasco_indata %>%
  select(-ecology) %>%
  scale() %>%
  dist() %>%
  hclust(method="ward.D2")

genome_phylo <-
  genome_hclust %>%
  as.phylo()

color_def <-
  color_scale %>%
  as.data.frame() %>%
  set_colnames("color") %>%
  rownames_to_column("ecology")

raxml_phylo <-
  keep.tip(raxml_phylo_fungi, rownames(allecoasco_indata))

tip_ecology <-
  allecoasco_indata %>%
  rownames_to_column("spcode") %>%
  select(spcode, ecology) %>%
  left_join(color_def) %>%
  mutate(spcode=ordered(spcode, levels=raxml_phylo$tip.label)) %>%
  arrange(spcode)

name_raxml_phylo <-
  raxml_phylo %>%
  as_tibble() %>%
  left_join(select(all_genome_tb, label=spcode, Name)) %>%
  select(-label) %>%
  rename(label=Name) %>%
  as.phylo()

name_genome_phylo <-
  genome_phylo %>%
  as_tibble() %>%
  left_join(select(all_genome_tb, label=spcode, Name)) %>%
  select(-label) %>%
  rename(label=Name) %>%
  as.phylo()

name_tip_ecology <-
  tip_ecology %>%
  left_join(select(all_genome_tb, spcode, Name)) %>%
  select(-spcode) %>%
  rename(label=Name)

cophyloplot(name_genome_phylo, name_raxml_phylo, assoc=cbind(name_raxml_phylo$tip.label,name_raxml_phylo$tip.label), show.tip.label = T, space=1000, length.line=0, col=pull(name_tip_ecology, color), font=1, gap=150)


text(x = 10, y=-1, labels = "Phylogeny")
text(x = 150, y=-1, labels = "Pfam clustering")
dev.copy(pdf, "pfam cophylo.pdf")
dev.off()

### heatmap
pfam_hclust <-
  allecoasco_indata %>%
  rownames_to_column("spcode") %>%
  select(-ecology) %>%
  transpose_tibble(index_col = "spcode", "pfam") %>%
  column_to_rownames("pfam") %>%
  scale() %>%
  dist() %>%
  hclust(method="ward.D2")

heat_annotation_color_method <-
  list(Wilcoxon=c(`TRUE`="black", `FALSE`="white"),
       PLR=c(`TRUE`="black", `FALSE`="white"),
       `Random Forest`=c(`TRUE`="black", `FALSE`="white")
  )

heat_annotation_color_ecology <-
  list(ecology=color_scale)

heat_annotation_color <-
  c(heat_annotation_color_method, heat_annotation_color_ecology)

allecoasco_heat <-
  allecoasco_indata %>%
  select(-ecology) %>%
  pheatmap(cluster_rows = genome_hclust, cluster_cols = pfam_hclust, annotation_row = select(allecoasco_indata, ecology), annotation_col = pfam_methods, scale = "column", annotation_colors = heat_annotation_color, fontsize=8)


ggsave("allasco_pfam_heatmap.pdf", allecoasco_heat, width = 6, height=8)

allasco_heat_tree <-
  keep.tip(raxml_phylo_fungi, rownames(allecoasco_indata)) %>%
  as_tibble() %>%
  left_join(eco_tb, by=c("label"="spcode")) %>%
  mutate(ecology=fct_relevel(ecology, "poplar endophyte", "endophyte")) %>%
  as.treedata() %>%
  ggtree() +
  geom_tiplab(aes(color=ecology))

asco_phyloheatmap <-
  gheatmap(allasco_heat_tree,
           data = scale(allecoasco_indata[,-1]),
           offset = 0.2,
           colnames_offset_y=-0.4,
           colnames_angle=90,
           font.size=1.2) +
  scale_color_manual(values=heat_annotation_color_ecology$ecology,  name="ecology")

ggsave("asco_phyloheatmap.pdf", asco_phyloheatmap)

# CAZyme
## parse cazyme to hierarchy
### define cazyme names
cazyme_fam <- c("AA", "CBM", "CE", "GH", "GT", "PL", "EXPN") #main cazyme family(=level1)

#### level 1
all_cazyme <- cazyme_rawtb %>% pull(cazyme)
cazyme_def <- tibble(cazyme_1=cazyme_fam)

#### level 2
cazyme_def <- cazyme_def %>% #names having number after main family, e.g., AA1, AA2...
  mutate(cazyme_2=map(cazyme_1, ~grep(paste0("^", .x, "\\d+$"),all_cazyme, value=TRUE))) 

cazyme_def <- cazyme_def %>% #if no such names, inherit from level1
  rowwise() %>% mutate(cazyme_2=ifelse(length(cazyme_2)==0, list(cazyme_1), list(cazyme_2)))
cazyme_def <- cazyme_def %>% unnest(cazyme_2)

#### level 3
cazyme_def <- cazyme_def %>% #names having underline after level2, e.g., AA1_1, AA1_dist...
  mutate(cazyme_3=map(cazyme_2, ~grep(paste0("^", .x, "_.+$"), all_cazyme, value=TRUE)))

cazyme_def <- cazyme_def %>% #if no such names, inherit from level2
  rowwise() %>% mutate(cazyme_3=ifelse(length(cazyme_3)==0, list(cazyme_2), list(cazyme_3)))
cazyme_def <- cazyme_def %>% unnest(cazyme_3)

#### cazyme group
c_lv1 <- cazyme_def %>% pull(cazyme_1) %>% unique()
c_lv2 <- cazyme_def %>% pull(cazyme_2) %>% unique()
c_lv3 <- cazyme_def %>% pull(cazyme_3) %>% unique()

cazyme_indata <-
  cazyme_tb %>%
  select(spcode, c_lv2) %>%
  mutate_all(~ifelse(is.na(.x), 0, .x)) %>%
  standarlize(col_excluded = "spcode", scale = F) %>%
  inner_join(contrast_asco_tb) %>%
  select(spcode, ecology, everything()) %>%
  arrange(ecology, spcode) %>%
  column_to_rownames("spcode")

cazyme_vars <- colnames(select(cazyme_indata, -ecology))

## wilcoxon test
w_test_cazyme_tb <-
  map_dfr(cazyme_vars, ~wilcox_wrapper(.data = cazyme_indata, pre_col = "ecology", .x)) %>%
  arrange(p)

w_cazyme_summary <-
  w_test_cazyme_tb %>%
  arrange(p) %>%
  filter(p<0.01) %>%
  mutate_at(-1, ~round(.x, digits = 5))

w_cazyme_vars <- w_cazyme_summary %>% pull(var)

write_csv(w_cazyme_summary, "cazyme wilcoxon.csv")

## random forest
set.seed(0)
fitControl <- trainControl(method = "LOOCV")

rf_cazyme_vita <- train(ecology ~ .,
                        data=cazyme_indata,
                        method="ranger",
                        trControl=fitControl,
                        importance="impurity_corrected"
)

rf_cazyme_tb <-
  rf_cazyme_vita$finalModel %>%
  importance_pvalues() %>%
  as_tibble(rownames = "var") %>%
  arrange(pvalue, desc(importance))

rf_cazyme_summary <-
  rf_cazyme_tb %>%
  filter(pvalue<0.01)

rf_cazyme_vars <-
  rf_cazyme_summary %>% pull(var)

write_csv(rf_cazyme_summary, "cazyme random forest.csv")

## PLR
cazyme_treedt <-
  select(cazyme_indata, all_of(cazyme_vars)) %>%
  rownames_to_column("spcode") %>%
  left_join(contrast_asco_tb) %>%
  arrange(ecology) %>%
  sub_treeXdata("spcode", raxml_phylo_fungi)

cazyme_treedt$data <-
  cazyme_treedt$data %>%
  as_tibble() %>%
  standarlize(col_excluded = c("spcode", "ecology"), scale = F) %>%
  as.data.frame()

cazyme_vars <-
  cazyme_treedt$data %>%
  select(-spcode, -ecology) %>%
  colnames()

cazyme_phylolm_fits <-
  map(cazyme_vars, ~phylolm_wrapper(.data = cazyme_treedt$data, id_col = "spcode", pre_cols = "ecology", res_cols = .x, phy = cazyme_treedt$tree, model="OUfixedRoot"))

cazyme_phylolm_tb <-
  map_dfr(cazyme_phylolm_fits, extract_phylolm_fit) %>%
  filter(term!="(Intercept)") %>%
  select(-term) %>%
  arrange(`p.value`)

cazyme_phylolm_summary <-  
  cazyme_phylolm_tb %>%
  mutate_at(-1, ~round(.x, 3))

cazyme_phylolm_vars <- cazyme_phylolm_summary %>% pull(var)

write_csv(cazyme_phylolm_summary, "cazyme PLR.csv")

