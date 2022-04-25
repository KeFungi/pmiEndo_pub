source("shared_setup.R")
# Pfam
# read data table
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

# wilcoxon test

w_test_tb <- tibble()
for (var in pfam_vars){
  w_results <-
    wilcox.exact(pull(pfam_indata, var) ~ pull(pfam_indata, "ecology"))
  out_tb <- 
    tibble(var=var, W=w_results$statistic, p=w_results$p.value)
  w_test_tb <- bind_rows(w_test_tb, out_tb)
}

w_summary <-
  w_test_tb %>%
  arrange(p) %>%
  filter(p<0.0001) %>%
  mutate_at(-1, ~round(.x, digits = 5))

w_vars <- w_test_tb %>% arrange(p) %>% filter(p<0.0001) %>% pull(var)

write_csv(w_summary, "results/pfam_w.csv")

## PLR
pfam_treedata <-
  pfam_indata %>%
  rownames_to_column("spcode") %>%
  sub_treeXdata("spcode", raxml_phylo_fungi)

pfam_treedata$tree$node.label <- NULL

plr_tb <- tibble()
for (var in pfam_vars){
  in_formula <-
    formulate_vars(res_vars=var, pre_vars="ecology")
  fit <-
    phylolm(in_formula, select(pfam_indata, "ecology", var), pfam_treedata$tree, model="OUfixedRoot")
  out_tb <-
    summary(fit) %>%
    .$coefficients %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(var=var) %>%
    filter(term!="(Intercept)") %>%
    select(var, everything())
  
  plr_tb <-
    bind_rows(plr_tb, out_tb)
}

plr_summary_tb <-
  plr_tb %>%
  filter(`p.value`<0.0001) %>%
  select(var, t.value, p.value) %>%
  mutate_at(-1, ~round(.x, digits = 5))

plr_vars <-
  plr_summary_tb %>% pull(var)

write_csv(plr_summary_tb, "results/pfam_plr.csv")

## random forest
set.seed(1)
fitControl <- trainControl(method = "LOOCV")

rf_loocv_vita <-
  train(
    select(pfam_indata, -ecology),
    pull(pfam_indata, ecology),
    method="ranger",
    trControl=fitControl,
    importance="impurity_corrected"
  )

rf_loocv_vita$bestTune

rf_summary_tb <-
  rf_loocv_vita$finalModel %>%
  importance_pvalues() %>%
  as_tibble(rownames = "var") %>%
  arrange(pvalue) %>%
  filter(pvalue<0.0001)

write_csv(rf_summary_tb, "results/pfam_rf.csv")

## summary
w_summary_tb <- read_csv("results/pfam_w.csv")
plr_summary_tb <- read_csv("results/pfam_plr.csv")
rf_summary_tb <- read_csv("results/pfam_rf.csv")

w_vars <- w_summary_tb %>% pull(var)
plr_vars <- plr_summary_tb %>% pull(var)
rf_vars <- rf_summary_tb %>% pull(var)

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

write_csv(pfam_summary, "results/pfam_stat_sum.csv")

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

ggsave("plots/pfam_Venn.pdf", ven_core_plot)

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

ggsave("plots/Asco_core_pfam.png", asco_core_heat, width = 12, height = 8)

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
dev.copy(pdf, "plots/pfam_cophylo.pdf")
dev.off()

## heatmap
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

ggsave("plots/pfam_allasco_heatmap.pdf", allecoasco_heat, width = 6, height=8)

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

ggsave("plots/pfam_asco_phyloheatmap.pdf", asco_phyloheatmap)

## PCA
pfam_pca <-
  allecoasco_indata %>%
  do_PCA("ecology") %>%
  ggplot(aes(x=PC1, y=PC2, color=ecology)) +
  geom_point() +
  scale_color_manual(values=heat_annotation_color_ecology$ecology,  name="ecology") +
  theme_classic()

ggsave("plots/pfam_pca.pdf", pfam_pca)
