source("shared_setup.R")
# Orthogroup
raw_ortho_tb <-
  read_tsv("orthofinder/Orthogroups.GeneCount.tsv")

ortho_tb <-
  raw_ortho_tb %>%
  transpose_tibble(index_col = "Orthogroup", new_var = "spcode.aa") %>%
  mutate(spcode=str_replace(spcode.aa, "\\.aa$", "")) %>%
  select(-spcode.aa) %>%
  select(spcode, everything()) %>%
  filter(spcode!="Total")

ortho_raw_indata <-
  ortho_tb %>%
  right_join(contrast_asco_tb) %>%
  arrange(ecology, spcode) %>%
  column_to_rownames("spcode")

ortho_indata <-
  ortho_raw_indata %>%
  standarlize("ecology", scale=F)

ortho_vars <- colnames(select(ortho_indata, -ecology))

## wilcoxon test
w_test_tb <- tibble()
for (var in ortho_vars){
  w_results <-
    wilcox.exact(pull(ortho_indata, var) ~ pull(ortho_indata, "ecology"))
  
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

write_csv(w_summary, "results/ortho_w.csv")

## PLR
ortho_treedata <-
  ortho_indata %>%
  rownames_to_column("spcode") %>%
  sub_treeXdata("spcode", raxml_phylo_fungi)

ortho_treedata$tree$node.label <- NULL

plr_tb <- tibble()
for (var in ortho_vars){
  in_formula <-
    formulate_vars(res_vars=var, pre_vars="ecology")
  
  fit <-
    phylolm(in_formula, select(ortho_indata, "ecology", var), ortho_treedata$tree, model="OUfixedRoot")
  
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

write_csv(plr_summary_tb, "results/ortho_plr.csv")

## random forest
set.seed(1)
fitControl <- trainControl(method = "LOOCV")

rf_loocv_vita <-
  train(
    select(ortho_indata, -ecology),
    pull(ortho_indata, ecology),
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

rf_vars <-
  rf_summary_tb %>% pull(var)

write_csv(rf_summary_tb, "results/ortho_rf.csv")


## summary
w_summary_tb <- read_csv("results/ortho_w.csv")
plr_summary_tb <- read_csv("results/ortho_plr.csv")
rf_summary_tb <- read_csv("results/ortho_rf.csv")

w_vars <- w_summary_tb %>% pull(var)
plr_vars <- plr_summary_tb %>% pull(var)
rf_vars <- rf_summary_tb %>% pull(var)

ortho_summary <-
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
  )

write_csv(ortho_summary, "results/ortho_stat_sum.csv")
