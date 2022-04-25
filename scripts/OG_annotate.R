library(tidyverse)

remove_na <-
  function(x){
    x[x!="-"&!is.na(x)]
  }

parse_IPR <-
  function(.tb){
    .tb %>%
      rename(OGID=X1) %>%
      group_by(OGID) %>%
      summarise(
        Pfam=list(ifelse(X4=="Pfam", X5, NA)),
        IPR=list(X12) 
      ) %>%
        rowwise() %>%
        mutate_at(-1, ~list(remove_na(.x))) %>%
        mutate_at(-1, ~list(unique(.x))) %>%
        mutate_at(-1, ~paste(.x, collapse = ";"))
  }

longer_OGpfam <-
  function(.tb){
    .tb %>%
      select(pfam=Pfam, OGID) %>%
      mutate(pfam=str_split(pfam, ";")) %>%
      unnest(pfam) %>%
      group_by(pfam) %>%
      summarise(OGID=paste0(OGID, collapse = ";"))
  }
  
read_tsv("ortho_anno/interproscan.tsv", col_names = FALSE) %>%
  parse_IPR() %>%
  write_csv("results/ortho_InterProScan.csv")

read_tsv("meta_tables/arabidopsis_interproscan.tsv", col_names = FALSE) %>%
  parse_IPR() %>%
  write_csv("results/arabidopsis_InterProScan.csv")

ortho_pfam <-
  read_csv("results/ortho_InterProScan.csv") %>%
  longer_OGpfam()

ara_pfam <-
  read_csv("results/arabidopsis_InterProScan.csv") %>%
  longer_OGpfam()

pfam_tb <-
  read_csv("results/pfam_stat_sum.csv")

left_join(pfam_tb, ortho_pfam) %>%
  write_csv("results/Pfam_sum_OG.csv")

left_join(pfam_tb, ortho_pfam) %>%
  left_join(rename(ara_pfam, araOG=OGID)) %>%
  write_csv("results/Pfam_sum_araOG.csv")

