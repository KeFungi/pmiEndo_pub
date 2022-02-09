library(tidyverse)
library(magrittr)
filelist <- list.files("raw_pfam")
csvlist <- grep(".+\\.csv", filelist, value = TRUE)
csv_talbe_list <- map(csvlist, ~read_csv(.x))
csv_talbe_list %<>% map(~select(.x, -starts_with("Total")))
csv_talbe_list %<>% map(~filter(.x, `Annotations/Genomes`!="PFAMs"))

meta_csv <- reduce(csv_talbe_list, full_join,
                   by=c("Annotations/Genomes"="Annotations/Genomes",
                        "Annotation Description"="Annotation Description"))
remove(csv_talbe_list)

remove_join_duplicate_col <- function(.tibble){
  pre_cols <- colnames(.tibble) %>% str_match("[^\\.]+") %>% .[,1]
  dup_cols <- pre_cols %>% (function(x) {duplicated(x) | duplicated(x, fromLast=TRUE)})(.)
  dup_pres <- pre_cols[dup_cols] %>% unique()
  new_cols <- lst()
  for (dup in dup_pres){
    treat_col <- colnames(.tibble) %>% grep(paste0("^", dup), .)
    coalesce_table <- .tibble %>% select(all_of(treat_col))
    new_cols %<>% append(tibble(!!dup:=pmap_dbl(coalesce, .l = coalesce_table)))
  }
  
  .tibble %>% .[, !dup_cols] %>% bind_cols(as_tibble(new_cols))
}

remove_join_duplicate_row <- function(.tibble){
  dup_rows <- meta_csv[,1] %>% (function(x) {duplicated(x) | duplicated(x, fromLast=TRUE)})(.)
  treat_rows <- .tibble[dup_rows, ]
  treat_rows %<>% group_by_at(1) %>% nest()
  treat_rows %<>% mutate(data=map(data, ~summarise_all(.x, function(x) x[which(!is.na(x))[1]])))
  treat_rows %<>% unnest(data)
  
  .tibble %>% .[!dup_rows, ] %>% bind_rows(treat_rows)
}
  
meta_csv %<>% remove_join_duplicate_col()
meta_csv %<>% remove_join_duplicate_row()

meta_csv %<>% arrange(str_sub(`Annotations/Genomes`,3))
write_csv(meta_csv, "../pfam.csv", na="")
