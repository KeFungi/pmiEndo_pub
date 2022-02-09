library(tidyverse)

make_rgi_path <-
  function(spcode) {
    paste0("rgi/", spcode, ".txt")
  }

read_rgi <-
  function(in_path) {
    read_tsv(in_path)
  }

spcode_vec <- scan("meta_tables/rel_genome_list.txt", what = character())

rgi_tb <-
  map(spcode_vec, ~mutate(read_rgi(make_rgi_path(.x)), spcode=.x)) %>%
  map(~select(.x, spcode, `AMR Gene Family`)) %>%
  map(~as.data.frame(table(.x))) %>%
  map_dfr(~pivot_wider(.x, names_from = 2, values_from=3)) %>%
  rowwise() %>%
  mutate(total=sum(c_across(2:ncol(.)), na.rm = TRUE))

write_tsv(rgi_tb, "rgi/rgi_tb.tsv")
