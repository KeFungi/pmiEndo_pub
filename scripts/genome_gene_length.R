library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]

in_file <- "meta_tables/rel_genome_list.txt"
out_file <- "meta_tables/genome_gene_length.tsv"


average_gene_length <- function(in_file) {
  read_tsv(in_file, col_names = FALSE) %>%
    pull(X2) %>% mean()
}

in_file_v <- scan(in_file, what = character())

gene_length <-
  in_file_v %>%
  map_chr(~paste0("proteins/", .x, ".aa.fasta.fai")) %>%
  map_dbl(average_gene_length) %>%
  tibble(gene_length=.) %>%
  mutate(spcode=in_file_v) %>%
  select(spcode, everything())

write_tsv(gene_length, out_file)
