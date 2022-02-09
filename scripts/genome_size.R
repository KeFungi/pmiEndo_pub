library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]

in_file <- "meta_tables/rel_genome_list.txt"
out_file <- "meta_tables/genome_size.tsv"


sum_genome_bp <- function(in_file) {
  read_tsv(in_file, col_names = FALSE) %>%
    pull(X2) %>% sum()
}

in_file_v <- scan(in_file, what = character())

genome_size <-
  in_file_v %>%
  map_chr(~paste0("ge/", .x, ".fasta.fai")) %>%
  map_dbl(sum_genome_bp) %>%
  tibble(genome_size=.) %>%
  mutate(spcode=in_file_v) %>%
  select(spcode, everything())

write_tsv(genome_size, out_file)
