library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]


in_file <- "meta_tables/rel_genome_list.txt"
out_file <- "meta_tables/genome_gene_count.tsv"

count_fasta_items <- function(in_file) {
  read_lines(in_file) %>%
    length()
}

in_file_v <- scan(in_file, what = character())

gene_count_tb <-
  in_file_v %>%
  map_chr(~paste0("proteins/", .x, ".aa.fasta.fai")) %>%
  map_dbl(count_fasta_items) %>%
  tibble(gene_count=.) %>%
  mutate(spcode=in_file_v) %>%
  select(spcode, everything())

write_tsv(gene_count_tb, out_file)
