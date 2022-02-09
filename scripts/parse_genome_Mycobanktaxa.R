#' search genome classification through mycobank
#'
#' search the first two words of strain names in "SMC/JGI_SMC.csv" against mycobank
#' by search_mycobank.R
library(tidyverse)
library(magrittr)
source("scripts/search_mycobank.R")

db_path="MBList.xlsx"

in_database <- read_xlsx(
  db_path,
  sheet=1
)

in_genome <- read_csv("meta_tables/genome_relatives.csv") %>%
  select(spcode=portal_ID, species)

mycobank_results <-
  in_genome %>%
  pull(search_taxa) %>%
  search_mycobank(in_database=in_database)

outable <-
  in_genome %>%
  select(-search_taxa) %>%
  bind_cols(select(mycobank_results, -in_taxa))

write_csv(outable, "meta_tables/genome_Mycobanktaxa.csv")
