library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]

read_busco_summary <- function(infile){
  read_lines(infile) %>%
    grep("^\t[0-9]", ., value = TRUE) %>%
    read_tsv(col_names = FALSE) %>%
    select(X2, X3) %>%
    pivot_wider(values_from = "X2", names_from = "X3")
}  

in_file_v <- scan(in_file, what = character())

busco_tb <-
  in_file_v %>%
  map_chr(~paste0("busco/", .x, "_fungi/run_fungi_odb10/short_summary.txt")) %>%
  map_dfr(read_busco_summary) %>%
  mutate(spcode=in_file_v) 

busco_tb <-
  busco_tb %>%
  mutate(completeness=`Complete BUSCOs (C)`/`Total BUSCO groups searched`) %>%
  select(spcode, completeness, everything())

write_tsv(busco_tb, out_file)
