library(tidyverse)

make_repeat_path <-
  function(spcode) {
    paste0("repeatmasker/", spcode, ".fasta.tbl")
  }

read_repeat_percent <-
  function(in_path) {
    in_lines <- read_lines(in_path)
    
    total_percent <-
      grep("bases masked:", in_lines, value = TRUE) %>%
      str_match("\\( (.+) %\\)") %>%
      .[,2] %>%
      drop() %>%
      as.double()
    return(total_percent)
  }
    
read_repeat <-
  function(in_path, main_class_only=TRUE) {
    in_lines <- read_lines(in_path)
    table_start <- grep("-{2,}", in_lines)
    table_end <- grep("={2,}", in_lines)[3]
    
    in_lines <- in_lines[(table_start+1) : (table_end-1)]
    in_lines <- in_lines[!grepl("^$", in_lines)]
    
    out_tb <- tibble()
    main_class <- character()
    
    for (line in in_lines) {
      sep_line <-
        line %>%
        str_replace_all("^\\W+", "") %>%
        str_replace_all("\\ {2,}", "\t") %>%
        paste0("\n") %>%
        str_split("\t", simplify = TRUE) %>%
        drop()
      
      if (grepl(":$", sep_line[1])) {
        sep_line[1] <- str_replace(sep_line[1], ":$", "")
        main_class <- sep_line[1]
      }
      
      n_element <-
        sep_line %>%
        .[grep("^[0-9]+$", .)]
      
      len <-
        sep_line %>%
        .[grep(" bp$", .)] %>%
        str_replace(" bp$", "") %>%
        drop()
      
      percentage <-
        sep_line %>%
        .[grep(" %\n", .)] %>%
        str_replace(" %\n", "") %>%
        drop()
      
      out_row <- tibble(main_class=main_class,
                        class=sep_line[1],
                        element=n_element,
                        length=len,
                        percentage=percentage
                        )
      out_row
      out_tb <- bind_rows(out_tb, out_row)
    }
    if (main_class_only) {
      out_tb <- out_tb[out_tb[["main_class"]]==out_tb[["class"]],]
      out_tb <- select(out_tb, -class)
    }
    return(out_tb)
  }

read_repeat_percent

spcode_vec <- scan("meta_tables/rel_genome_list.txt", what = character())

repeat_percent <- 
  spcode_vec %>%
  map(make_repeat_path) %>%
  map_dbl(read_repeat_percent)

repeat_tb <-
  spcode_vec %>%
  map(make_repeat_path) %>%
  map(read_repeat) %>%
  map(~pivot_longer(.x, element:percentage, )) %>%
  map(~pivot_wider(.x, names_from=`main_class`:`name`)) %>%
  map2_dfr(spcode_vec, ~mutate(.x, spcode=.y)) %>%
  mutate(total_percent=repeat_percent) %>%
  select(spcode, total_percent, everything())

write_tsv(repeat_tb, "repeatmasker/repeat_tb.tsv")
