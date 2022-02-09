library(tidyverse)
library(rlang)

make_signalp5_path <-
  function(spcode) {
    paste0("secretome/signalp5/", spcode, "_summary.signalp5")
  }

read_signalp5 <-
  function(in_path) {
    header_lines <- read_lines(in_path, n_max = 2)
    col_names <- header_lines[2] %>% str_replace("# ", "") %>%
      str_split("\t") %>% .[[1]]
    out_tb <- read_tsv(in_path, skip=2, col_names = col_names)
    return(out_tb)
  }

make_tmhmm_path <-
  function(spcode) {
    paste0("secretome/tmhmm/", spcode, ".tmhmm")
  }

read_tmhmm <-
  function(in_path) {
    col_names <- c("ID", "len", "ExpAA", "First60", "PredHel", "Topology")
    read_tsv(in_path, col_names = col_names) %>%
      mutate_at(2:6,
                funs({ in_col <- quo(.)
                       col_name <- as_name(in_col)
                       str_replace(., paste0(col_name, "="), "")}
                     )
                )
  }

make_WoLFPSort_path <-
  function(spcode) {
    paste0("secretome/WoLFPSort/", spcode, ".WoLFPSort")
  }

parse_WoLFPSort_content <-
  function(WoLFPSort_content) {
    WoLFPSort_content %>%
      str_split(", ") %>%
      .[[1]] %>%
      str_split(" ") %>%
      map_dfc(~lst(!!.x[1]:=as.double(.x[2])))
  }

read_WoLFPSort <-
  function(in_path) {
    in_tb <-
      read_tsv(in_path, skip = 1, col_names = FALSE) %>%
      mutate(ID=str_match(X1,"^[\\S]+")[,1]) %>%
      mutate(content=str_replace(X1,"^[\\S]+\\s", "")) %>%
      select(-X1) %>%
      mutate(content=map(content, parse_WoLFPSort_content)) %>%
      unnest(content)
    
    in_tb[is.na(in_tb)] <- 0
    col_names <- colnames(in_tb)[-1]
    
    in_tb %>%
      rowwise() %>%
      mutate(max_score=max(c_across(-ID))) %>%
      mutate(prediction=col_names[which.max(c_across(-ID))]) %>%
      select(ID, prediction, max_score, everything())
  }

make_targetp2_path <-
  function(spcode) {
    paste0("secretome/targetp2/", spcode, "_summary.targetp2")
  }

read_targetp2 <-
  function(in_path) {
    header_lines <- read_lines(in_path, n_max = 2)
    col_names <- header_lines[2] %>% str_replace("# ", "") %>%
      str_split("\t") %>% .[[1]]
    out_tb <- read_tsv(in_path, skip=2, col_names = col_names)
    return(out_tb)
  }

make_psscan_path <-
  function(spcode, AC=NULL) {
    if (is.null(AC)) paste0("secretome/ps_scan/", spcode, ".scan")
    else paste0("secretome/ps_scan/", spcode, "_", AC, ".scan")
  }

read_psscan <-
  function(in_path) {
    read_lines(in_path) %>%
      grep(">", . , value=TRUE) %>%
      str_replace("^>", "") %>%
      str_split(" : ", simplify = TRUE) %>%
      as_tibble() %>%
      rename(ID=V1, PROSITE=V2) %>%
      separate(PROSITE, "PROSITE", " ")
  }

make_length_path <-
  function(spcode) {
    paste0(paste0("proteins/", spcode, ".aa.fasta.fai"))
  }

read_length <-
  function(in_path) {
    read_tsv(in_path, col_names = FALSE) %>%
      rename(ID=X1, length=X2)
  }

secretome_pipe <-
  function(spcode) {
    signalp5_tb <-
      spcode %>%
      make_signalp5_path() %>%
      read_signalp5()
    
    tmhmm_tb <-
      spcode %>%
      make_tmhmm_path() %>%
      read_tmhmm()
      
    targetp2_tb <-
      spcode %>%
      make_targetp2_path() %>%
      read_targetp2()
    
    wolf_tb <-
      spcode %>%
      make_WoLFPSort_path() %>%
      read_WoLFPSort()
    
    psscan_tb <-
      spcode %>%
      make_psscan_path("PS00014") %>%
      read_psscan()
    
    length_tb <-
      spcode %>%
      make_length_path() %>%
      read_length()
    
    full_secrotome_tb <-
      select(signalp5_tb, ID, signalp_prediction=Prediction) %>%
      full_join(select(tmhmm_tb, ID, PredHel)) %>%
      full_join(select(targetp2_tb, ID, targetp2_prediction=Prediction)) %>%
      full_join(select(wolf_tb, ID, wolf_prediction=prediction)) %>%
      full_join(select(psscan_tb, ID, PROSITE)) %>%
      full_join(select(length_tb, ID, length)) %>%
      mutate(PROSITE=ifelse(is.na(PROSITE), "none", PROSITE))

    
    message("missing data ", spcode)
    full_secrotome_tb[, 2:7] %>% is.na() %>% which(arr.ind = TRUE) %>% print()
    
    out_tb <-
      full_secrotome_tb %>%
      rowwise() %>%
      mutate(
        secreted=
          ifelse(signalp_prediction != "OTHER" &&
                 PredHel <= 1 &&
                 targetp2_prediction == "SP" &&
                 wolf_prediction == "extr" &&
                 PROSITE == "none" &&
                 length <= 300,
                 1,
                 0
          )
        )
  }

in_file_list <- "meta_tables/rel_genome_list.txt"
spcode_v <- scan(in_file_list, what=character())

summary_tb <-
  tibble(`spcode`=character(), `n_secreted`=integer())

for (spcode in spcode_v) {
  out_tb <-
    secretome_pipe(spcode) %>%
    filter(secreted==1)
  
  write_tsv(out_tb, paste0("secretome/", spcode, ".tsv"))
  summary_tb <- bind_rows(summary_tb, tibble(spcode=spcode, n_secreted=nrow(out_tb)))
}

write_tsv(summary_tb, "secretome/secretome_count.tsv")
