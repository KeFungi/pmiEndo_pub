#' VERY SLOW hierarchical one-by-one search mycobank
require(tidyverse)
require(readxl)

rank_name_switch <-
  function(rank_name){
    rank_abr <-
      c(
        "p" = "div.",
        "c" = "cl.",
        "o" = "ordo",
        "f" = "fam.",
        "g" = "gen.",
        "s" = "sp."
        )
    
    left_tb <- tibble(input=rank_name)
    right_tb <- as_tibble(rank_abr, rownames="input")
    
    left_join(left_tb, right_tb, by=c("input"="input")) %>%
      mutate(value=ifelse(is.na(value), input, value)) %>%
      pull(value)
  }

search_mycobank <-
  function(v_data, ranks=c("p", "c", "o", "f", "g", "s"), in_database){
    all_ranks <-
      in_database %>%
      pull(`Rank`) %>%
      unique()
    
    s_ranks <- rank_name_switch(ranks)
    
    if (!all(s_ranks %in% all_ranks)) {
      stop("unrecognized rank name")
    }
    
    wk_tb <- tibble(in_taxa=v_data)
    
    c_name_tb <-
      select(
        in_database,
        s_name = `Taxon_name`,
        c_name = `Current name`,
        status = `Name_status`
      ) %>%
      distinct() %>%
      mutate(c_name=ifelse(c_name=="-", s_name, c_name)) %>%
      mutate(status=fct_other(status, keep="Legitimate")) %>%
      group_by(s_name) %>%
      arrange(status) %>%
      slice_head(n=1) %>%
      select(-status)
    
    wk_tb <- wk_tb %>%
      left_join(
        c_name_tb,
        by= c("in_taxa" = "s_name")
      )
    
    s_database <-
      filter(in_database, Name_status == "Legitimate") %>%
      group_by(Taxon_name) %>%
      slice(n=1) %>%
      select(Taxon_name, Classification, rank_name=`Rank`)
    
    search_higher_taxa <-
      function(t_name){
        s_database[match(t_name, s_database[["Taxon_name"]]),] %>%
          pull(Classification) %>%
          str_match(",\\ ([^\\ ]+$)") %>%
          .[, 2] %>%
          drop()
      }
    
    recur_search_taxa <-
      function(start_t){
        t_hier <- lst()
        c_t <- start_t
        
        if(is.na(c_t)) {
          return(lst(!!s_ranks[1]:=NA))
        }
        
        while (!is.na(c_t)){
          c_rank <-
            s_database[match(c_t, s_database[["Taxon_name"]]),] %>%
            pull(rank_name)
          
          if (!c_rank %in% names(t_hier)) {
            t_hier <- append(t_hier, lst(!!c_rank:=c_t))
          }
          
          n_t <- search_higher_taxa(c_t)
          c_t <- n_t
        }
        return(t_hier)
      }
    
    out_col_tb <- tibble()
    for (col_id in s_ranks){
      out_col_tb <- bind_cols(out_col_tb, tibble(!!col_id:=character()))
    }
    
    results_tb <-
      wk_tb %>%
        rowwise() %>%
        mutate(data=map_dfc(c_name, recur_search_taxa)) %>%
        .$data
    
    full_re_tb <- bind_rows(out_col_tb, results_tb)
    
    out_tb <-
      select(wk_tb, in_taxa, current_name=c_name) %>%
      bind_cols(select(full_re_tb, all_of(s_ranks)))
    
    return(out_tb)
  }
