library(tidyverse)
read_busco <-
    function(filepath){
        in_tb <- read_tsv(filepath, skip=3, col_names=c("Busco id","Status","Sequence","Score", "Length","OrthoDB url","Description"))
        return(in_tb)
    }

concat_busco_status <-
	function(busco_tb) {
		busco_tb %>%
                    select(busco_id=`Busco id`, stat=`Status`, gene=`Sequence`) %>%
		    filter(stat=="Complete") %>%
		    rowwise() %>%
		    mutate(cstat=list(c(stat, gene))) %>%
		    select(busco_id, cstat)
	}


spcodel <- scan("../meta_tables/rel_genome_list.txt", what=character())
buscol <- map(paste0(spcodel, "_fungi/run_fungi_odb10/full_table.tsv"), read_busco)
buscocl <- map2(buscol, spcodel, function(x, y){concat_busco_status(x) %>% mutate(spcode=y)})
buscoc_tb <- buscocl %>% reduce(bind_rows) %>% pivot_wider(names_from="spcode", values_from="cstat")
uni_gene <- buscoc_tb %>% select(-busco_id) %>% rowwise() %>% mutate_all(.funs=function(x){is.null(x[[1]])}) %>% rowSums() %>% `==`(0) %>% which()
busco_uni_tb <- buscoc_tb[uni_gene, ]
busco_gene_tb <- busco_uni_tb %>% rowwise() %>% mutate_at(vars(!busco_id), .funs=function(x){x[[2]]})

write_tsv(busco_gene_tb, "busco_fungi_gene_tb.tsv")
