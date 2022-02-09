require(ape)
sub_treeXdata <-
  function(.data.frame, label_col, .phylo) {
    label_col <- rlang::ensym(label_col)
    label_order <- .phylo$tip.label
    data_label <- pull(.data.frame, !!label_col)
    
    if (anyDuplicated(data_label)) {
      stop("label_col has duplicated value")
    }
    
    if (length(label_order) != length(data_label)) {
      warning("labels in data.frame and phylo do not match, try inner join")
    } else {
      if (!all(sort(label_order)==sort(data_label))){
      warning("labels in data.frame and phylo do not match, try inner join")
      }
    }
    
    inner_label <- intersect(label_order, data_label)
    
    inner_tree <- keep.tip(.phylo, inner_label)
    inner_label_order <- inner_tree$tip.label
    
    inner_tb <- filter(.data.frame, !!label_col %in% inner_label)
    inner_tb <- arrange(inner_tb, match(!!label_col, inner_label_order))
    
    return(list(tree=inner_tree, data=inner_tb, label_col=as.character(label_col)))
  }

pic_treeXdata <-
  function(treeXdata) {
    label_col <- treeXdata$label_col
    map_dfc(dplyr::select(treeXdata$data, -!!label_col), pic, phy=treeXdata$tree)
  }

to_treedata <-
  function(treeXdata) {
    intree_tb <- as_tibble(treeXdata$tree)
    full_join(intree_tb, treeXdata$data, by=c("label" = treeXdata$label_col)) %>%
      as.treedata()
}
