#' transpose tibble through pivot()
#'
#' @param .table input table
#' @param index_col a column to become new column names
#' @new_var new_var column name for the new index column
transpose_tibble <- function(.table, index_col, new_var) {
  table_long <- pivot_longer(.table, -index_col,  names_to = new_var)
  table_wide <-
    pivot_wider(table_long, names_from=index_col, values_from = value)
  return(table_wide)
}

#' filter (near)zero and scale to zscore
#'
#' @param .table input table
#' @param col_excluded col
#' @param freqCut arg to nearZeroVar()
#' @param uniqueCut arg to nearZeroVar()
standarlize <- function(.table, col_excluded, freqCut=Inf, uniqueCut=0, scale=TRUE){
  table_exclude <- select(.table, all_of(col_excluded))
  table_include <- select(.table, -all_of(col_excluded))
  nzv_index <- nearZeroVar(table_include, freqCut=freqCut, uniqueCut=uniqueCut)
  table_include <- select(table_include, -all_of(nzv_index))
  
  if (scale) {
    table_include <- scale(table_include) %>% as_tibble()
  }
  
  return(bind_cols(table_exclude, table_include))
}

#' plot heatmap
plot_pheatmap <-
  function(
    .table, id_col, dep_col,
    cluster_rows = F,
    cluster_cols = T,
    annotation_names_row = F,
    show_colnames = F
  ) {
    id_col <- rlang::sym(id_col)
    dep_col <- rlang::sym(dep_col)
    
    input_matrix <-
      select(.table, -!!dep_col) %>% column_to_rownames(as.character(id_col))
    
    input_color_code <-
      .table %>% select(!!id_col, !!dep_col) %>%
      column_to_rownames(as.character(id_col))
    
    pheatmap(input_matrix,
             annotation_row = input_color_code,
             scale = "column",
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             annotation_names_row = annotation_names_row,
             show_colnames = show_colnames
    )
  }

#' calculate PCA score
#'
#' @param data data matrix
#' @inheritParams table_to_lda
do_PCA <- function(.data, id_cols){
  input_matrix <- select(.data, -all_of(id_cols))
  meta_matrix <- select(.data, all_of(id_cols))
  pca <- prcomp(input_matrix)
  pca_score <-
    bind_cols(meta_matrix, as.data.frame(pca$x))
  return(pca_score)
}

#' plot PCA score and grouping
#'
#' @param pcasocre pca$x
#' @param shape_col second grouping column represented as shapes
#' @inheritParams table_to_lda
plot_PCA <- function(pca_score, id_col, dep_col, shape_col=NA){
  id_col <- rlang::sym(id_col)
  dep_col <- rlang::sym(dep_col)
  if (is.na(shape_col)) {
    PCplot <-
      ggplot(pca_score, aes(x=PC1, y=PC2, group=!!id_col,
                            color=!!dep_col)) +
      geom_point(alpha=0.6, size=0.5) + coord_equal()
  } else {
    shape_col <- rlang::ensym(shape_col)
    PCplot <-
      ggplot(pca_score, aes(x=PC1, y=PC2, group=!!id_col,
                            color=!!dep_col, shape=!!shape_col)) +
      geom_point(alpha=0.6, size=2) +
      scale_shape_manual(values = c(1, 4, 2, 3, 5, 6:20)) +
      coord_equal()
  }
  
  pca_percent <-
    pca_score %>%
    select(matches("^PC[0-9]+$")) %>%
    summarise_all(.funs = var) %>%
    (function(x) x/sum(x))
  
  PCplot <-
    PCplot +
    xlab(paste0("PC1 ", round(pca_percent[1 ,1]*100, 1), "%")) +
    ylab(paste0("PC2 ", round(pca_percent[1, 2]*100, 1), "%"))
  
  return(PCplot)
}

#' print interactive table
#'
#' @param data.frame a data frame
inter_table <-
  function(.df, scrollX=600, scrollY=400, ...) {
    datatable(
      .df,
      extensions = 'Buttons',
      options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv'),
        scrollX=scrollX,
        scrollY=scrollY,
        paging=FALSE,
        searching=FALSE,
        lengthChange=FALSE
      ),
      ...
    )
  }

#' make formula from vectors of variables
#'
#' @param res_vars a vector of response variables
#' @param pre_vars a vector of predict variables
formulate_vars <- function(res_vars, pre_vars){
  pre_term <- paste(pre_vars, collapse ="+")
  res_term <- paste(res_vars, collapse ="+")
  as.formula(paste(res_term, "~", pre_term))
}

univariate_plots <-
  function(.tb, id_col, group_col, var_col) {
    id_col <- rlang::sym(id_col)
    group_col <- rlang::sym(group_col)
    var_col <- rlang::sym(var_col)
    
    #bar plot
    gg1 <- ggplot(.tb, aes(x=reorder(!!id_col, as.integer(as.factor(!!group_col))), y=!!var_col, fill=!!group_col, group=!!id_col)) +
      geom_col() +
      theme_classic()
    
    gg2 <- ggplot(.tb, aes(x=!!group_col, y=!!var_col, fill=!!group_col)) +
      geom_boxplot() +
      theme_classic() +
      ylim(0, NA)
    
    return(list(gg1, gg2))
  }

## phylogenetic least square regression
phylolm_wrapper <-
  function(.data, id_col, pre_cols, res_cols, phy, ...) {
    .data <- column_to_rownames(.data, id_col)
    in_formula <-
      formulate_vars(res_vars=res_cols, pre_vars=pre_cols)
    phylolm(in_formula, .data, phy, ...)
  }

extract_phylolm_fit <-
  function(fit) {
    var_name <- fit$formula %>% terms() %>% attr("variables") %>% .[[2]] %>% as.character()
    summary(fit) %>%
      .$coefficients %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      mutate(var=var_name) %>%
      select(var, term, everything())
  }

wilcox_wrapper <-
  function(.data, pre_col, res_col){
    pre_col <- rlang::sym(pre_col)
    res_col <- rlang::sym(res_col)
    
    .data <- filter(.data, !is.na(!!pre_col), !is.na(!!res_col))
    x <- pull(.data, !!pre_col) %>% as.factor()
    y <- pull(.data, !!res_col)
    
    w_results <- wilcox.exact(y ~ x)
    
    tibble(var=as.character(res_col), W=w_results$statistic, p=w_results$p.value)
  }
