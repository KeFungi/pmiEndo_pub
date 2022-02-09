comparative.data <-
  function (phy, data, names.col, vcv = FALSE, vcv.dim = 2, na.omit = TRUE, 
          force.root = FALSE, warn.dropped = FALSE, scope = NULL) 
{
  phy.name <- deparse(substitute(phy))
  data.name <- deparse(substitute(data))
  if (!is.data.frame(data)) 
    stop("'data' must be an object of class 'data.frame'.")
  namesInd <- match(names.col, names(data))
  if (is.na(namesInd)) {
    stop("Names column '", names.col, "' not found in data frame '", 
         data.name, "'")
  }
  rownames(data) <- as.character(data[, namesInd])
  data <- data[, -namesInd, drop = FALSE]
  if (!inherits(phy, "phylo")) 
    stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
  if (!is.rooted(phy)) {
    if (force.root) {
      phy$root.edge <- 1
    }
    else {
      stop("'", deparse(substitute(phy)), "' is not rooted or has a basal polytomy.")
    }
  }
  if (any(duplicated(phy$tip.label))) 
    stop("Duplicate tip labels present in phylogeny")
  if (any(duplicated(c(phy$tip.label, phy$node.label)))) 
    stop("Labels duplicated between tips and nodes in phylogeny")
  origTips <- with(phy, max(edge) - Nnode)
  origData <- nrow(data)
  in.both <- intersect(rownames(data), phy$tip.label)
  if (length(in.both) == 0) 
    stop("No tips are common to the dataset and phylogeny")
  row.in.tree <- match(rownames(data), in.both)
  row.not.in.tree <- rownames(data)[is.na(row.in.tree)]
  data <- subset(data, !is.na(row.in.tree))
  tip.in.data <- match(phy$tip.label, in.both)
  to.drop <- phy$tip.label[is.na(tip.in.data)]
  if (length(to.drop) > 0) 
    matchedPhy <- drop.tip(phy, to.drop)
  else matchedPhy <- phy
  root <- length(matchedPhy$tip.label) + 1
  tip.order <- match(matchedPhy$tip.label, rownames(data))
  if (any(is.na(tip.order))) 
    stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
  data <- data[tip.order, , drop = FALSE]
  rownames(data) <- matchedPhy$tip.label
  IntNd <- root:max(matchedPhy$edge)
  if (is.null(matchedPhy$node.label)) {
    matchedPhy$node.label <- IntNd
  }
  else {
    matchedPhy$node.label <- ifelse(matchedPhy$node.label == 
                                      "", NA, matchedPhy$node.label)
    if (any(duplicated(na.omit(matchedPhy$node.label)))) 
      stop("Duplicate node labels present in phylogeny")
    matchedPhy$node.label <- ifelse(is.na(matchedPhy$node.label), 
                                    IntNd, matchedPhy$node.label)
  }
  RET <- list(phy = matchedPhy, data = data, data.name = data.name, 
              phy.name = phy.name, dropped = list(tips = to.drop, unmatched.rows = row.not.in.tree))
  class(RET) <- "comparative.data"
  if (vcv) {
    RET$vcv <- VCV.array(matchedPhy, dim = vcv.dim)
    RET$vcv.dim <- vcv.dim
  }
  if (na.omit) {
    before.drop.rows <- rownames(RET$data)
    RET <- na.omit(RET, scope)
    if (!identical(rownames(RET$data), before.drop.rows)) 
      RET$dropped$NA.rows <- before.drop.rows
  }
  if (warn.dropped) {
    if (any(sapply(RET$dropped, length) > 0)) 
      warning("Data dropped in compiling comparative data object")
  }
  return(RET)
}