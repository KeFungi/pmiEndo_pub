phylolmstep_fixed <-
function (full.formula, starting.formula, data, phy, trace = 2, k = 2, direction=c("both", "backward", "forward"), ...) 
{
  direction = match.arg(direction)
  response = full.formula[[2]]
  covariates = attr(terms(full.formula), "term.labels")
  p = length(covariates)
  if (p == 0) 
    stop("Your full model only has intercept. Model selection is not needed.")
  create.formula <- function(plm) {
    on = which(plm == 1)
    str = paste(response, " ~ 1", sep = "")
    if (length(on) > 0) 
      for (i in 1:length(on)) str = paste(str, " + ", covariates[on[i]], 
                                          sep = "")
    return(str)
  }
  
  fit <- function(plm) {
    return(
      phylolm(create.formula(plm), data=data, phy=phy, ...)
    )
  }
  
  covariates.current = attr(terms(starting.formula), "term.labels")
  plm.full=rep(1, p)
  plm.current = rep(0, p)
  position = match(covariates.current, covariates)
  
  if (any(is.na(position))) 
    stop("The starting model is not a submodel of the full model.")
  plm.current[position] = 1
  
  fit.current = fit(plm.current)

  if (trace > 0) {
    cat("----------\n")
    cat(paste("Starting model: ", create.formula(plm.current), 
              "\n", sep = ""))
    cat(paste("AIC(k=", k, "): ", AIC(fit.current, k), "\n", 
              sep = ""))
  }
  flag = 0
  count = 1
  while (flag == 0) {
    flag = 1
    plm.best = plm.current
    fit.best = fit.current
    terms.add = add.scope(as.formula(create.formula(plm.current)), 
                          as.formula(create.formula(plm.full)))
    terms.drop = drop.scope(formula(create.formula(plm.current)))
    for (i in 1:p) {
      plm.propose = plm.current
      do.update = FALSE
      if ((plm.current[i] == 1) &&
          (direction %in% c("both","backward")) &&
          (covariates[i] %in% terms.drop)) {
        plm.propose[i] = 0
        do.update = TRUE
      }
      if ((plm.current[i] == 0) &&
          (direction %in% c("both","forward")) &&
          (covariates[i] %in% terms.add)) {
        plm.propose[i] = 1
        do.update = TRUE
      }
      if (do.update) {
        fit.propose = fit(plm.propose)
        if (trace > 1) {
          cat(paste0("\tProposed: ", create.formula(plm.propose), 
                     "\n"))
          cat(paste0("\tAIC(k=", k, "): ", AIC(fit.propose, 
                                               k), "\n"))
        }
        if (AIC(fit.propose, k) < AIC(fit.best, k)) {
          plm.best = plm.propose
          fit.best = fit.propose
          flag = 0
        }
      }
    }
    plm.current = plm.best
    fit.current = fit.best
    if (trace > 0) {
      cat("----------\nStep ", count, "\n", sep = "")
      cat(paste("Current model: ", create.formula(plm.current), 
                "\n", sep = ""))
      cat(paste("AIC(k=", k, "): ", AIC(fit.current, k), 
                "\n", sep = ""))
    }
    count = count + 1
  }
  if (trace > 0) 
    cat("---END\n")
  return(fit.current)
}
