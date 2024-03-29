#' Print a generalized additive model object.
#'
#' @param x An object of class \code{gamlps}.
#' @param ... Further arguments to be passed to print.
#'
#' @details Prints informative output of a fitted generalized additive model.
#'  In particular, the model formula, sample size, number of B-splines in basis,
#'  number of smooth terms, the chosen penalty order,
#'  the latent field dimension, model degrees of freedom, the estimated
#'  coefficients of the linear part and
#'  the estimated effective degrees of freedom of the smooth terms are
#'  displayed.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @seealso \code{\link{gamlps}}, \code{\link{gamlps.object}},
#'  \code{\link{plot.gamlps}}
#'
#' @export

print.gamlps <- function(x, ...) {

  # Function to display significance codes of p-value
  pvalstars <- function(pval){
    char.vec <- c()
    for(j in 1 :length(pval)){
      if(pval[j] <= 0.001) char.vec[j] <- intToUtf8(c(42,42,42))
      else if(pval[j] <= 0.01 && pval[j] > 0.001) char.vec[j] <- intToUtf8(c(42,42))
      else if(pval[j] <= 0.05 && pval[j] > 0.01) char.vec[j] <- intToUtf8(42)
      else if(pval[j] <= 0.1 && pval[j] > 0.05) char.vec[j] <- "."
      else char.vec[j] <- ""
    }
    return(char.vec)
  }

  if(is.null(x$data)) {
    mf <- stats::model.frame(x$formula) # Extract model frame from formula
    X  <- stats::model.matrix(mf)     # Full design matrix
    colXnames <- colnames(X)
    smterms <- grepl("sm(", colnames(X), fixed = TRUE)
    X <- cbind(X[, as.logical(1 - smterms)], X[, smterms])
    colnames(X) <- colXnames
  } else{
    mf <- stats::model.frame(x$formula, data = x$data)
    X <- stats::model.matrix(mf, data = x$data)
    colXnames <- colnames(X)
    smterms <- grepl("sm(", colnames(X), fixed = TRUE)
    X <- cbind(X[, as.logical(1 - smterms)], X[, smterms])
    colnames(X) <- colXnames
  }
  q  <- x$q                    # Number of smooth terms in model
  p  <- ncol(X) - q   # Number of regression coefficients in linear part
  edf.matrix <- matrix(0, nrow = q, ncol = 6)
  edf.matrix[, 1] <- format(round(x$EDf, 4), nsmall = 4)
  edf.matrix[, 2] <- format(round(x$EDfHPD.95[, 1], 4), nsmall = 4)
  edf.matrix[, 3] <- format(round(x$EDfHPD.95[, 2], 4), nsmall = 4)
  edf.matrix[, 4] <- format(round(matrix(x$Approx.signif, ncol = 2)[, 1],
                                  4), nsmall = 4)
  edf.matrix[, 5] <- format.pval(matrix(x$Approx.signif, ncol = 2)[, 2],
                                 digits = 4, nsmall = 4)
  edf.matrix[, 6] <- pvalstars(matrix(x$Approx.signif, ncol = 2)[, 2])
  rownames(edf.matrix) <- colnames(X)[(p + 1):(p + q)]
  colnames(edf.matrix) <- c("edf", "lower.95", "upper.95", "Tr", "p-value", " ")
  cat("Formula: \n")
  print(x$formula)
  cat("\n")
  cat("Family:                      ", x$family, "\n")
  cat("Link function:               ", x$link, "\n")
  cat("Sample size:                 ", x$n, "\n")
  cat("Number of B-splines in basis:", x$K, "\n")
  cat("Number of smooth terms:      ", x$q, "\n")
  cat("Penalty order:               ", x$penalty.order, "\n")
  cat("Latent vector dimension:     ", x$latfield.dim, "\n")
  cat("Model degrees of freedom:    ",
      format(x$ED, nsmall = 2, digits = 2), "\n")
  cat("\n")
  cat("Linear coefficients: \n")
  print.table(format(round(x$linear.coeff, 4), nsmall = 4), right = TRUE)
  cat("--- \n")
  cat("Effective degrees of freedom of smooth terms: \n")
  print.table(format(edf.matrix, nsmall = 4), right = FALSE)
  cat("\n")
  cat(paste0("Signif. codes: 0 ","'",intToUtf8(c(42,42,42)),"' 0.001 ", "'",
             intToUtf8(c(42,42)), "' 0.01 ", "'",intToUtf8(42),
             "' 0.05 '.' 0.1 ' ' 1", sep=" "), "\n")
  cat("--- \n")
  cat("Posterior interval corresponds to a 95% HPD interval \n")
  cat("\n")
  cat("Adjusted R-squared:", format(round(x$r2.adj, 4), nsmall = 4))
}
