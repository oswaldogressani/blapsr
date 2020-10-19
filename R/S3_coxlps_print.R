#' Print a coxlps object.
#'
#' @description Print method for a \code{coxlps} object.
#'
#' @param x An object of class \code{coxlps}.
#' @param ... Further arguments passed to print.
#'
#' @details Prints informative output of a fitted Cox proportional hazards model
#'  with Laplace-P-splines.
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @seealso \code{\link{coxlps}}
#'
#' @export

# Print method for an object of class coxlaps
print.coxlps <- function(x, ...) {
  cat("Formula: \n")
  print(x$formula)
  cat("Object class: \"coxlps\" \n")
  cat("\n")
  cat("Number of B-splines in basis:", x$K, "\n")
  cat("Number of parametric coeffs.:", x$p, "\n")
  cat("Latent vector dimension     :", x$latfield.dim, "\n")
  cat("Penalty order               :", x$penalty.order, "\n")
  cat("Sample size                 :", x$n, "\n")
  cat("Number of events:           :", x$num.events, "\n")
  cat("Effective dimension (ED)    :", format(round(x$ED, digits = 2),
                                              nsmall = 2), sep = " ", "\n")
  cat("\n")
  cat("Estimated model coefficients: \n")
  cat("\n")
  print.table(format(round(as.matrix(as.data.frame(x$regcoeff)[,1:4]),
                           digits = 4), nsmall = 4), right = TRUE)
  cat("\n")
  print.table(format(round(as.matrix(as.data.frame(x$regcoeff)[,5:8]),
                           digits = 4), nsmall = 4), right = TRUE)
  cat("--- \n")
  cat("AIC.p =",  format(round(x$AIC.p, 4), nsmall = 4),
      "AIC.ED =", format(round(x$AIC.ED, 4), nsmall = 4), "\n")
  cat("BIC.p =", format(round(x$BIC.p, 4), nsmall = 4),
      "BIC.ED =", format(round(x$BIC.ED, 4), nsmall = 4), "\n")
}
