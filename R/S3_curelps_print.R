#' Print the fit of a promotion time cure model.
#'
#' @description Print method for a \code{curelps} object.
#'
#'
#' @param x An object of class \code{curelps}.
#' @param ... Further arguments to be passed to print.
#'
#' @details {Prints informative output of a fitted promotion time cure model
#'  with the Laplace-P-spline approach. In particular, the model formula,
#'  number of B-splines in basis, chosen penalty order, latent field
#'  dimension, sample size, number of events and effective model dimension
#'  are provided. The estimated model coefficients related to the
#'  cure probability (long-term survival) and the population hazard dynamics
#'  (short-term survival) are also provided, where \code{coef} is
#'  the point estimate, \code{sd.post} the posterior standard deviation,
#'  \code{z} is the Wald test statistic and \code{lower.95} and
#'  \code{upper.95} the lower, respectively upper bounds of the approximate
#'  95\% pointwise credible interval.}
#'
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @seealso \code{\link{curelps}}, \code{\link{curelps.extract}},
#'  \code{\link{plot.curelps}}
#'
#' @export

print.curelps <- function(x, ...){
  cat("Formula: \n")
  print(x$formula)
  cat("Object class: \"curelps\" \n")
  cat("\n")
  cat("Number of B-splines in basis:", x$K, sep = " ", "\n")
  cat("Number of parametric coeffs.:", x$p, "\n")
  cat("Latent vector dimension:     ", x$latfield.dim, sep = " ", "\n")
  cat("Penalty order:               ", x$penalty.order, sep = " ", "\n")
  cat("Sample size:                 ", x$n, sep = " ", "\n")
  cat("Number of events:            ", x$num.events, sep = " ", "\n")
  cat("Effective model dimension:   ", format(round(x$ED, digits = 2),
                                              nsmall = 2), sep = " ", "\n")

  cat("\n")
  cat("Coefficients influencing the cure probability (long-term survival): \n")
  print.table(format(round(x$coeff.probacure, digits = 4), nsmall = 4),
              right = TRUE)
  cat("\n")
  cat("Coefficients affecting the population hazard dynamics (short-term survival): \n")
  print.table(format(round(as.matrix(as.data.frame(x$coeff.cox)[, 1:4]),
                           digits = 4), nsmall = 4), right = TRUE)
  cat("--- \n")
  print.table(format(round(as.matrix(as.data.frame(x$coeff.cox)[, 5:8]),
                           digits = 4), nsmall = 4), right = TRUE)
  cat("--- \n")
  cat("AIC.p =",  format(round(x$AIC.p, 4), nsmall = 4),
      "AIC.ED =", format(round(x$AIC.ED, 4), nsmall = 4), "\n")
  cat("BIC.p =", format(round(x$BIC.p, 4), nsmall = 4),
      "BIC.ED =", format(round(x$BIC.ED, 4), nsmall = 4), "\n")
}
