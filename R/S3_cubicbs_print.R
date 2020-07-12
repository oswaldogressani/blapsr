#' @method print cubicbs
#' @export


# Print method for cubicbs class
print.cubicbs <- function(x, ...) {
  cat("Cubic B-spline basis:\n")
  cat("\n")
  cat("Basis lower bound = ", x$lower, "\n")
  cat("Basis upper bound = ", x$upper, "\n")
  cat("Number of B-splines in basis = ", x$K, "\n")
  cat("Number of knots = ", x$nknots, "\n")
  cat("\n")
  cat("Knots:\n")
  cat(format(round(x$knots,4), nsmall = 4), "\n")
  cat("\n")
  cat("B-spline basis matrix has", x$dimbasis[1], "row(s) and",
      x$dimbasis[2], "columns. \n")
  cat("\n")
  cat("B-spline matrix:\n")
  print(round(x$Bmatrix, 5))
}



