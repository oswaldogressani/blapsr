#' @method print simsurvdata
#' @export

# Print method for simsurvdata class
print.simsurvdata <- function(x, ...) {
  cat("Sample size:         ", x$sample.size, "\n")
  cat("Censoring:           ", x$censoring, "\n")
  cat("Number of events:    ", x$num.events, "\n")
  cat("Censoring percentage:", x$censoring.percentage, "\n")
  cat("Weibull mean:        ", format(round(x$Weibull.mean, 2),
                                      nsmall =2), "\n")
  cat("Weibull variance:    ", format(round(x$Weibull.variance, 2),
                                      nsmall = 2), "\n")
}
