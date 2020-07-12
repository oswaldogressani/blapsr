#' @method plot cubicbs
#' @export

# Plot method for cubicbs class
plot.cubicbs <- function(x, ...) {
  xdom <- seq(x$lower, x$upper, length = 1000)
  internal.knots <- x$knots[x$knots >= x$lower &  x$knots <= x$upper]
  Bsmatrix <- cubicbs(xdom, x$lower, x$upper, x$K)$Bmatrix
  Bsplot <- graphics::plot(xdom, Bsmatrix[, 1], type = "l", col = "black",
                 ylim = c(0, max(Bsmatrix) + 0.2), ylab = "",
                 xlab = "", main = "Cubic B-spline basis"
  )
  for (k in 2:x$K) {
    graphics::lines(xdom, Bsmatrix[, k], type = "l", col = "black")
  }
  graphics::abline(v = internal.knots, lty = 2)
  graphics::rug(x$x, col = "blue")
}
