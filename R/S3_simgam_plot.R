#' @method plot simgam
#' @export

# Plot method for an object of class simgam
plot.simgam <- function(x, ...){

  if(x$dist == "poisson") {
    graphics::plot(x$data$y, xlab = "", ylab = "y",
         main = paste("Poisson data, n =", length(x$data$y)),
         ylim = c(0, max(x$data$y) + 1),
         pch = 16)
  } else if(x$dist == "gaussian") {
    graphics::plot(x$data$y, xlab = "", ylab = "y",
         main = paste("Gaussian data, n =", length(x$data$y)),
         ylim = c(min(x$data$y) - 1, max(x$data$y) + 1),
         pch = 16)
  } else if(x$dist == "bernoulli") {
    graphics::plot(x$data$y, xlab = "", ylab = "y",
         main = paste("Binary data, n =", length(x$data$y)),
         ylim = c(0, 1), pch = 16)
  } else if(x$dist == "binomial") {
    graphics::plot(x$data$y, xlab = "", ylab = "y",
         main = paste("Binomial data, n =", length(x$data$y)),
         ylim = c(0, max(x$data$y) + 1), pch = 16)
  }
}
