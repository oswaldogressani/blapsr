#' @method plot simsurvdata
#' @export

# Plot method for simsurvdata class
plot.simsurvdata <- function(x, ...) {
  n <- x$sample.size
  if (n <= 25) obs.plot <- x$sample.size else obs.plot <- 25
  yindex <- matrix(rep(seq_len(obs.plot), 2), ncol = 2, byrow = FALSE)
  tindex <- matrix(c(rep(0, obs.plot), x$survdata$time[1:obs.plot]),
                   ncol = 2, byrow = FALSE)
  graphics::plot(tindex[1, ], yindex[1, ], type ="l",
       ylim = c(0, obs.plot + 10),
       xlim = c(0, max(x$survdata$time[1:obs.plot])), ylab = "",
       xlab = "time")
  if(x$survdata$delta[1:obs.plot][1] == 1) {
    graphics::lines(tindex[1, 2], yindex[1,2], type ="p", pch = 19)
  } else {
    graphics::lines(tindex[1, 2], yindex[1,2], type ="p", pch = 4, lwd = 2)
  }
  if (obs.plot > 1) {
    for (i in 2:obs.plot) {
      graphics::lines(tindex[i, ], yindex[i, ], type = "l")
      if(x$survdata$delta[1:obs.plot][i] == 1) {
        graphics::lines(tindex[i, 2], yindex[i, 2], type ="p", pch = 19)
      } else {
        graphics::lines(tindex[i, 2], yindex[i,2], type ="p", pch = 4, lwd = 2)
      }
    }
  }
  graphics::legend("topright", c("Right censored", "Event occured"),
         lty = c("blank", "blank"), pch = c(4, 19), lwd = c(2, 1), bty = "n")
}

