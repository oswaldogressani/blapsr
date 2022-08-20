#' Plot baseline hazard and survival curves from a coxlps object.
#'
#' @description Produces a plot of the baseline hazard and/or survival based
#'  on a coxlps object.
#'
#'
#' @usage
#' \method{plot}{coxlps}(x, S0 = TRUE, h0 = TRUE, cred.int = 0.95, overlay.km = FALSE,
#'      plot.cred = FALSE, np = 50, show.legend = TRUE, ...)
#'
#' @param x An object of class \code{coxlps}.
#' @param S0 Logical. Should the estimated baseline survival be plotted?
#' @param h0 Logical. Should the estimated baseline hazard be plotted?
#' @param cred.int The level for an approximate pointwise credible interval
#'  to be computed for the baseline curves. Default is 0.95.
#' @param overlay.km  A logical value indicating whether the Kaplan-Meier
#'   estimate should be plotted together with the smooth baseline survival
#'   curve. The default is \code{FALSE}.
#' @param plot.cred Logical. Should the credible intervals be plotted ?
#'   Default is \code{FALSE}.
#' @param np The number of points used to plot the smooth baseline
#'   functions. Default is 50 and allowed values are
#'   between 20 and 200.
#' @param show.legend Logical. Should a legend be displayed?
#' @param ... Further arguments to be passed to plot.
#'
#' @details Plots for the baseline hazard and survival curves are computed on
#'   a grid (of length \code{np}) between 0 and the 99th percentile of follow-up
#'   times. When \code{plot.cred} is \code{FALSE}, the fit omits to compute the
#'   approximate pointwise credible intervals for plotting and hence is less
#'   computationally intensive. Vertical ticks on the x-axis correspond to the
#'   observed follow-up times.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @examples
#'
#'
#' ## Simulate survival data
#' set.seed(6)
#' betas <- c(0.35, -0.20, 0.05, 0.80) # Regression coefficients
#' data <- simsurvdata(a = 1.8, b = 2, n = 200, betas = betas, censperc = 25)
#' simdat <- data$survdata
#'
#' # Fit model
#' fit <- coxlps(Surv(time, delta) ~ x1 + x2 + x3 + x4, data = simdat)
#' plot(fit, h0 = FALSE, S0 = TRUE, overlay.km = FALSE, show.legend = FALSE)
#' domt <- seq(0, 5.5, length = 500)
#' lines(domt, data$S0(domt), type = "l", col = "red")
#' legend("topright", c("Bayesian LPS", "Target"), col = c("black", "red"),
#'        lty = c(1, 1), bty = "n", cex = 0.8)
#'
#' @seealso \code{\link{coxlps}} \code{\link{coxlps.object}}
#'
#' @export

plot.coxlps<- function(x, S0 = TRUE, h0 = TRUE, cred.int = 0.95,
          overlay.km = FALSE, plot.cred = FALSE, np = 50, show.legend = TRUE,
          ...){

  if (!inherits(x, "coxlps"))
    stop("Object must be of class coxlps")
  if (!is.vector(cred.int, mode = "numeric") ||
      length(cred.int) > 1 ||
      is.na(cred.int) || is.infinite(cred.int) ||
      cred.int <= 0 || cred.int >= 1)
    stop("cred.int must be between 0 and 1")
  if (!is.logical(overlay.km) || length(overlay.km) > 1)
    stop("overlay.km must be either TRUE or FALSE")
  if (!is.logical(plot.cred) || length(plot.cred) > 1)
    stop("plot.cred must be either TRUE or FALSE")
  if (np < 20 || np > 200)
    stop("Specify n between 20 and 200")

  # Fine grid of time values
  q99 <- stats::quantile(x$event.times, prob = 0.99)
  tfine <- seq(0, q99, length = np)
  ttobs <- x$event.times[x$event.times <= q99] * x$sd.time
  df.student <- x$n - x$ED

  # Fitted baseline hazard as a function of t
  basehaz.fit <- function(t) {
    y <- as.numeric(exp(crossprod(t(cubicbs(t, lower = 0, upper = x$tup,
                                K = x$K)$Bmatrix), x$spline.estim)))
    return(y)
  }

  # Fitted baseline survival as a function of t
  basesurv.fit <- function(t) {
    y <- exp(-stats::integrate(basehaz.fit, lower = 0, upper = t)$value)
    return(y)
  }

  # Credible interval at level cred.int for baseline survival at time t
  credS0 <- function(t) {
    if (t == 0) {
      CI.lb <- 1
      CI.ub <- 1
    } else{
      G0 <-function(t) log(stats::integrate(basehaz.fit, lower = 0, upper = t)$value)
      DG0 <- function(t) {
        bbashaz <- function(s, k) {
          cubs <- cubicbs(s, lower = 0, upper = x$tup,
                          K = x$K)$Bmatrix
          y <- as.numeric(exp(cubs %*% x$spline.estim)) * cubs[, k]
          return(y)
        }

        DG0.vec <- c()

        for (k in 1:x$K) {
          DG0.vec[k] <- stats::integrate(bbashaz, k = k, lower = 0, upper = t)$value
        }

        val <- DG0.vec / (stats::integrate(basehaz.fit, lower = 0, upper = t)$value)
        return(val)
      }

      # Posterior mean
      postG0.mean <- G0(t)
      # Posterior standard deviation
      postG0.sd <- sqrt(sum((DG0(t) * x$Covthetamix) %*% DG0(t)))
      # Approximate credible interval
      # z.quant <- qnorm(.5 * (1 - cred.int), lower.tail = FALSE) # Normal quantile
      t.quant <- stats::qt(.5 * (1 - cred.int), df = df.student,
                           lower.tail = FALSE) # Student quantile
      CI.lb <- exp(-exp(postG0.mean + t.quant *
                          sqrt(df.student/(df.student - 2)) * postG0.sd))
      CI.ub <- exp(-exp(postG0.mean - t.quant *
                          sqrt(df.student/(df.student - 2)) * postG0.sd))
    }
    return(c(CI.lb, CI.ub))
  }

  # Credible interval at level cred.int for baseline hazard at time t
  credh0 <- function(t) {
    Bspline.eval <- as.numeric(cubicbs(t, lower = 0, upper = x$tup,
                                       K = x$K)$Bmatrix)
    # Posterior mean
    postG0.mean <- as.numeric(crossprod(Bspline.eval, x$spline.estim))
    # Posterior standard deviation
    postG0.sd <- sqrt(sum((Bspline.eval * x$Covthetamix) %*%
                            Bspline.eval))
    # Approximate credible interval
    z.quant <- stats::qnorm(.5 * (1 - cred.int), lower.tail = FALSE) # Normal quantile
    CI.lb <- exp(postG0.mean - z.quant * postG0.sd)
    CI.ub <- exp(postG0.mean + z.quant * postG0.sd)
    return(c(CI.lb, CI.ub))
  }

  # Fitted baseline hazard and survival on tfine
  basehaz.estim <- sapply(tfine, basehaz.fit)
  basesurv.estim <- sapply(tfine, basesurv.fit)

  if (plot.cred == TRUE) {
    # Credible intervals for baseline survival on tfine
    CIS0.mat <- matrix(sapply(tfine, credS0), ncol = 2, byrow = TRUE)
    ylimsurv.lb <- min(CIS0.mat[, 1])
    # Credible intervals for baseline hazard on tfine
    CIh0.mat <- matrix(sapply(tfine, credh0), ncol = 2, byrow = TRUE)
    ylimhaz.lb <- min(CIh0.mat[, 1])
    ylimhaz.ub <- max(CIh0.mat[, 2])
  } else{
    ylimsurv.lb <- min(basesurv.estim)
    ylimhaz.lb <- min(basehaz.estim)
    ylimhaz.ub <- max(basehaz.estim)
  }

  # Fit with coxph function
  n.regcoeff <- x$latfield.dim - x$K
  dataframe <- data.frame(scale(x$X, center = TRUE, scale = FALSE))
  colnames(dataframe) <- paste0(rep("x", n.regcoeff), seq_len(n.regcoeff))
  newdata <- data.frame(matrix(0 , ncol = n.regcoeff))
  colnames(newdata) <- colnames(dataframe)
  formula.coxph <- stats::as.formula(paste(
    "Surv((x$event.times * x$sd.time) ,x$event.indicators)~",
    paste(colnames(dataframe), collapse = "+")))
  fit.coxph <- survival::coxph(formula.coxph, data = dataframe)
  tnatscale <- tfine * x$sd.time # Time in natural scale

  # Plot baseline estimations
  if(h0 == TRUE && S0 == TRUE){
    opar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))
    graphics::par(mfrow = c(1,2))
  }
  if(h0 == TRUE){
    # Baseline hazard
    graphics::plot(tnatscale, basehaz.estim, type = "l", col = "black",
         ylab = expression(h[0](t)), xlab = "t",
         ylim = c(ylimhaz.lb, ylimhaz.ub), ...)
    if(show.legend == TRUE){
      graphics::legend("topleft", "Bayesian LPS", lty = 1, col = "black",
             cex = .8, bty = "n")
    }
    graphics::rug(ttobs)
    if(plot.cred == TRUE){
      graphics::polygon(x = c(tnatscale, rev(tnatscale)),
              y = c(CIh0.mat[,1], rev(CIh0.mat[,2])),
              col = "gray", border = NA)
      graphics::lines(tnatscale, basehaz.estim, type = "l", col = "black", lwd = 2)
      if(show.legend == TRUE){
        graphics::legend("topleft", c("Bayesian LPS",
              paste(cred.int * 100, "% CI", sep = "")), lty = c(1, 1),
               col=c("black","gray"), cex=0.8, bty="n")}
    }
  }

  if(S0 == TRUE){
    # Baseline survival
    if(plot.cred == FALSE && overlay.km == FALSE){
      graphics::plot(tnatscale, basesurv.estim, type = "l", col = "black",
           ylab = expression(S[0](t)), xlab = "t", ylim = c(ylimsurv.lb, 1),
           ...)
      if(show.legend == TRUE){
        graphics::legend("topright", "Bayesian LPS", lty = 1, col = "black",
               cex = .8, bty = "n")
      }
      graphics::rug(ttobs)
    }

    if(plot.cred == FALSE && overlay.km == TRUE){
      graphics::plot(tnatscale, basesurv.estim, type = "l", col = "black",
           ylab = expression(hat(S)[0](t)), xlab = "t", ylim=c(ylimsurv.lb ,1),
           ...)
      graphics::lines(survival::survfit(fit.coxph, newdata = newdata),
                      col = "blue", conf.int = FALSE)
      if(show.legend == TRUE){
        graphics::legend("topright", c("Bayesian LPS", "Kaplan-Meier"), lty = c(1, 1),
               col = c("black", "blue"), cex = .8, bty = "n")
      }
      graphics::rug(ttobs)
    }

    if(plot.cred == TRUE && overlay.km == FALSE){
      graphics::plot(tnatscale, basesurv.estim, type = "l", col = "black",
           ylab = expression(hat(S)[0](t)), xlab = "t", ylim = c(ylimsurv.lb, 1),
           ...)
      graphics::polygon(x = c(tnatscale, rev(tnatscale)),
              y = c(CIS0.mat[, 1], rev(CIS0.mat[, 2])),
              col = "gray", border=NA)
      graphics::lines(tnatscale, basesurv.estim, type = "l", col = "black", lwd = 2)
      if(show.legend == TRUE){
        graphics::legend("topright", c("Bayesian LPS",
                             paste(cred.int * 100, "% CI", sep = "")),
               lty = c(1, 1), col=c("black","gray"), cex=0.8, bty="n")
      }
      graphics::rug(ttobs)
    }

    if(plot.cred == TRUE && overlay.km == TRUE){
      graphics::plot(tnatscale, basesurv.estim, type = "l", col = "black",
           ylab = expression(hat(S)[0](t)), xlab = "t", ylim = c(ylimsurv.lb, 1),
           ...)
      graphics::polygon(x = c(tnatscale, rev(tnatscale)),
              y = c(CIS0.mat[, 1], rev(CIS0.mat[, 2])),
              col = "gray", border = NA)
      graphics::lines(tnatscale, basesurv.estim, type = "l", col = "black", lwd = 2)
      graphics::lines(survival::survfit(fit.coxph, newdata = newdata, conf.type="log-log",
                    conf.int = cred.int), col = "blue")
      if(show.legend == TRUE){
        graphics::legend("topright", c("Bayesian LPS",
                    paste("Kaplan-Meier+(", cred.int * 100, "% CI)", sep = ""),
                    paste(cred.int * 100, "% CI", sep = "")), lty = c(1, 1, 1),
               col = c("black", "blue", "gray"),
               cex = 0.8, bty = "n")
      }
      graphics::rug(ttobs)
    }
  }

}
