#' Plot estimated survival functions and cure probability for the promotion
#' time cure model.
#'
#' The routine takes as input an object of class \code{curelps} and plots
#' the estimated baseline survival curve, the population survival
#' curve for a specific covariate profile and a a smooth curve for the cure
#' probability. Approximate pointwise credible intervals are available.
#'
#' @usage
#' \method{plot}{curelps}(x, cred.int = 0.95, curvetype = c("baseline", "population",
#'      "probacure"), covar.profile, fit.col = "black", shade.col = "gray75",
#'       plot.cred = FALSE, np = 50, show.legend = TRUE, ...)
#'
#' @param x An object of class \code{curelps}.
#' @param cred.int The level for an approximate pointwise credible interval
#'   to be computed for the smooth curves. Default is 0.95.
#' @param curvetype The curve to be plotted; \code{baseline} (the default) will
#'   plot the estimated baseline survival, \code{population} the population
#'   survival function for a profile of covariates given in
#'   \code{covar.profile}, and \code{probacure} the probability to be cured (for
#'   a profile of covariates given in \code{covar.profile}) given
#'   that the subject has survived until time t.
#' @param covar.profile A numeric vector of the same length as the number
#'   of covariates in the model. This corresponds to the profile of covariates
#'   for which to compute the population survival function  and cure probability
#'   estimates. The order of the covariates in \code{covar.profile} is the same
#'   as the order specified in \code{formula} of the \code{curelps} routine.
#'   Each component of \code{covar.profile} should be in the range of the
#'   observed values for the corresponding covariate. If \code{covar.profile}
#'   is left unspecified by the user, the default will be to take the median
#'   covariate values.
#' @param fit.col The color used for the estimated survival curve.
#' @param shade.col The color used for the shading of the credible intervals.
#' @param plot.cred Logical. Should the credible intervals be plotted?
#'   Default is \code{FALSE}.
#' @param np The number of points used to plot the smooth curves. Default is 50
#'  and allowed values are between 20 and 200.
#' @param show.legend Show the legend? Default is TRUE.
#' @param ... Further arguments to be passed to \code{plot}.
#'
#' @details When \code{plot.cred} is \code{FALSE}, the routine omits to compute
#'  the approximate pointwise credible intervals for plotting and hence is
#'  less computationally intensive.
#'
#' @seealso \code{\link{curelps}}, \code{\link{curelps.object}},
#'  \code{\link{curelps.extract}}, \code{\link{print.curelps}}.
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @examples
#'
#' # Example on phase III clinical trial e1684 on melanoma data
#'
#' data(ecog1684)
#'
#' # Kaplan-Meier curve
#' plot(survfit(Surv(time, status) ~ 1, data = ecog1684), mark.time = TRUE)
#' fit <- curelps(Surv(time, status) ~ lt(age + trt + sex) +
#'              st(age + trt + sex), data = ecog1684, K = 20, penorder = 2)
#' fit
#' profile1 <- c(0, 1, 1, 0, 1, 1) # Mean age, trt = IFN, sex = Female.
#' profile2 <- c(0, 0, 1, 0, 0, 1) # Mean age, trt = control, sex = Female.
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 2))
#' plot(fit, curvetype = "probacure", plot.cred = TRUE, ylim = c(0,1),
#'      covar.profile = profile1, cred.int = 0.90,
#'      main = "Mean age, trt = IFN, sex = Female", cex.main = 0.8,
#'      show.legend = FALSE)
#' plot(fit, curvetype = "probacure", plot.cred = TRUE, ylim = c(0,1),
#'      covar.profile = profile2, cred.int = 0.90,
#'      main = "Mean age, trt = control, sex = Female", cex.main = 0.8,
#'      show.legend = FALSE)
#' par(opar)
#' @export

plot.curelps <- function(x, cred.int = 0.95,
                         curvetype = c("baseline", "population", "probacure"), covar.profile,
                         fit.col = "black", shade.col ="gray75",
                         plot.cred = FALSE, np = 50, show.legend = TRUE, ...){
  if (!inherits(x, "curelps"))
    stop("Object must be of class curelps")
  if (!is.vector(cred.int, mode = "numeric") ||
      length(cred.int) > 1 ||
      is.na(cred.int) || is.infinite(cred.int) ||
      cred.int <= 0 || cred.int >= 1)
    stop("cred.int must be between 0 and 1")
  if (!is.logical(plot.cred) || length(plot.cred) > 1)
    stop("plot.cred must be either TRUE or FALSE")
  if (np < 20 || np > 200)
    stop("np should be between 20 and 200")
  XZ <- cbind(x$X, x$Z)[,-1]
  ncovar <- ncol(XZ)
  if(!missing(covar.profile)){
    if(!is.vector(covar.profile, mode = "numeric") || length(covar.profile) != ncovar )
      stop("covar.profile must be a numeric vector with as much components as the
           number of covariates in the model formula")
    check.range <- c()
    for (j in 1:ncovar){
      check.range[j] <- (covar.profile[j] >= min(XZ[, j]) && covar.profile[j] <= max(XZ[, j]))
    }
    if(sum(check.range) != ncovar)
      stop("covariates in covar.profile are not in the range of observed covariate values")
  }


  tup <- x$tup
  K <- x$K
  thetahat <- x$spline.estim[-K]
  df.student <- x$n - x$ED
  plotcurve <- match.arg(curvetype)
  if (plotcurve == "baseline") {
    type <- "baseline"
  }
  else if (plotcurve == "population") {
    type <- "population"
  }
  else if (plotcurve == "probacure") {
    type <- "probacure"
  }
  else stop("Specify a correct curve type")

  # Fine grid of time values
  tfine <- seq(0, stats::quantile(x$event.times, prob = 0.99),
               length = np)

  # Baseline hazard function
  basehaz.fit <- function(t) {
    y <- exp(as.numeric(cubicbs(t, lower = 0, upper = tup,
                                K = (K - 1))$Bmatrix %*% thetahat))
    return(y)
  }

  # Baseline survival function
  basesurv.fit <- function(t) {
    y <- exp(-stats::integrate(basehaz.fit, lower = 0, upper = t)$value)
    return(y)
  }

  # Estimated baseline survival on tfine
  basesurv.estim <- sapply(tfine, basesurv.fit)

  # Propose a new grid of time values to truncate where baseline is nearly 0
  if (any(basesurv.estim < 1e-5)) {
    t.cut <- round(tfine[utils::head(which(basesurv.estim < 1e-5), 1)], 2)
    tfine <- seq(0, t.cut, length = np)
    basesurv.estim <- sapply(tfine, basesurv.fit)
  }

  # Credible interval at level cred.int for baseline survival at time t
  credS0 <- function(t) {
    if (t == 0) {
      CI.lb <- 1
      CI.ub <- 1
    } else {
      G0 <- function(t) {
        log(stats::integrate(basehaz.fit, lower = 0, upper = t)$value)
      }
      DG0 <- function(t) {
        factor <- ((-1) * log(basesurv.fit(t))) ^ (-1)
        bbashaz <- function(s, k) {
          cubs <- cubicbs(s, lower = 0, upper = tup, K = (K - 1))$Bmatrix
          y <- as.numeric(exp(cubs %*% thetahat)) * cubs[, k]
          return(y)
        }
        DG0.vec <- c()

        for (k in 1:(K - 1)) {
          DG0.vec[k] <- stats::integrate(bbashaz, k = k, lower = 0, upper = t)$value
        }
        val <- factor * DG0.vec
        return(val)
      }

      # Posterior mean at t
      postG0.mean <- G0(t)
      # Posterior standard deviation
      DG0t <- DG0(t)
      postG0.sd <- sqrt(sum((DG0t * x$Covtheta.map) %*% DG0t))
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

  # Credible intervals for baseline survival on tfine
  if (plot.cred == TRUE && type == "baseline") {
    CIS0.mat <- matrix(sapply(tfine, credS0), ncol = 2, byrow = TRUE)
    ylimsurv.lb <- min(CIS0.mat[, 1])
  }

  if(type == "baseline") {
    # Plot of baseline survival
    graphics::plot(tfine, basesurv.estim, type = "l", xlab = "t",
         ylab = expression(S[0](t)), col = fit.col, ...)
    if (plot.cred == TRUE) {
      graphics::polygon(x = c(tfine, rev(tfine)),
              y = c(CIS0.mat[, 1], rev(CIS0.mat[, 2])),
              col = shade.col, border = NA)
    }
    graphics::lines(tfine, basesurv.estim, col = fit.col, lwd = 2)
    if(show.legend == TRUE){
      if (plot.cred == TRUE) {
        graphics::legend("topright", c("Estimated baseline survival",
                             paste0((round(cred.int * 100 , 2)), "% CI")),
               col = c(fit.col, shade.col), lwd = c(2, 2), bty = "n")
      } else{
        graphics::legend("topright", c("Estimated baseline survival"), col = fit.col,
               bty = "n", lwd = 2)
      }
    }
  }

  # Population survival function at median covariate values
  if (missing(covar.profile)) {
    covar.profile <- c(1, apply(XZ, 2, stats::median))
  } else{
    covar.profile <- c(1, covar.profile)
  }
  beta.hat <- as.numeric(x$coeff.probacure[, 1])
  gamma.hat <- as.numeric(x$coeff.cox[, 1])
  Xbeta <- as.numeric(covar.profile[1:ncol(x$X)] %*% beta.hat)
  Zgamma <- as.numeric(covar.profile[(ncol(x$X)+1):(ncol(XZ) + 1)] %*% gamma.hat)
  popsurv.fit <- function(t) {
    y <- exp(-exp(Xbeta) * (1 - ((basesurv.fit(t)) ^ (exp(Zgamma)))))
    return(y)
  }
  popsurv.estim <- sapply(tfine, popsurv.fit)

  # Credible interval for population survival function
  credS0.pop <- function(t) {
    if (t == 0) {
      CI.lb <- 1
      CI.ub <- 1
    } else {
      G0.pop <- function(t) {
        Xbeta + log(1 - ((basesurv.fit(t)) ^ (exp(Zgamma))))
      }
      DG0.pop <- function(t) {
        S0t <- basesurv.fit(t)
        if (S0t == 0) S0t <- 1e-10
        v.inv <- (1 - (S0t ^ (exp(Zgamma)))) ^ (-1)
        factor1 <- v.inv * exp(Zgamma) * (S0t ^ (exp(Zgamma)))
        factor2 <- (-1) * log(S0t)
        bbashaz <- function(s, k) {
          cubs <- cubicbs(s, lower = 0, upper = tup, K = (K - 1))$Bmatrix
          y <- as.numeric(exp(cubs %*% thetahat)) * cubs[, k]
          return(y)
        }
        DG0.vec <- c()
        for (k in 1:(K - 1)) {
          DG0.vec[k] <- stats::integrate(bbashaz, k = k, lower = 0, upper = t)$value
        }
        DG0.vec <- DG0.vec
        DG0.popval <- c(factor1 * DG0.vec,
                        covar.profile[1:ncol(x$X)],
                        factor1 * factor2 * covar.profile[(ncol(x$X)+1):(ncol(XZ) + 1)])
        return(DG0.popval)
      }
      # Posterior mean at t
      postG0.mean <- G0.pop(t)
      # Posterior standard deviation
      DG0t <- DG0.pop(t)
      postG0.sd <- sqrt(sum((DG0t * x$Covlatc.map) %*% DG0t))
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

  # Credible intervals for population survival
  if (plot.cred == TRUE && type == "population") {
    CIS0pop.mat <- matrix(sapply(tfine, credS0.pop), ncol = 2, byrow = TRUE)
    ylimsurvpop.lb <- min(CIS0pop.mat[, 1])
  }

  if(type == "population") {
    # Plot population survival
    graphics::plot(tfine, popsurv.estim, type = "l", xlab = "t",
         ylab = expression(S[pop](t)), col = fit.col, ...)
    if (plot.cred == TRUE) {
      graphics::polygon(x = c(tfine, rev(tfine)),
              y = c(CIS0pop.mat[, 1], rev(CIS0pop.mat[, 2])),
              col = shade.col,
              border = NA)
    }
    graphics::lines(tfine, popsurv.estim, col = fit.col, lwd = 2)
    if(show.legend == TRUE){
      if (plot.cred == TRUE) {
        graphics::legend("topright", c("Population survival",
                             paste0((round(cred.int * 100 , 2)), "% CI")),
               col = c(fit.col, shade.col), lwd = c(2, 2), bty = "n")
      } else{
        graphics::legend("topright", c("Population survival"), col = fit.col, bty = "n",
               lwd = 2)
      }
    }
  }

  # Probability to be cured at t for a given profile of covariates
  probcure <- function(t) {
    y <- exp(-exp(Xbeta) * ((basesurv.fit(t)) ^ (exp(Zgamma))))
    return(y)
  }

  # Credible interval at time t for probability to be cured
  credprobcure <- function(t) {
    G0.probcure <- function(t) {
      S0t <- basesurv.fit(t)
      if (S0t == 0) S0t <- 1e-15
      Xbeta + exp(Zgamma) * log(S0t)
    }
    DG0.probcure <- function(t) {
      S0t <- basesurv.fit(t)
      if (S0t == 0) S0t <- 1e-15
      bbashaz <- function(s, k) {
        cubs <- cubicbs(s, lower = 0, upper = tup, K = (K - 1))$Bmatrix
        y <- as.numeric(exp(cubs %*% thetahat)) * cubs[, k]
        return(y)
      }
      DG0.vec <- c()
      for (k in 1: (K - 1)) {
        DG0.vec[k] <- stats::integrate(bbashaz, k = k, lower = 0, upper = t)$value
      }
      DG0.vec <- DG0.vec
      DG0val.probcureval <- c(-exp(Zgamma) * DG0.vec,
                              covar.profile[1:ncol(x$X)],
                              covar.profile[(ncol(x$X)+1):(ncol(XZ) + 1)] * exp(Zgamma)  * log(S0t))
      return(DG0val.probcureval)
    }

    # Posterior mean at t
    postG0.mean <- G0.probcure(t)

    # Posterior standard deviation
    DG0t <- DG0.probcure(t)
    postG0.sd <- sqrt(sum((DG0t * x$Covlatc.map) %*% DG0t))
    # Approximate credible interval
    # z.quant <- qnorm(.5 * (1 - cred.int), lower.tail = FALSE) # Normal quantile
    t.quant <- stats::qt(.5 * (1 - cred.int), df = df.student,
                         lower.tail = FALSE) # Student quantile
    CI.lb <- exp(-exp(postG0.mean + t.quant *
                        sqrt(df.student/(df.student - 2)) * postG0.sd))
    CI.ub <- exp(-exp(postG0.mean - t.quant *
                        sqrt(df.student/(df.student - 2)) * postG0.sd))
    return(c(CI.lb, CI.ub))
  }

  if(plot.cred == TRUE) {
    CIprobcure.mat <- matrix(sapply(tfine, credprobcure), ncol = 2, byrow = T)
  }

  if(type == "probacure"){
    # Plot probability to be cured over time
    probcure.estim <- sapply(tfine, probcure)

    graphics::plot(tfine, probcure.estim, col = fit.col, type = "l", xlab = "t",
         ylab = "Proba. to be cured", ...)
    if (plot.cred == TRUE) {
      graphics::polygon(x = c(tfine, rev(tfine)),
              y = c(CIprobcure.mat[, 1], rev(CIprobcure.mat[, 2])),
              col = shade.col,
              border = NA)
    }
    graphics::lines(tfine, probcure.estim, col = fit.col, lwd = 2)
    if(show.legend == TRUE){
      if (plot.cred == TRUE) {
        graphics::legend("topleft", c("Cure probability", paste0((round(cred.int * 100 , 2)),
                                                       "% CI")),
               col = c(fit.col, shade.col), lwd = c(2, 2), bty = "n")
      } else{
        graphics::legend("topleft", c("Cure probability"), col = fit.col, bty = "n", lwd = 2)
      }
    }
  }

}
