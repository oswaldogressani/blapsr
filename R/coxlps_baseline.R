#' Extract estimated baseline quantities from a fit with coxlps.
#'
#' The routine takes as input an object of class \code{coxlps} and computes
#' point estimates and credible intervals for the baseline hazard and survival
#' on a user-specified time vector.
#'
#' @usage coxlps.baseline(object, time = NULL, compute.cred = TRUE, cred.int = 0.95,
#'                 verbose = TRUE)
#'
#' @param object An object of class \code{coxlps}.
#' @param time A vector of time values on which to compute the estimated
#'   baseline quantities. Each component of \code{time} must be between 0 and
#'   the largest observed follow-up time. If time is \code{NULL} (the default),
#'   then only the baseline median lifetime (if available) is computed.
#' @param compute.cred Should the credible intervals be computed? Default is
#'  TRUE.
#' @param cred.int The level for an approximate pointwise credible interval
#'   to be computed for the baseline hazard and survival curves. Default
#'   is 0.95.
#' @param verbose Should the table of estimated values be printed to console?
#'   Default is TRUE.
#'
#' @return A list with the following components:
#'
#'   \item{fit.time}{A matrix with point and set estimates of
#'   the baseline hazard and survival curves for values provided in \code{time}.
#'   Only available if \code{time} is not \code{NULL}. Column \emph{Time}
#'   summarizes the provided values in \code{time}. Columns named \emph{h0},
#'   \emph{S0}, are the point estimates of the baseline hazard and baseline
#'   survival respectively. \emph{low} and \emph{up} give the lower and
#'   upper bound respectively of the approximate pointwise credible interval.}
#'
#'   \item{median.lifetime \verb{ }}{The estimated baseline median lifetime.}
#'
#'   \item{cred.int \verb{ }}{The chosen level to construct credible intervals.}
#'
#' @seealso \code{\link{coxlps}}, \code{\link{coxlps.object}}
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @examples
#'
#' ## Simulate survival data
#' set.seed(2)
#' betas <- c(0.15, 0.82, 0.41) # Regression coefficients
#' data <- simsurvdata(a = 1.8, b = 2, n = 300, betas = betas, censperc = 15)
#' simdat <- data$survdata
#'
#' # Fit model
#' fit <- coxlps(Surv(time, delta) ~ x1 + x2 + x3, data = simdat, K = 20)
#' coxlps.baseline(fit, time = seq(0, 2, by = 0.5), cred.int = 0.90)
#'
#' @export

coxlps.baseline <- function(object, time = NULL, compute.cred = TRUE, cred.int = 0.95,
                            verbose = TRUE){

  if (!inherits(object, "coxlps"))
    stop("Object must be of class coxlps")
  if (!is.null(time)) {
    if (any(is.na(time)) || any(is.infinite(time)) ||
        !is.vector(time, mode = "numeric") || any(time < 0))
      stop("time must be a finite numeric vector and cannot contain
           negative values")
    if (any((time / object$sd.time) > object$tup))
      stop("time must be between 0 and largest observed event time")
  }
  if (!is.vector(cred.int, mode = "numeric") ||
      length(cred.int) > 1 ||
      is.na(cred.int) || is.infinite(cred.int) ||
      cred.int <= 0 || cred.int >= 1)
    stop("cred.int must be between 0 and 1")
  df.student <- object$n - object$ED

  # Fitted baseline hazard as a function of t
  basehaz.fit <- function(t) {
    y <- as.numeric(exp(crossprod(t(cubicbs(t, lower = 0, upper = object$tup,
                    K = object$K)$Bmatrix), object$spline.estim)))
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
          cubs <- cubicbs(s, lower = 0, upper = object$tup,
                          K = object$K)$Bmatrix
          y <- as.numeric(exp(cubs %*% object$spline.estim)) * cubs[, k]
          return(y)
        }

        DG0.vec <- c()

        for (k in 1:object$K) {
          DG0.vec[k] <- stats::integrate(bbashaz, k = k, lower = 0, upper = t)$value
        }

        val <- DG0.vec / (stats::integrate(basehaz.fit, lower = 0, upper = t)$value)
        return(val)
      }

      # Posterior mean
      postG0.mean <- G0(t)
      # Posterior standard deviation
      postG0.sd <- sqrt(sum((DG0(t) * object$Covthetamix) %*% DG0(t)))
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
    Bspline.eval <- as.numeric(cubicbs(t, lower = 0, upper = object$tup,
                                       K = object$K)$Bmatrix)
    # Posterior mean
    postG0.mean <- as.numeric(crossprod(Bspline.eval, object$spline.estim))
    # Posterior standard deviation
    postG0.sd <- sqrt(sum((Bspline.eval * object$Covthetamix) %*%
                            Bspline.eval))
    # Approximate credible interval
    z.quant <- stats::qnorm(.5 * (1 - cred.int), lower.tail = FALSE) # Normal quantile
    CI.lb <- exp(postG0.mean - z.quant * postG0.sd)
    CI.ub <- exp(postG0.mean + z.quant * postG0.sd)
    return(c(CI.lb, CI.ub))
  }

  # Median baseline lifetime
  if (basesurv.fit(object$tup) > 0.5) {
    median.lifetime <- NA
  } else{
    median.baseline <- function(t) basesurv.fit(t) - .5
    median.lifetime <- stats::uniroot(median.baseline, lower = 0, upper = object$tup,
                               tol = 1e-5)$root
    median.lifetime <- median.lifetime * object$sd.time
    median.lifetime <- round(median.lifetime, 4)
  }

  if (!is.null(time)) {
    if (compute.cred == TRUE) {
      n.time <- length(time)
      fit.time <- matrix(0, nrow = n.time, ncol = 7)
      colnames(fit.time) <- c("Time", "h0", "h0.low", "h0.up", "S0", "S0.low",
                              "S0.up")
      credintsh0 <- t(sapply((time / object$sd.time), credh0))
      credintsS0 <- t(sapply((time / object$sd.time), credS0))
      fit.time[, 1] <- time
      fit.time[, 2] <- sapply((time / object$sd.time), basehaz.fit)
      fit.time[, 3] <- credintsh0[, 1]
      fit.time[, 4] <- credintsh0[, 2]
      fit.time[, 5] <- sapply((time / object$sd.time), basesurv.fit)
      fit.time[, 6] <- credintsS0[, 1]
      fit.time[, 7] <- credintsS0[, 2]
      fit.time <- round(fit.time, 4)
      outputlist <- list(fit.time = fit.time, median.lifetime = median.lifetime,
                         cred.int = cred.int)
      if (verbose == TRUE) {
        cat("Estimated baseline hazard and survival at specified time points (*): \n")
        cat("\n")
        print.table(format(fit.time, nsmall = 4), right = TRUE)
        cat("--- \n")
        cat("* Bounds correspond to a",
            paste(format(round(cred.int  * 100, 2), nsmall = 2), "%", sep = ""),
            "credible interval. \n")
        cat("\n")
        cat("Median lifetime:", format(median.lifetime, nsmall = 4))
      }
    } else{
      n.time <- length(time)
      fit.time <- matrix(0, nrow = n.time, ncol = 3)
      colnames(fit.time) <- c("Time", "h0", "S0")
      fit.time[, 1] <- time
      fit.time[, 2] <- sapply((time / object$sd.time), basehaz.fit)
      fit.time[, 3] <- sapply((time / object$sd.time), basesurv.fit)
      outputlist <- list(fit.time = fit.time, median.lifetime = median.lifetime,
                         cred.int = cred.int)
      if (verbose == TRUE) {
        cat("Estimated baseline hazard and survival at specified time points : \n")
        cat("\n")
        print.table(format(fit.time, nsmall = 4), right = TRUE)
        cat("--- \n")
        cat("\n")
        cat("Median lifetime:", format(median.lifetime, nsmall = 4))
      }
     }
    }else{outputlist <- list(median.lifetime = median.lifetime,
                            cred.int = cred.int)
  }
  return(invisible(outputlist))
}

