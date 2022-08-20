#' Extract estimates of survival functions and cure probability for the
#' promotion time cure model.
#'
#' The routine takes as input an object of class \code{curelps} and computes
#' estimates of the baseline survival curve, the population survival
#' curve and the cure probability on a specified time vector. Approximate
#' pointwise credible intervals are available.
#'
#' @usage curelps.extract(object, time = NULL, curvetype = c("baseline", "population", "probacure"),
#'                 covar.profile, compute.cred = TRUE, cred.int = 0.95, verbose = TRUE)
#'
#' @param object An object of class \code{curelps}.
#' @param time A vector of time values on which to compute the estimates.
#'   Each component of \code{time} must be between 0 and the largest observed
#'   follow-up time.
#' @param curvetype The curve on which estimates are computed ; \code{baseline}
#'   (the default) is for the baseline survival, \code{population} is for the
#'   population survival function for a profile of covariates given in
#'   \code{covar.profile}, and \code{probacure} is for the probability to be
#'   cured (for a profile of covariates given in \code{covar.profile}) given
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
#' @param compute.cred Should credible intervals be computed? Default is TRUE.
#' @param cred.int The level for an approximate pointwise credible interval.
#'   Default is 0.95.
#' @param verbose Should estimates be printed to console?
#'
#' @return A list with the following components:
#'
#' \item{fit.time}{Estimates on the time values provided in \code{time}.}
#'
#' \item{cred.int}{The chosen level to construct approximate pointwise
#'   credible intervals.}
#'
#' \item{covar.profile}{The chosen profile of covariates.}
#'
#' @seealso \code{\link{curelps}}, \code{\link{curelps.object}},
#'  \code{\link{plot.curelps}}, \code{\link{print.curelps}}.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @examples
#'
#' # Example on phase III clinical trial e1684 on melanoma data
#'
#' data(ecog1684)
#'
#' # Kaplan-Meier curve
#' plot(survfit(Surv(time, status) ~ 1, data = ecog1684), mark.time = TRUE)
#' fit <- curelps(Surv(time, status) ~ lt(age + trt+ sex) +
#'              st(age + trt + sex), data = ecog1684, K = 20, penorder = 2)
#' fit
#' profile1 <- c(0, 1, 1, 0, 1, 1) # Mean age, trt = IFN, sex = Female.
#' profile2 <- c(0, 0, 1, 0, 0, 1) # Mean age, trt = control, sex = Female.
#'
#' # Extract cure probabilities
#' curelps.extract(fit, time = c(0, 1, 2, 3), curvetype = "probacure",
#'                 covar.profile = profile1, cred.int = 0.90)
#' curelps.extract(fit, time = c(0, 1, 2, 3), curvetype = "probacure",
#'                 covar.profile = profile2, cred.int = 0.90)
#' @export

curelps.extract <- function(object, time = NULL,
                  curvetype = c("baseline", "population", "probacure"),
      covar.profile, compute.cred = TRUE, cred.int = 0.95, verbose = TRUE){


  if (!inherits(object, "curelps"))
    stop("Object must be of class curelps")
  if (!is.null(time)) {
    if (any(is.na(time)) || any(is.infinite(time)) ||
        !is.vector(time, mode = "numeric") || any(time < 0))
      stop("time must be a finite numeric vector and cannot contain negative,
           NA/NaN values")
    if(any(time > object$tup))
      stop("time must be between 0 and the largest observed event time")
  }
  if (!is.vector(cred.int, mode = "numeric") ||
      length(cred.int) > 1 ||
      is.na(cred.int) || is.infinite(cred.int) ||
      cred.int <= 0 || cred.int >= 1)
    stop("cred.int must be between 0 and 1")
  XZ <- cbind(object$X, object$Z)[,-1]
  ncovar <- ncol(XZ)
  if(!missing(covar.profile)){
    if(!is.vector(covar.profile, mode = "numeric") || length(covar.profile) != ncovar)
      stop("covar.profile must be a numeric vector with as much components as the
           number of covariates in the model formula")
    check.range <- c()
    for (j in 1:ncovar){
      check.range[j] <- (covar.profile[j] >= min(XZ[, j]) &&
                           covar.profile[j] <= max(XZ[, j]))
    }
    if(sum(check.range) != ncovar)
      stop("covariates in covar.profile are not in the range of observed covariate values")
  }

  tup <- object$tup
  K <- object$K
  thetahat <- object$spline.estim[-K]
  df.student <- object$n - object$ED
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
      postG0.sd <- sqrt(sum((DG0t * object$Covtheta.map) %*% DG0t))
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

  # Population survival function at median covariate values
  if (missing(covar.profile)) {
    covar.profile <- c(1, apply(XZ, 2, stats::median))
  } else{
    covar.profile <- c(1, covar.profile)
  }
  beta.hat <- as.numeric(object$coeff.probacure[, 1])
  gamma.hat <- as.numeric(object$coeff.cox[, 1])
  Xbeta <- as.numeric(covar.profile[1:ncol(object$X)] %*% beta.hat)
  Zgamma <- as.numeric(covar.profile[(ncol(object$X)+1):(ncol(XZ) + 1)] %*% gamma.hat)

  popsurv.fit <- function(t) {
    y <- exp(-exp(Xbeta) * (1 - ((basesurv.fit(t)) ^ (exp(Zgamma)))))
    return(y)
  }

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
                        covar.profile[1:ncol(object$X)],
                        factor1 * factor2 * covar.profile[(ncol(object$X)+1):(ncol(XZ) + 1)])
        return(DG0.popval)
      }
      # Posterior mean at t
      postG0.mean <- G0.pop(t)
      # Posterior standard deviation
      DG0t <- DG0.pop(t)
      postG0.sd <- sqrt(sum((DG0t * object$Covlatc.map) %*% DG0t))
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
                              covar.profile[1:ncol(object$X)],
                              covar.profile[(ncol(object$X)+1):(ncol(XZ) + 1)] * exp(Zgamma)  * log(S0t))
      return(DG0val.probcureval)
    }

    # Posterior mean at t
    postG0.mean <- G0.probcure(t)

    # Posterior standard deviation
    DG0t <- DG0.probcure(t)
    postG0.sd <- sqrt(sum((DG0t * object$Covlatc.map) %*% DG0t))
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

  if (!is.null(time)) {
    if(type == "baseline") {
      # Output fitted results on supplied time vector
      n.time <- length(time)
      if(compute.cred == TRUE) ncols = 4 else ncols = 2
      fit.time <- matrix(0, nrow = n.time, ncol = ncols)
      if(compute.cred == TRUE){
        colnames(fit.time) <- c("Time", "S0", "S0.low", "S0.up")
      } else colnames(fit.time) <- c("Time", "S0")

      fit.time[, 1] <- time
      fit.time[, 2] <- sapply(time, basesurv.fit)

      if(compute.cred == TRUE){
        credintS0 <- t(sapply(time, credS0))
        fit.time[, 3] <- credintS0[, 1]
        fit.time[, 4] <- credintS0[, 2]
      }
      fit.time <- round(fit.time, 4)
      outputlist <- list(fit.time = fit.time, cred.int = cred.int)
      if(verbose == TRUE) {
        if(compute.cred == TRUE){
          cat("Estimated baseline survival at specified time points (*): \n")
          cat("\n")
          print.table(format(fit.time, nsmall = 4), right = TRUE)
          cat("--- \n")
          cat("* Bounds correspond to a", paste(format(round(cred.int  * 100, 2),
                                                       nsmall = 2), "%", sep = ""), "credible interval. \n")
        } else{
          cat("Estimated baseline survival at specified time points: \n")
          cat("\n")
          print.table(format(fit.time, nsmall = 4), right = TRUE)
          cat("--- \n")
        }
      }
      return(invisible(outputlist))
    } else if (type == "population") {
      # Output fitted results on supplied time vector
      n.time <- length(time)
      if(compute.cred == TRUE) ncols = 4 else ncols = 2
      fit.time <- matrix(0, nrow = n.time, ncol = ncols)
      if(compute.cred == TRUE){
        colnames(fit.time) <- c("Time", "Spop(**)", "Spop.low", "Spop.up")
      } else colnames(fit.time) <- c("Time", "Spop(**)")

      fit.time[, 1] <- time
      fit.time[, 2] <- sapply(time, popsurv.fit)

      if(compute.cred == TRUE){
        credintS0pop <- t(sapply(time, credS0.pop))
        fit.time[, 3] <- credintS0pop[, 1]
        fit.time[, 4] <- credintS0pop[, 2]
      }
      fit.time <- round(fit.time, 4)
      outputlist <- list(fit.time = fit.time,
                         cred.int = cred.int,
                         covar.profile = covar.profile[-1])
      if(verbose == TRUE) {
        if(compute.cred == TRUE){
          cat("Estimated population survival at specified time points (*): \n")
          cat("\n")
          print.table(format(fit.time, nsmall = 4), right = TRUE)
          cat("--- \n")
          cat("* Bounds correspond to a", paste(format(round(cred.int  * 100, 2),
                                                       nsmall = 2), "%", sep = ""), "credible interval. \n")
          cat("** Population survival for covariate profile:",
              paste(round(covar.profile[-1], 4), collapse =", "),".")
        } else {
          cat("Estimated population survival at specified time points: \n")
          cat("\n")
          print.table(format(fit.time, nsmall = 4), right = TRUE)
          cat("--- \n")
          cat("** Population survival for covariate profile:",
              paste(round(covar.profile[-1], 4), collapse =", "),".")
        }
      }
      return(invisible(outputlist))
    } else if (type == "probacure"){
      # Output fitted results on supplied time vector
      n.time <- length(time)
      if(compute.cred == TRUE) ncols = 4 else ncols = 2
      fit.time <- matrix(0, nrow = n.time, ncol = ncols)
      if(compute.cred == TRUE){
        colnames(fit.time) <- c("Time", "Cure.prob(**)", "Cure.low", "Cure.up")
      } else colnames(fit.time) <- c("Time", "Cure.prob(**)")
      fit.time[, 1] <- time
      fit.time[, 2] <- sapply(time, probcure)

      if(compute.cred == TRUE){
        cred.probcure <- t(sapply(time, credprobcure))
        fit.time[, 3] <- cred.probcure[, 1]
        fit.time[, 4] <- cred.probcure[, 2]
      }
      fit.time <- round(fit.time, 4)
      outputlist <- list(fit.time = fit.time,
                         cred.int = cred.int,
                         covar.profile = covar.profile[-1])
      if(verbose == TRUE){
        if(compute.cred == TRUE){
          cat("Estimated cure prediction at specified time points (*): \n")
          cat("\n")
          print.table(format(fit.time, nsmall = 4), right = TRUE)
          cat("--- \n")
          cat("* Bounds correspond to a", paste(format(round(cred.int  * 100, 2),
                                                       nsmall = 2), "%", sep = ""), "credible interval. \n")
          cat("** Cure prediction for covariate profile:",
              paste(round(covar.profile[-1], 4), collapse =", "),".")
        } else {
          cat("Estimated cure prediction at specified time points: \n")
          cat("\n")
          print.table(format(fit.time, nsmall = 4), right = TRUE)
          cat("--- \n")
          cat("** Cure prediction for covariate profile:",
              paste(round(covar.profile[-1], 4), collapse =", "),".")
        }
      }
      return(invisible(outputlist))
    }
  }
}

