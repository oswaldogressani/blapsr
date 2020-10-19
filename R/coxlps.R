#' Fit a Cox proportional hazards regression model with Laplace-P-splines.
#'
#' Fits a Cox proportional hazards regression model for right censored data by
#' combining Bayesian P-splines and Laplace approximations.
#'
#'
#' @param formula A formula object where the ~ operator separates the response
#'   from the covariates. In a Cox model, it takes the form
#'   \emph{response ~ covariates}, where \emph{response} is a survival object
#'   returned by the \code{\link[survival]{Surv}} function of the \emph{survival}
#'   package.
#' @param data Optional. A data frame to match the variable names provided in
#'   \code{formula}.
#' @param K A positive integer specifying the number of cubic B-spline
#'   functions in the basis. Default is \emph{K = 30} and allowed values
#'   are \code{10 <= K <= 60}.
#' @param penorder The penalty order used on finite differences of the
#'  coefficients of contiguous B-splines. Can be either 2 for a second-order
#'  penalty (the default) or 3 for a third-order penalty.
#' @param tmax A user-specified value for the upper bound of the B-spline basis.
#'  The default is NULL, so that the B-spline basis is specified in the interval
#'  \emph{[0, tup]}, where \emph{tup} is the  upper bound of the follow-up
#'  times. It is required that \emph{tmax} > \emph{tup}.
#'
#'
#' @details The log-baseline hazard is modeled as a linear combination of
#'   \code{K} cubic B-splines as obtained from \code{\link{cubicbs}}. The
#'   B-spline basis is specified in the interval \emph{[0, tup]}, where
#'   \emph{tup} is the upper bound of the follow-up times,
#'   i.e. the largest observed follow-up time. Following Jullion
#'   and Lambert (2007), a robust Gamma prior is imposed on the roughness
#'   penalty parameter. A grid-based approach is used to explore the posterior
#'   penalty space and the resulting quadrature points serve to compute the
#'   approximate (joint) marginal posterior of the latent field vector. Point
#'   and set estimates of latent field elements are obtained from a finite
#'   mixture of Gaussian densities. The routine centers the columns of the
#'   covariate matrix around their mean value for numerical stability.
#'
#' @return An object of class \code{coxlps} containing various components from
#'   the fit. Details can be found in \code{\link{coxlps.object}}. Plot of
#'   estimated smooth hazard and survival curves can be obtained using
#'   \code{\link{plot.coxlps}}. If required, estimated baseline quantities
#'   on specific time values can be obtained with
#'   \code{\link{coxlps.baseline}}.
#'
#' @import survival
#' @examples
#'
#' ### Example 1 (Simulated survival data)
#'
#' set.seed(3)
#'
#' # Simulate survival data  with simsurvdata
#' betas <- c(0.13, 0.52, 0.30)
#' simul <- simsurvdata(a = 3.8, b = 2.2, n = 250, betas = betas , censperc = 20)
#' simul
#' simdat <- simul$survdata
#' plot(simul) # Plot survival data
#'
#' # Estimation with coxlps
#' fit <- coxlps(Surv(time, delta) ~ x1 + x2 + x3, data = simdat, K = 15)
#' # Compare coxlps and coxph
#' fit
#' summary(coxph(Surv(time, delta) ~ x1 + x2 + x3, data = simdat))
#'
#' # Fitted baseline survival vs target
#' plot(fit, h0 = FALSE, cred.int = 0.95, overlay.km = TRUE)
#' domt <- seq(0, 4, length = 100)
#' lines(domt, simul$S0(domt), type = "l", col = "red")
#' legend("topright", col=c("black", "blue", "red"), lty = rep(1,3),
#'       c("Bayesian LPS", "Kaplan-Meier", "Target"), cex = 0.8, bty = "n")
#'
#' ### Example 2 (Kidney transplant data)
#'
#' data(kidneytran)
#' Surv.obj <- Surv(kidneytran$time, kidneytran$delta)
#' fit <- coxlps(Surv.obj ~ age + gender + race, data = kidneytran)
#' coxphfit <- coxph(Surv.obj ~ age + gender + race, data = kidneytran)
#' ## Compare coxph and coxlps results
#' summary(coxphfit)
#' fit
#' ## Plot Kaplan-Meier curve vs Laplace-P-spline fit
#' plot(fit, h0 = FALSE, overlay.km = TRUE, plot.cred = FALSE)
#'
#' ### Example 3 (Laryngeal cancer data)
#'
#' data(laryngeal)
#' fit <- coxlps(Surv(time, delta) ~ age + diagyr + as.factor(stage),
#'                data = laryngeal)
#' coxphfit <- coxph(Surv(time, delta) ~ age + diagyr + as.factor(stage),
#'                   data = laryngeal)
#' ## Compare coxph and coxlps results
#' summary(coxphfit)
#' fit
#' ## Plot Kaplan-Meier curve vs Laplace-P-spline fit
#' plot(fit, h0 = FALSE, overlay.km = TRUE, plot.cred = FALSE)
#'
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @seealso \code{\link{Surv}}, \code{\link{coxph}}, \code{\link{simsurvdata}},
#'  \code{\link{coxlps.object}}, \code{\link{coxlps.baseline}}.
#'
#' @references Cox, D.R. (1972). Regression models and life-tables.
#'   \emph{Journal of the Royal Statistical Society: Series B (Methodological)}
#'   \strong{34}(2): 187-202. \url{https://doi.org/10.1111/j.2517-6161.1972.tb00899.x}
#' @references Gressani, O. and Lambert, P. (2018). Fast Bayesian inference
#'   using Laplace approximations in a flexible promotion time cure model based
#'   on P-splines. \emph{Computational Statistical & Data Analysis} \strong{124}:
#'   151-167. \url{https://doi.org/10.1016/j.csda.2018.02.007}
#' @references Jullion, A. and Lambert, P. (2007). Robust specification of the
#'   roughness penalty prior distribution in spatially adaptive Bayesian
#'   P-splines models. \emph{Computational Statistical & Data Analysis}
#'   \strong{51}(5): 2542-2558. \url{https://doi.org/10.1016/j.csda.2006.09.027}
#'
#' @export

coxlps <- function(formula, data, K = 30, penorder = 2, tmax = NULL){

  if(!inherits(formula, "formula"))
    stop("Incorrect model formula")
  if (missing(data)) {
    mf <- stats::model.frame(formula)  # Extract model frame from formula
  } else{
    mf <- stats::model.frame(formula, data = data)
  }
  if (!survival::is.Surv(mf[, 1]))
    stop("Response in formula is not an object of class Surv")
  if(missing(data)){
    X  <- stats::model.matrix(mf) # Design matrix
  } else{
    X <- stats::model.matrix(mf, data = data)
  }
  colnames.X <- colnames(X)
  if (colnames(X)[1] == "(Intercept)") {
    X <- as.matrix(X[,-1]) # Design matrix
    colnames(X) <- colnames.X[-1]
  }
  if (any(is.infinite(X)))
    stop("Data contains infinite covariates")
  if (ncol(X) == 0)
    stop("The model has no covariates")
  y <- matrix(stats::model.extract(mf, "response"), ncol = 2) # Response
  if (any(is.infinite(y[, 1]))  || any(y[, 1] < 0))
    stop("Event times contain Inf/negative values")
  sd.time <- stats::sd(y[, 1])     # Standard deviation of event times
  time <- y[, 1] / sd.time  # Standardized event times
  event <- y[, 2]           # Event indicator
  tlow <- 0                 # Lower bound of follow-up
  tup <- max(time)          # Upper bound of follow-up
  if(!is.null(tmax)){
    if(tmax <=0) stop("tmax must be positive")
    else if((tmax / sd.time) < tup)
      stop("tmax must be larger than the largest observed event time")
    tup <- tmax / sd.time
  }
  nbeta <- ncol(X)          # Number of covariates in the Cox model
  n <- length(time)         # Sample size
  Xdata <- X                # Design matrix
  X <- scale(X, center = TRUE, scale = FALSE) # Mean centered design matrix
  if (!is.vector(K, mode = "numeric") || length(K) > 1  || is.na(K))
    stop("K must be numeric of length 1")
  if (K < 10 || K > 60)
    stop("Number K of B-splines in basis should be between 10 and 60")
  H <- K + nbeta  # Latent field dimension
  if (n < H)
    warning("Number of coefficients to be estimated is larger than sample size")
  penorder <- floor(penorder)
  if(penorder < 2 || penorder > 3)
    stop("Penalty order must be either 2 or 3")

  # Bins and midpoints to estimate baseline survival
  bins <- 300   # Total number of bins
  partition <- seq(tlow, tup, length = bins + 1)  # Partition domain
  width <- partition[2] - partition[1]            # Bin width
  middleBins <- partition[1:bins] + (width / 2)   # Middle points of the bins
  upto <- as.integer(time / width) + 1   # Indicates in which bin time falls
  upto[which(upto == bins + 1)] <- bins  # If time == tup (bin 301 to 300)
  # Matrix matsumtheta required for computing Hessian of log likelihood
  fill.ones <- function(vecidx) c(rep(1, vecidx[1]), rep(0, vecidx[2]))
  m.idx <- matrix(cbind(upto, (bins - upto)), ncol = 2)
  matsumtheta <- matrix(apply(m.idx, 1, fill.ones), nrow = n, ncol = bins,
                        byrow = TRUE) * width

  # B-spline basis
  Bobs <- cubicbs(time, tlow, tup, K)$Bmatrix
  Bmiddle <- cubicbs(middleBins, tlow, tup, K)$Bmatrix
  Bdesign <- cbind(Bobs, X)     # Full design matrix
  eventBobs <- event %*% Bobs   # Product of event indicator and spline basis

  # Penalty matrix
  D <- diag(K)
  for(k in 1:penorder){D <- diff(D)}
  P <- t(D) %*% D           # Penalty matrix of dimension c(K,K)
  P <- P + diag(1e-06, K)   # Perturbation to ensure P full rank

  # Robust prior following Jullion and Lambert (2007)
  a.delta <- 1e-04
  b.delta <- 1e-04
  nu <- 3
  prec.betas <- 1e-05 # Precision for regression coefficients

  # Log-likelihood function
  loglik <- function(latent) {
    val <- sum(event * (Bdesign %*% latent)) -
      sum(cumsum(exp((Bmiddle %*% latent[1:K])) * width)[upto] *
            exp(X %*% latent[(K + 1):H]))
    return(val)
  }

  # Latent field prior precision matrix as a function of v = log(lambda)
  Q.v <- function(v) {
    Qmatrix <- matrix(0, nrow = H, ncol = H)
    Qmatrix[1:K, 1:K] <- exp(v) * P
    Qmatrix[(K + 1):H, (K + 1):H] <- diag(prec.betas, nbeta)
    return(Qmatrix)
  }

  # Log conditional posterior of latent field given v
  logplat <- function(latent, v) {
    Qv <- Q.v(v)
    y <- as.numeric(loglik(latent) - .5 * sum((latent * Qv) %*% latent))
    return(y)
  }

  # Gradient of the log conditional posterior latent field
  grad.logplat <- function(latent, v) {
    # Gradient of the log-likelihood
    theta.vector <- latent[1:K]
    beta.vector <- latent[(K + 1):H]
    expBmidtheta <- as.numeric(exp(Bmiddle %*% theta.vector))
    expbetaX <- as.numeric(exp(X %*% beta.vector))
    prodwidth <- expBmidtheta * Bmiddle * width
    grad.theta <- eventBobs - as.numeric(expbetaX %*%
                                           apply(prodwidth, 2, cumsum)[upto, ])
    grad.beta <-
      as.numeric((event - cumsum(expBmidtheta * width)[upto] *
                    expbetaX) %*% X)
    grad.loglik <- c(grad.theta, grad.beta)

    # Gradient of the log conditional posterior latent field given v
    Qv <- Q.v(v)
    grad.vec <- as.numeric(grad.loglik - (Qv %*% latent))
    return(grad.vec)
  }

  # Hessian of the log conditional posterior latent field
  Hess.logplat <- function(latent, v) {
    # Hessian of the log-likelihood
    theta.vector <- latent[1:K]
    beta.vector <- latent[(K + 1):H]
    expBmidtheta <- as.numeric(exp(Bmiddle %*% theta.vector))
    expbetaX <- as.numeric(exp(X %*% beta.vector))

    # Upper left block
    colvecsum <- as.numeric(expbetaX %*% matsumtheta)
    upleft.block <- (-1) * (t(Bmiddle) %*% (expBmidtheta * colvecsum * Bmiddle))

    # Upper right block
    matprod.upright <- expBmidtheta * Bmiddle * width
    upright.block <- (-1) * (t(apply(matprod.upright, 2, cumsum)[upto,] *
                                 expbetaX) %*% X)
    # Lower left block
    lowleft.block <- t(upright.block)

    # Lower right block
    lowright.block <- (-1) * t(X * (cumsum(expBmidtheta * width)[upto] *
                                      expbetaX)) %*% X

    # Hessian matrix
    Hess.loglik <- matrix(0, nrow = H, ncol = H)
    Hess.loglik[1:K, 1:K] <- upleft.block
    Hess.loglik[1:K, (K + 1):H] <- upright.block
    Hess.loglik[(K + 1):H, 1:K] <- lowleft.block
    Hess.loglik[(K + 1):H, (K + 1):H] <- lowright.block

    # Hessian of the log conditional posterior latent field given v
    Qv <- Q.v(v)
    Hess.mat <- Hess.loglik - Qv
    return(list(Hess.mat = Hess.mat, Hess.loglik = Hess.loglik))
  }

  # Log-posterior of log-penalty v
  logpost.v <- function(v, latent0) {
    # 1) Given log-penalty v, compute posterior mode of Laplace approximation
    tol <- 1e-3            # Tolerance
    dist <- 3              # Initial distance
    iter <- 0              # Iteration counter

    while (dist > tol) {
      dlat <- (-1) * solve(Hess.logplat(latent0, v)$Hess.mat,
                           grad.logplat(latent0, v))
      lat.new <- latent0 + dlat
      step <- 1
      iter.halving <- 1
      logplat.current <- logplat(latent0, v)
      while (logplat(lat.new, v) <= logplat.current) {
        step <- step * .5
        lat.new <- latent0 + (step * dlat)
        iter.halving <- iter.halving + 1
        if (iter.halving > 10) {
          break
        }
      }
      dist <- sqrt(sum((lat.new - latent0) ^ 2))
      iter <- iter + 1
      latent0 <- lat.new
    }

    latentstar <- latent0   # Posterior mode of Laplace approximation given v
    Htilde <- Hess.logplat(latentstar, v)$Hess.loglik # Hessian of log-lik.

    # 2) Given latentstar, compute the log-posterior of v
    Qv <- Q.v(v)
    L <- t(chol(Qv - Htilde))

    y <- loglik(latentstar) -
      (.5 * sum((Qv * latentstar) %*% latentstar)) +
      (.5 * (K + nu) * v) -
      sum(log(diag(L))) -
      (a.delta + .5 * nu) * log(b.delta + .5 * nu * exp(v))
    Covmatrix <- solve(Qv - Htilde)

    return(list(y = y, # Log-posterior of v
                latentstar = latentstar, # Posterior mode of Laplace approx.
                Sigstar.theta = Covmatrix[1:K, 1:K], # Covar. of spline coeffs.
                varbetas = diag(Covmatrix)[(K + 1):H])) # Var. of beta coeffs.
  }

  # Grid construction
  grid.construct <- function() {
    v.grid <- c()
    logpv.v <- c()
    iter <- 1
    jump <- 3
    v.0 <- 20
    penlow <- 0
    latent0 <- rep(sum(event) / sum(time), H) # Initial latent field vector
    logpv <- logpost.v(v.0, latent0)
    logpv.v[iter] <- logpv$y
    latent0 <- logpv$latentstar
    v.grid[iter] <- v.0

    while (jump > 0) {
      iter <- iter + 1
      v.0 <- v.0 - 1
      v.grid[iter] <- v.0
      logpv <- logpost.v(v.0, latent0)
      logpv.v[iter] <- logpv$y
      latent0 <- logpv$latentstar
      jump <- logpv.v[iter] -  logpv.v[iter - 1]
      if (v.grid[iter] == 0) {
        penlow <- 1
        break
      }
    }

    v.grid <- rev(v.grid)
    logpv.v <- rev(logpv.v)
    crit.val <- exp(-.5 * stats::qchisq(0.9, df = 1, lower.tail = TRUE))

    if (penlow == 1) {
      vmap <- 0
      logpvmap <- logpv.v[1]
      vquadrature <- v.grid[which((logpv.v / logpvmap) > crit.val)]
      dvquad <- vquadrature[2] - vquadrature[1]
      list.star <- lapply(vquadrature, logpost.v, latent0 = latent0)
      logpvquad <- unlist(lapply(list.star, "[[", 1))
      nquad <- length(vquadrature)
      omega.m <- exp(logpvquad - logpvmap) * dvquad /
        sum(exp(logpvquad - logpvmap) * dvquad)
    } else{
      # Finer nesting of vmap
      nest1.lb <- v.grid[1]
      nest1.ub <- v.grid[3]
      jump <- 3
      iter <- 1
      v.0 <- nest1.ub
      v.grid <- c()
      logpv.v <- c()
      v.grid[iter] <- v.0
      logpv <- logpost.v(v.0, latent0)
      logpv.v[iter] <- logpv$y
      latent0 <- logpv$latentstar

      while (jump > 0) {
        iter <- iter + 1
        v.0 <- v.0 - 0.1
        v.grid[iter] <- v.0
        logpv <- logpost.v(v.0, latent0)
        logpv.v[iter] <- logpv$y
        latent0 <- logpv$latentstar
        jump <- logpv.v[iter] -  logpv.v[iter - 1]
      }

      v.grid <- rev(v.grid)
      logpv.v <- rev(logpv.v)

      vmap <- (v.grid[1] + v.grid[3]) / 2
      logpvmap <- logpost.v(vmap, latent0)$y
      nquad <- 10
      left.step.j  <- (-1)
      right.step.j <- 1

      # Grid construction
      v.left.j <- vmap
      ratio <- 1
      while (ratio > crit.val) {
        v.left.j <- v.left.j + left.step.j
        ratio <- exp(logpost.v(v.left.j, latent0)$y - logpvmap)
      }
      v.right.j <- vmap
      ratio <- 1
      while (ratio > crit.val) {
        v.right.j <- v.right.j + right.step.j
        ratio <- exp(logpost.v(v.right.j, latent0)$y - logpvmap)
      }
      lb.j <- v.left.j
      ub.j <- v.right.j
      vquadrature <- seq(lb.j, ub.j, length = nquad)
      dvquad <- vquadrature[2] - vquadrature[1]
      list.star <- lapply(vquadrature, logpost.v, latent0 = latent0)
      logpvquad <- unlist(lapply(list.star, "[[", 1))
      omega.m <- exp(logpvquad - logpvmap) * dvquad /
        sum(exp(logpvquad - logpvmap) * dvquad)
    }
    listout <- list(vmap = vmap,
                    logpvmap = logpvmap,
                    nquad = nquad,
                    vquadrature = vquadrature,
                    logpvquad = logpvquad,
                    omega.m = omega.m,
                    list.star = list.star,
                    latent.0 = latent0)
    return(invisible(listout))
  }

  gridpoints <- grid.construct()
  vmap <- gridpoints$vmap
  vquadrature <- gridpoints$vquadrature
  nquad <- gridpoints$nquad
  logpvquad <- gridpoints$logpvquad
  omega.m <- gridpoints$omega.m
  list.star <- gridpoints$list.star

  # Posterior inference
  latent.matrix <- matrix(unlist(lapply(list.star, "[[", 2)),
                          ncol = H, byrow = TRUE)
  varbetas <- matrix(unlist(lapply(list.star, "[[", 4)),
                     ncol = nbeta, byrow = TRUE)
  # Posterior point estimate of latent field vector
  latent.hat <- colSums(omega.m * latent.matrix)
  theta.hat <- latent.hat[1:K]
  betas.hat <- latent.hat[(K + 1):H]

  # Posterior standard deviation of regression coefficients
  sdbeta.post <- sqrt(colSums(omega.m * varbetas) +
                        colSums(omega.m * (latent.matrix[, (K + 1):H] -
                        matrix(rep(betas.hat, nquad), nrow = nquad,
                         byrow = TRUE)) ^ 2))

  z.Wald <- betas.hat / sdbeta.post # Wald test statistic
  Sigstar.list <- lapply(list.star, "[[", 3)
  diffthetas <- list()

  for(m in 1 : nquad){
    Sigstar.list[[m]] <- omega.m[m] * Sigstar.list[[m]]
    diffthetas[[m]] <- omega.m[m] * outer(latent.matrix[m, 1:K] - theta.hat,
                                          latent.matrix[m, 1:K] - theta.hat)
  }

  Covthetamix <- Reduce("+", Sigstar.list) + Reduce("+", diffthetas)

  # Posterior credible intervals for regression coefficients
  postbetas.CI95 <- matrix(0, nrow = nbeta, ncol = 2)

  for (j in 1:nbeta) {
    means.betaj <- latent.matrix[, (K + j)]
    variances.betaj <- varbetas[, j]
    betaj.dom <- seq(betas.hat[j] - 4 * sdbeta.post[j],
                     betas.hat[j] + 4 * sdbeta.post[j],
                     length = 500)
    dj <- betaj.dom[2] - betaj.dom[1]

    post.betaj <- function(betaj) {
      mixture <- sum(omega.m * stats::dnorm(betaj, mean = means.betaj,
                                     sd = sqrt(variances.betaj)))
      return(mixture)
    }
    mix.image <- sapply(betaj.dom, post.betaj)
    cumsum.image <- cumsum(mix.image * dj)
    lb95.j <- betaj.dom[sum(cumsum.image < 0.025)]
    ub95.j <- betaj.dom[sum(cumsum.image < 0.975)]
    postbetas.CI95[j, ] <- c(lb95.j, ub95.j)
  }

  # Output

  # Regression coefficients
  regcoeff.output <- matrix(c(betas.hat, exp(betas.hat), sdbeta.post,
                              z.Wald, exp(betas.hat), exp( - betas.hat),
                              exp(postbetas.CI95[, 1]),
                              exp(postbetas.CI95[,2])), ncol = 8)
  colnames(regcoeff.output) <- c("coef", "exp(coef)","sd.post", "z","exp(coef)",
                                 "exp(-coef)", "lower.95", "upper.95")
  rownames(regcoeff.output) <- colnames(X)
  penvalues <- exp(vquadrature) # Penalty values
  splinecoeffs <- theta.hat     # Estimated B-spline coefficients
  I.loglik <- (-1) * Hess.logplat(latent.hat, vmap)$Hess.loglik
  EDmat <- solve(I.loglik + Q.v(vmap)) %*% I.loglik
  edf <- diag(EDmat)
  names(edf) <- c(paste0(rep("spline.coeff"), seq_len(K)), colnames(X))
  ED <- sum(edf) # Effective dimension of entire model
  p <- ncol(X)
  loglik.value <- loglik(latent.hat)
  AIC.p <-  (-2 * loglik.value) + 2 * p
  AIC.ED <- (-2 * loglik.value) + 2 * ED
  BIC.p <-  (-2 * loglik.value) + p * log(sum(event))
  BIC.ED <- (-2 * loglik.value) + ED * log(sum(event))

   listout <- list(formula = formula,           # Model formula
                   K = K,                       # Number of B-splines in basis
                   penalty.order = penorder,    # Penalty order
                   latfield.dim = H,            # Latent field dimension
                   n = n,                       # Sample size
                   num.events = sum(event),     # Number of events
                   tup = tup,                   # Upper bound of follow-up
                   event.times = time,          # Standardized event times
                   sd.time = sd.time,           # Stand. dev. of original times
                   event.indicators = event,    # Event indicators
                   regcoeff = regcoeff.output,  # Estimated regression coeff.
                   penalty.vector = penvalues,  # Grid of penalty values
                   vmap = vmap,                 # Maximum a posteriori for v
                   spline.estim = splinecoeffs, # Estimated B-spline coeffs.
                   edf = round(edf, 4),         # Individual degrees of freedom
                   ED = ED,                     # Effective model dimension
                   Covthetamix = Covthetamix,   # Covariance of spline coeffs.
                   X = Xdata,                   # Design matrix
                   loglik = loglik.value,       # Log-likelihood at estimated parameters
                   p = p,                       # Number of parametric terms
                   AIC.p  = AIC.p,              # AIC with penalty term p
                   AIC.ED = AIC.ED,             # AIC with penalty term ED
                   BIC.p  = BIC.p,              # BIC with penalty term p
                   BIC.ED = BIC.ED)             # BIC with penalty term ED
  attr(listout, "class") <- "coxlps"
  listout
}








