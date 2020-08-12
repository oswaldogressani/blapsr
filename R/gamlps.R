#' Bayesian generalized additive modeling with Laplace-P-splines.
#'
#' @description {Fits a generalized additive model (GAM) to data using an
#' approximate Bayesian inference technique based on penalized regression
#' splines and Laplace approximations. Smooth additive terms are specified as a
#' linear combination of a large number of cubic B-splines. To counterbalance
#' the roughness of the fit, a discrete penalty on neighboring spline
#' coefficients is imposed in the spirit of Eilers and Marx (1996). The
#' effective degrees of freedom of the smooth terms are also estimated.
#'
#' The optimal amount of smoothing is determined by a grid-based exploration of
#' the posterior penalty space when the number of smooth terms is small to
#' moderate. When the dimension of the penalty space is large, the optimal
#' smoothing parameter is chosen to be the value that maximizes the
#' (log-)posterior of the penalty vector. Approximate Bayesian credible
#' intervals for latent model variables and functions of latent model variables
#' are available.
#' }
#'
#' @usage gamlps(formula, data, K = 30, family = c("gaussian", "poisson", "bernoulli", "binomial"),
#'        gauss.scale, nbinom, penorder = 2, cred.int = 0.95)
#'
#' @param formula A formula object where the ~ operator separates the response
#'   from the covariates of the linear part \code{z1,z2,..} and the smooth
#'   terms. A smooth term is specified by using the notation \code{sm(.)}.
#'   For instance, the formula \code{y ~ z1+sm(x1)+sm(x2)} specifies a
#'   generalized additive model of the form
#'   \emph{g(mu) = b0+b1z1+f1(x1)+f2(x2)}, where \emph{b0, b1} are the
#'   regression coefficients of the linear part and \emph{f1(.)} and
#'   \emph{f2(.)} are smooth functions of the continuous covariates \emph{x1}
#'   and \emph{x2} respectively. The function \emph{g(.)} is the canonical
#'   link function.
#' @param data Optional. A data frame to match the variables names provided in
#'   formula.
#' @param K A positive integer specifying the number of cubic B-spline
#'   functions in the basis used to model the smooth terms. Default is
#'   \code{K = 30} and allowed values are \code{15 <= K <= 60}. The same basis
#'   dimension is used for each smooth term in the model. Also, the
#'   computational cost to fit the model increases with \code{K}.
#' @param family The error distribution of the model. It is a character string
#'   that must partially match either \code{"gaussian"} for Normal data with
#'   an identity link, \code{"poisson"} for Poisson data with a log link,
#'   \code{"bernoulli"} or \code{"binomial"} for Bernoulli or Binomial data
#'   with a logit link.
#' @param gauss.scale The scale parameter to be specified when
#'   \code{family = "gaussian"}. It corresponds to the variance of the response.
#' @param nbinom The number of experiments in the Binomial family.
#' @param penorder The penalty order used on finite differences of the
#'  coefficients of contiguous B-splines. Can be either 2 for a second-order
#'  penalty (the default) or 3 for a third-order penalty.
#' @param cred.int The level of the pointwise credible interval to be computed
#'  for the coefficients in the linear part of the model.
#'
#' @details {The B-spline basis used to approximate a smooth additive component
#'  is computed with  the function \code{\link{cubicbs}}. The lower (upper)
#'  bound of the B-spline basis is taken to be the minimum (maximum) value of
#'  the covariate associated to the smooth. For identifiability
#'  purposes, the B-spline matrices (computed over the observed covariates)
#'  are centered. The centering consists is subtracting from each column of the
#'  B-spline matrix, the corresponding column average of another B-spline matrix
#'  computed on a fine grid of equidistant values in the domain of the smooth
#'  term.
#'
#'  A hierarchical Gamma prior is imposed on the roughness penalty vector. A
#'  Newton-Raphson algorithm is used to compute the posterior
#'  mode of the (log-)posterior penalty vector. The latter algorithm uses
#'  analytically derived versions of the gradient and Hessian. When the number
#'  of smooth terms in the model is smaller or equal to 4, a grid-based strategy
#'  is used for posterior exploration of the penalty space. Above that
#'  threshold, the optimal amount of smoothness is determined by the posterior
#'  maximum value of the penalty vector. This strategy allows to keep the
#'  computational burden to fit the model relatively low and conserve good
#'  statistical performance.}
#'
#' @return An object of class \code{gamlps} containing several components from
#'   the fit. Details can be found in \code{\link{gamlps.object}}. Details on
#'   the output printed by \code{gamlps} can be found in
#'   \code{\link{print.gamlps}}. Fitted smooth terms can be visualized with the
#'   \code{\link{plot.gamlps}} routine.
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @seealso \code{\link{cubicbs}}, \code{\link{gamlps.object}},
#'  \code{\link{print.gamlps}}, \code{\link{plot.gamlps}}
#'
#' @examples
#'
#' set.seed(14)
#' sim <- simgamdata(n = 300, setting = 2, dist = "binomial", scale = 0.25)
#' dat <- sim$data
#' fit <- gamlps(y ~ z1 + z2 + sm(x1) + sm(x2), K = 15, data = dat,
#'               penorder = 2, family = "binomial", nbinom = 15)
#' fit
#'
#' # Check fit and compare with target (in red)
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 2))
#' domx <- seq(-1, 1 ,length = 200)
#' plot(fit, smoo.index = 1, cred.int = 0.95, ylim = c(-2, 2))
#' lines(domx, sim$f[[1]](domx), type= "l", lty = 2, lwd = 2, col = "red")
#' plot(fit, smoo.index = 2, cred.int = 0.95, ylim = c(-3, 3))
#' lines(domx, sim$f[[2]](domx), type= "l", lty = 2, lwd = 2, col = "red")
#' par(opar)
#'
#' @references Hastie, T. J. and Tibshirani., RJ (1990). Generalized additive
#'   models. \emph{Monographs on statistics and applied probability}, 43, 335.
#' @references Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with
#'  B-splines and penalties. \emph{Statistical Science}, \strong{11}(2): 89-121.
#' @references Gressani, O. and Lambert, P. (2018). Fast Bayesian inference
#'  using Laplace approximations in a flexible promotion time cure model based
#'  on P-splines. \emph{Computational Statistical & Data Analysis} \strong{124}:
#'  151-167.\url{https://doi.org/10.1016/j.csda.2018.02.007}

#'
#' @export

gamlps <- function(formula, data, K = 30, family = c("gaussian", "poisson",
                    "bernoulli", "binomial"), gauss.scale, nbinom,
                    penorder = 2, cred.int = 0.95){

  if (!inherits(formula, "formula"))
    stop("Incorrect model formula")
  if(missing(data)) {
    mf <- stats::model.frame(formula) # Extract model frame from formula
    X  <- stats::model.matrix(mf)     # Full design matrix
    colXnames <- colnames(X)
    smterms <- grepl("sm(", colnames(X), fixed = TRUE)
    X <- cbind(X[, as.logical(1 - smterms)], X[, smterms])
    colnames(X) <- colXnames
  } else{
    mf <- stats::model.frame(formula, data = data)
    X <- stats::model.matrix(mf, data = data)
    colXnames <- colnames(X)
    smterms <- grepl("sm(", colnames(X), fixed = TRUE)
    X <- cbind(X[, as.logical(1 - smterms)], X[, smterms])
    colnames(X) <- colXnames
  }
  if(any(is.infinite(X)))
    stop("Covariates contain Inf values")
  q  <- sum(smterms) # Number of smooth terms in model
  if(q == 0)
    stop("Model does not contain any smooth terms")
  p  <- ncol(X) - q # Number of regression coefficients in linear part
  n  <- nrow(X)     # Sample size
  Z  <- scale(X[, 1:p], center = TRUE, scale = FALSE)  # Centered Z matrix
  Z[, 1] <- rep(1, n) # Column for intercept
  if(ncol(Z)==1) colnames(Z) <- "(Intercept)"
  y  <- as.numeric(stats::model.extract(mf, "response")) # Response vector
  if(any(is.infinite(y)) || any(is.na(y)))
    stop("Response contains Inf, NA or NaN values")
  if(!is.vector(K, mode = "numeric") || length(K) > 1 || is.na(K))
    stop("K must be numeric of length 1")
  if(K < 15 || K > 60)
    stop("K must be between 15 and 60")
  if(penorder < 2 || penorder > 3)
    stop("Penalty order must be either 2 or 3")
  if(!missing(gauss.scale)){
    if(gauss.scale <= 0)
      stop("Gaussian scale parameter must be positive")
  }
  if(!missing(nbinom)){
    if(nbinom <= 1){
      stop("nbinom must be larger than one")
    }
  }

  # Chosen family
  EF.type <- match.arg(family)
  if (EF.type == "gaussian") {
    family <- "gaussian"
    link <- "identity"
  }
  if (EF.type == "poisson") {
    family <- "poisson"
    link <- "log"
  }
  if (EF.type == "bernoulli") {
    family <- "bernoulli"
    link <- "logit"
  }
  if (EF.type == "binomial") {
    family <- "binomial"
    link <- "logit"
  }

  B.list <- list() # List of B-spline bases for smooth terms
  # Centered B-spline matrices
  for(j in 1:q){
    xj <- as.numeric(X[, p + j]) # values of the jth covariate
    min.xj <- min(xj) # minimum
    max.xj <- max(xj) # maximum
    xj.fine <- seq(min.xj, max.xj, length = 1000) # fine grid
    Bj.fine <- cubicbs(xj.fine, lower = min.xj, upper = max.xj, K = K)$Bmatrix
    Bj.fine.mean <- colMeans(Bj.fine)
    Bj <- cubicbs(xj, lower = min.xj, upper = max.xj, K = K)$Bmatrix
    Bj.centered <- Bj - matrix(rep(Bj.fine.mean, n), nrow = n, byrow = TRUE)
    B.list[[j]] <- Bj.centered
  }
  B.list.trim <- lapply(B.list, function(x) x[, -K])
  B <- cbind(Z, do.call(cbind, B.list.trim)) # Design matrix ncol: q * (K-1) + p
  H <- p + (q * (K - 1))  # Latent field dimension
  if(n < H)
    stop("Number of coefficients to be estimated is larger than sample size")

  # Penalty matrix
  D <- diag(K) # Diagonal matrix
  for (k in 1:penorder) D <- diff(D) # Difference matrix of dimension K-r by K
  D <- D[, -K]                       # Delete last column for identifiability
  P <- t(D) %*% D                    # Penalty matrix of dimension K-1 by K-1
  P <- P + diag(1e-12, (K - 1))      # Perturbation to make P full rank

  # Robust prior on penalty parameters
  a.delta <- 0.5
  b.delta <- 0.5
  nu <- 1
  prec.betas <- diag(1e-05, p) # Precision for linear regression coefficients

  ## Matrices related to the chosen family

  if(family == "gaussian") {
    gam.scale <- gauss.scale
    invgam.scale <- (1 / gam.scale)
    s.gam <- function(latent) .5 * (as.numeric(B %*% latent) ^ 2)
    s.gam.deriv <- function(pred) pred
    s.gam.deriv2 <- function(pred) 1
    W.gam <- function(latent) (1 / gam.scale) * diag(1, n)
    mu.gam <- function(latent) as.numeric(B %*% latent)
  }

  if(family == "poisson") {
    gam.scale <- 1
    invgam.scale <- 1
    s.gam <- function(latent) as.numeric(exp(B %*% latent))
    s.gam.deriv <- function(pred) exp(pred)
    s.gam.deriv2 <- function(pred) exp(pred)
    W.gam <- function(latent) diag(exp(as.numeric(B %*% latent)))
    mu.gam <- function(latent) exp(as.numeric(B %*% latent))
  }

  if(family == "bernoulli") {
    gam.scale <- 1
    invgam.scale <- 1
    s.gam <- function(latent) log(1 + exp(as.numeric(B %*% latent)))
    s.gam.deriv <- function(pred) exp(pred) / (1 + exp(pred))
    s.gam.deriv2 <- function(pred) exp(pred) / ((1 + exp(pred)) ^ 2)
    W.gam <- function(latent){
      predictor <- as.numeric(B %*% latent)
      val <- diag(exp(predictor) / ((1 + exp(predictor)) ^ 2))
      return(val)
    }
    mu.gam <- function(latent) {
      predictor <- as.numeric(B %*% latent)
      val <- exp(predictor) / (1 + exp(predictor))
      return(val)
    }
  }

  if(family == "binomial") {
    gam.scale <- 1
    invgam.scale <- 1
    s.gam <- function(latent) nbinom * log(1 + exp(as.numeric(B %*% latent)))
    s.gam.deriv <- function(pred) nbinom * (exp(pred) / (1 + exp(pred)))
    s.gam.deriv2 <- function(pred) nbinom * (exp(pred) / ((1 + exp(pred)) ^ 2))
    W.gam <- function(latent){
      predictor <- as.numeric(B %*% latent)
      val <- diag(nbinom * exp(predictor) / ((1 + exp(predictor)) ^ 2))
      return(val)
    }
    mu.gam <- function(latent) {
      predictor <- as.numeric(B %*% latent)
      val <- nbinom * exp(predictor) / (1 + exp(predictor))
      return(val)
    }
  }

  # Prior precision matrix of latent field for a given v
  Qv <- function(v) as.matrix(Matrix::bdiag(prec.betas, diag(exp(v), q) %x% P))

  BWD <- t(B) %*% diag(invgam.scale, n) # To compute log-likelihood gradient

  # Log conditional posterior of latent field given v
  logplat <- function(latent, v) {
    post <- invgam.scale * sum(y * as.numeric(B %*% latent) -
                          s.gam(latent)) - .5 * t(latent) %*% Qv(v) %*% latent
    return(post)
  }

  # Gradient of the log conditional posterior latent field
  grad.logplat <- function(latent, v) {
    grad.loglik <- BWD %*% (y - mu.gam(latent))
    grad.vec <- as.numeric(grad.loglik - Qv(v) %*% latent)
    return(grad.vec)
  }

  # Hessian of the log conditional posterior latent field
  Hess.logplat <- function(latent, v) {
    Hess.loglik <- (-1) * t(B) %*% W.gam(latent) %*% B
    Hess.mat <- Hess.loglik - Qv(v)
    return(Hess.mat)
  }

  # Laplace approximation to conditional latent field
  Laplace <- function(latent0, v){

    epsilon <- 1e-05 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter

    for (k in 1:maxiter) {
      dlat <- as.numeric((-1) * solve(Hess.logplat(latent0, v),
                                      grad.logplat(latent0, v)))
      lat.new <- latent0 + dlat
      step <- 1
      iter.halving <- 1
      logplat.current <- logplat(latent0, v)
      while (logplat(lat.new, v) <= logplat.current) {
        step <- step * .5
        lat.new <- latent0 + (step * dlat)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((lat.new - latent0) ^ 2))
      iter <- iter + 1
      latent0 <- lat.new
      if(dist < epsilon) break
    }

    latentstar <- latent0   # Posterior mode of Laplace approximation given v
    Wtilde <- W.gam(latentstar)
    Sigmastar <- solve(t(B) %*% Wtilde %*% B + Qv(v))
    return(list(latentstar = latentstar,
                Sigmastar = Sigmastar,
                Wtilde = Wtilde))
  }

  # Log-posterior of log-penalty vector
  logp.vD0 <- function(v, laplace.approx) {
    Q <- Qv(v)
    latentstar <- laplace.approx$latentstar
    Sigmastar <- laplace.approx$Sigmastar
    Wtilde <- laplace.approx$Wtilde
    L <- t(chol((t(B) %*% (diag(Wtilde) * B)) + Q))
    Mv <- chol2inv(t(L))
    Blat <- as.numeric(B %*% latentstar)

    lpost <- (-1) * sum(log(diag(L))) +
      .5 * (nu + K - 1) * sum(v) +
      invgam.scale * sum(y * Blat) -
      invgam.scale * sum(s.gam(latentstar)) -
      .5 * as.numeric(t(latentstar) %*% Q %*% latentstar) -
      (.5 * nu + a.delta) * sum(log(b.delta + .5 * nu * exp(v)))

    listout <- list(value = lpost,
                    latentstar = latentstar,
                    Sigmastar = Sigmastar)

    return(listout)
  }

  # Gradient of logp.vD
  gradient <- function(v, laplace.approx){

    Q <- Qv(v)
    latentstar <- laplace.approx$latentstar
    #Sigmastar <- laplace.approx$Sigmastar
    Wtilde <- laplace.approx$Wtilde
    Blat <- as.numeric(B %*% latentstar)
    wbar <- as.numeric(t(B) %*% (invgam.scale * (y - mu.gam(latentstar)) +
                                   Wtilde %*% Blat))
    L <- t(chol((t(B) %*% (diag(Wtilde) * B)) + Q))
    Mv <- chol2inv(t(L))

    Pvj <- function(j){
      as.matrix(Matrix::bdiag(matrix(0, nrow = p, ncol = p),
                              diag(diag(q)[, j] * exp(v), q) %x% P))
    }
    partial.v <- function(j) {
      Gamma.vj <- Mv %*% Pvj(j) %*% Mv
      Gwbar <- Gamma.vj %*% wbar
      Bgw <- as.numeric(B %*% Gwbar)
      Mwbar <- Mv %*% wbar
      val <- (-.5) * sum(Mv[((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1)),
                            ((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1))] *
                           (exp(v[j]) * P)) +
        (.5 * (nu + (K - 1))) -
        invgam.scale * sum(y * Bgw) +
        invgam.scale * sum(s.gam.deriv(as.numeric(B %*% Mwbar)) * Bgw) +
        as.numeric(t(wbar) %*% Gamma.vj %*% Q %*% Mwbar) -
        as.numeric(.5 * (t(wbar) %*% Gwbar)) -
        (.5 * nu + a.delta) / (1 + (2 * b.delta) / (nu * exp(v[j])))
      return(val)
    }

    partial.vec <- sapply(seq_len(q), partial.v)
    return(partial.vec)
  }

  # Hessian of logp.vD
  Hessian <- function(v, laplace.approx) {
    Q <- Qv(v)
    latentstar <- laplace.approx$latentstar
    #Sigmastar <- laplace.approx$Sigmastar
    Wtilde <- laplace.approx$Wtilde
    Blat <- as.numeric(B %*% latentstar)
    wbar <- as.numeric(t(B) %*% (invgam.scale * (y - mu.gam(latentstar)) +
                                   Wtilde %*% Blat))
    L <- t(chol((t(B) %*% (diag(Wtilde) * B)) + Q))
    Mv <- chol2inv(t(L))
    Pvj <- function(j){
      as.matrix(Matrix::bdiag(matrix(0, nrow = p, ncol = p),
                              diag(diag(q)[, j] * exp(v), q) %x% P))
    }
    Bmw <- as.numeric(B %*% Mv %*% wbar)

    partial.diag <- function(j) {
      Gamma.vj <- Mv %*% Pvj(j) %*% Mv
      MPMP <- Mv %*% Pvj(j) %*% Mv %*% Pvj(j)
      MPMP.min <- 2 * (MPMP %*% Mv) - Gamma.vj
      Gwbar <- Gamma.vj %*% wbar
      Mwbar <- Mv %*% wbar

      val <-  .5 * sum(diag((Mv %*% Pvj(j) %*% Mv %*% Pvj(j) - Mv %*% Pvj(j)))) +
        invgam.scale * sum(y * (B %*% MPMP.min %*% wbar)) -
        invgam.scale * sum((s.gam.deriv(Bmw) *
                              as.numeric(B %*% MPMP.min %*% wbar)) + (s.gam.deriv2(Bmw) *
                                                                        as.numeric(B %*% Gwbar) ^ 2)) +
        as.numeric(-2 * (t(wbar) %*% MPMP %*% Mv %*% Q %*% Mwbar) +
                     (t(wbar) %*% Gamma.vj %*% (Q + Pvj(j)) %*% Mwbar) -
                     (t(wbar) %*% Gamma.vj %*% Q %*% Gwbar)) +
        as.numeric(t(wbar) %*% MPMP %*% Mwbar) -
        as.numeric(.5 * (t(wbar) %*% Gwbar)) -
        (b.delta * (1 + (2 * a.delta) / nu) * exp(-v[j])) /
        ((1 + (2 * b.delta) / (nu * exp(v[j]))) ^ 2)
      return(val)
    }
    deriv.diag <- sapply(seq_len(q), partial.diag)

    if(q > 1){
      partial.offdiag <- function(s, j) {
        Gamma.vj <- Mv %*% Pvj(j) %*% Mv
        Gamma.vs <- Mv %*% Pvj(s) %*% Mv
        MPsMPj <- Mv %*% Pvj(s) %*% Mv %*% Pvj(j)
        GPM <- (Gamma.vs %*% Pvj(j) %*% Mv) + (Mv %*% Pvj(j) %*% Gamma.vs)
        Mwbar <- Mv %*% wbar

        val <- .5 * sum(diag(MPsMPj)) +
          invgam.scale * sum(y * (B %*% GPM %*% wbar)) -
          invgam.scale * sum((s.gam.deriv(Bmw) * (B %*% GPM %*% wbar)) +
                               (s.gam.deriv2(Bmw) * (B %*% Gamma.vs %*% wbar) *
                                  (B %*% Gamma.vj %*% wbar))) -
          as.numeric(t(wbar) %*% Gamma.vs %*% Pvj(j) %*% Mv %*% Q %*% Mwbar) -
          as.numeric(t(wbar) %*% Mv %*% Pvj(j) %*% Gamma.vs %*% Q %*% Mwbar) +
          as.numeric(t(wbar) %*% Gamma.vj %*% Pvj(s) %*% Mwbar) -
          as.numeric(t(wbar) %*% Gamma.vj %*% Q %*% Gamma.vs %*% wbar) +
          as.numeric(t(wbar) %*% Gamma.vs %*% Pvj(j) %*% Mwbar)
        return(val)
      }

      deriv.offdiag <- matrix(0, nrow = q, ncol = q)

      col.idx <- 2
      for (i in 1:(q - 1)) {
        for (j in col.idx:q) {
          deriv.offdiag[j, i] <- partial.offdiag(i, j)
          deriv.offdiag[i, j] <- partial.offdiag(i, j)
        }
        col.idx <- col.idx + 1
      }

      hessmatrix <- adjustPD(-(deriv.offdiag + diag(deriv.diag)))$PD
      hess.mat <- (-1) * hessmatrix
    } else{hess.mat <- deriv.diag}

    return(hess.mat)

  }

  #############################################################################
  ###       Newton-Raphson with step-halving to find max of logp.vD         ###
  #############################################################################

  NRaphson <- function(v0, epsilon = 1e-05, maxiter = 100){

    NRoutput <- matrix(0, nrow = maxiter + 1, ncol = q + 3)
    colnames(NRoutput) <- c("Iteration", paste0(rep("v", q), seq_len(q)), "f",
                            "Distance")
    NRoutput[1, 2 : (q + 1)] <- v0
    L.approx <- Laplace(rep(0, H), v0)
    NRoutput[1, (q + 2)] <- logp.vD0(v0, L.approx)$value
    NRoutput[1, (q + 3)] <- NA

    iter <- 0
    v <- v0

    for (k in 1:maxiter) {
      gradf <- gradient(v, L.approx) # Gradient of function at v
      Hessf <- Hessian(v, L.approx)  # Hessian of function at v
      dv <- (-1) * solve(Hessf, gradf) # Ascent direction
      dv <- (dv / sqrt(sum(dv ^ 2))) * 5 # Step damping
      vnew <- v + dv # New point
      step <- 1 # Initial step length
      iter.halving <- 1 # Counter for step halving
      fv <- logp.vD0(v, L.approx)$value       # Function value at v
      fvnew <- logp.vD0(vnew, L.approx)$value # Function value at vnew
      # Step halving
      while (fvnew <= fv) {
        step <- step * .5 # Decrease step size
        vnew <- v + (step * dv)
        fvnew <- logp.vD0(vnew, L.approx)$value
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) break
      }
      dist <- sqrt(sum((vnew - v) ^ 2))
      iter <- iter + 1
      v <- vnew
      L.approx <- Laplace(L.approx$latentstar, v = v)
      NRoutput[k + 1, 1] <- iter
      NRoutput[k + 1, 2 : (q + 1)] <- vnew
      NRoutput[k + 1, (q + 2)] <- fvnew
      NRoutput[k + 1, (q + 3)] <- dist

      if(dist < epsilon){ # stop algorithm
        listout <- list(voptim = vnew,
                        info = NRoutput[1 : (iter + 1), ],
                        converge = TRUE,
                        iter = iter,
                        foptim = fvnew,
                        distance = dist,
                        L.approx = L.approx)
        attr(listout, "class") <- "NewtonRaphson"
        return(listout)
        break
      }
    }

    if(dist >= epsilon){
      listout <- list(info = NRoutput,
                      converge = FALSE,
                      iter = iter)
      attr(listout, "class") <- "NewtonRaphson"
      return(listout)
    }

  }

  print.NewtonRaphson <- function(obj){
    if(obj$converge == TRUE){
      cat("Newton-Raphson algorithm converged after", obj$iter,
          "iterations.", "\n")
      cat("\n")
      cat("Maximum point:", format(round(obj$voptim, 4), nsmall = 4), "\n")
      cat("Value of objective function at maximum:",
          format(round(obj$foptim,4), nsmall = 4), "\n")
      cat("\n")
      print(as.data.frame(round(obj$info, 6)))
    }
    else(cat("Newton-Raphson algorithm failed convergence"))
  }

  NR.start <- rep(6, q)  # Newton-Raphson starting point
  NR.counter <- 1        # Number of Newton-Raphson runs
  NR.fail <- 0           # Indicator if Newton-Raphson fails

  newton <- NRaphson(NR.start)

  while (sum(abs(gradient(newton$voptim,
                          Laplace(rep(0, H),newton$voptim))) < 0.5) != q) {
    NR.start <- NR.start - 2
    newton <- NRaphson(NR.start)
    NR.counter <- NR.counter + 1
    if (NR.counter >= 7) {
      NR.fail <- 1
      break
    }
  }
  if(NR.fail == 1)
    stop("Newton-Raphson did not converge")

  # Extracting info from newton
  L.approx <- newton$L.approx
  v.max <- newton$voptim
  hess.max <- Hessian(v.max, L.approx)
  logpvmax <- newton$foptim
  Cov.vmax <- solve(-hess.max)
  var.max <- diag(Cov.vmax)

  # Compute the covariance of theta_j at a v.max for all functions j=1,..,q
  post.max <-  logp.vD0(v.max, L.approx)
  latmaximum <- post.max$latentstar
  Covmaximum <- post.max$Sigmastar
  Wtilde.max <- L.approx$Wtilde
  BWB <- t(B) %*% Wtilde.max %*% B

  # Diagonal of hat matrix (as a function of v)
  diagHmatrix <- function(v) diag(solve(BWB + Qv(v)) %*% BWB)

  # Estimated effective degrees of freedom for each latent field component
  edf.full <- diagHmatrix(v.max)
  edf.smooths <- matrix(edf.full[-(1:p)], ncol = (K - 1), nrow = q, byrow = TRUE)
  edf <- c(edf.full[1:p], as.vector(t(edf.smooths)))
  names(edf) <- c(colnames(Z),
                  paste0(rep(colnames(X)[(p + 1):(p + q)], rep(K - 1, q)), ".",
                         seq_len(K - 1)))
  edf[1:p] <- round(edf[1:p], 3)

  # Estimated effective degrees of freedom of the smooth terms
  ED.functions <- c()
  for (j in 1:q) {
    fjdf <- edf[((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1))]
    edf.j <- sum(fjdf)
    if(edf.j < (penorder - 1)) edf.j <- sum(fjdf[fjdf > 0])
    if(edf.j < (penorder - 1)) edf.j <- (penorder - 1)
    ED.functions[j] <- edf.j
  }

  ED.global <- sum(edf[1:p]) + sum(ED.functions)

  # Log-posterior of log-penalty vector at Wtilde.max, wbar.max
  logp.vD <- function(v) {
    Q <- Qv(v)
    latentstar <- latmaximum
    Sigmastar <- Covmaximum
    Wtilde <- Wtilde.max
    L <- t(chol((t(B) %*% (diag(Wtilde) * B)) + Q))
    Mv <- chol2inv(t(L))
    Blat <- as.numeric(B %*% latentstar)

    lpost <- (-1) * sum(log(diag(L))) +
      .5 * (nu + K - 1) * sum(v) +
      (1 / gam.scale) * sum(y * Blat) -
      (1 / gam.scale) * sum(s.gam(latentstar)) -
      .5 * as.numeric(t(latentstar) %*% Q %*% latentstar) -
      (.5 * nu + a.delta) * sum(log(b.delta + .5 * nu * exp(v)))

    listout <- list(value = lpost,
                    latentstar = latentstar,
                    Sigmastar = Sigmastar)

    return(listout)
  }

  # Matrix to host sampled values of v for effective dimension HPD interval
  nv.sample <- 500 # Number of samples to draw for vector v
  v.sample <- matrix(0, nrow = nv.sample, ncol = q)

  #############################################################################
  ##########           Grid-based inference when q < 5                      ###
  #############################################################################
  pendist.params <- matrix(nrow = q, ncol = 3)

  if(q < 5){
    crit.val <- exp(-.5 * stats::qchisq(0.9, df = q, lower.tail = TRUE))
    mode.check <- 1

    mix.components <- function(){

      # Univariate grid construction
      if(q == 1){grid.size <- 10}
      if(q == 2){grid.size <- 8}
      if(q == 3){grid.size <- 5}
      if(q == 4){grid.size <- 5}
      univariate.grids <- list() # list in which to store univariate grids
      dvj <- c() # univariate grid width

      for(j in 1:q){
        # step size
        left.step.j  <- (-1) * diag(1, q)[j, ]
        right.step.j <- diag(1, q)[j, ]

        # Explore left direction in dimension j
        v.left.j <- v.max
        ratio <- 1
        while(ratio > crit.val){
          v.left.j <- v.left.j + left.step.j
          logpvleft <- logp.vD(v.left.j)$value
          if(logpvleft > logpvmax){
            mode.check <- 0
            break
          }
          ratio <- exp(logpvleft-logpvmax)
        }
        if(mode.check == 0){break}

        # Explore right direction in dimension j
        v.right.j <- v.max
        ratio <- 1
        while(ratio > crit.val){
          v.right.j <- v.right.j + right.step.j
          logpvright <- logp.vD(v.right.j)$value
          if(logpvright > logpvmax){
            mode.check <- 0
            break
          }
          ratio <- exp(logpvright-logpvmax)
        }
        if(mode.check == 0){break}

        lb.j <- v.left.j[j]  # lower bound in dimension j
        ub.j <- v.right.j[j] # upper bound in dimension j
        x.j <- seq(lb.j, ub.j, length = 20) # equidistant grid
        dj <- x.j[2]-x.j[1] # grid width
        vgrid.j <- matrix(rep(v.max, 20), ncol = q, byrow = TRUE)
        vgrid.j[, j] <- x.j
        y.j <- unlist(lapply(apply(vgrid.j, 1, logp.vD),'[[',1))
        y.j <- exp(y.j - max(y.j))
        cnorm <- 1 / sum(y.j * dj)
        y.j <- cnorm * y.j
        snfit.j <- snmatch(x.j, y.j)
        v.sample[, j] <- sn::rsn(nv.sample, xi = snfit.j$location,
                             omega = snfit.j$scale,
                             alpha = snfit.j$shape)
        pendist.params[j, ] <- c(snfit.j$location, snfit.j$scale, snfit.j$shape)
        colnames(pendist.params) <- c("location", "scale", "shape")
        rownames(pendist.params) <- paste0(rep("SN.", q), "v", seq_len(q))
        q0.025 <- snfit.j$quant[1] # Lower bound of grid
        q0.975 <- snfit.j$quant[3] # Upper bound of grid
        grid.vals <- seq(q0.025, q0.975, length = grid.size) # Starting grid

        # Shift grid so that it goes through v.max[j]
        abs.dist <- abs(grid.vals - v.max[j])
        min.dist <- min(abs.dist)
        pos.close <- which.min(abs.dist)
        if(grid.vals[pos.close] - v.max[j] >= 0){
          grid.vals <- grid.vals - min.dist
          dir <- "left"}else{
            grid.vals <- grid.vals + min.dist
            dir <- "right"}
        if(dir == "left"){
          while(utils::tail(grid.vals, 1) < q0.975){
            grid.vals <- c(grid.vals, utils::tail(grid.vals, 1) + diff(grid.vals)[1])
          }
        }
        if(dir == "right"){
          while(utils::head(grid.vals, 1) > q0.025){
            grid.vals <- c(utils::head(grid.vals, 1) - diff(grid.vals)[1], grid.vals)
          }
        }
        univariate.grids[[j]] <- grid.vals
        dvj[j] <- grid.vals[2]-grid.vals[1]
      }

      if(mode.check == 1){ #continue and construct grid
        # Construction of the mesh (cartesian product of univariate grids)
        mesh.cartesian <- expand.grid(univariate.grids)
        eval.mesh.cartesian <- apply(mesh.cartesian, 1, logp.vD)
        eval.target <- exp(unlist(lapply(eval.mesh.cartesian, '[[', 1)) -
                             logpvmax)
        low.contrib <- which(eval.target < crit.val) # Low contribution points

        ## Computation of the mixture components
        if(length(low.contrib) > 0){ # Filter grid
          M <- nrow(mesh.cartesian)-length(low.contrib) # No. of points after filter
          mesh.log.image <- log(eval.target[-low.contrib])
          xistar.mat <- matrix(unlist(lapply(
            eval.mesh.cartesian,'[[', 2)[-low.contrib]),
            ncol = H, nrow = M, byrow = TRUE)
          var.regcoeffs <- matrix(unlist(lapply(lapply(
            eval.mesh.cartesian, '[[', 3)[-low.contrib], diag)),
            ncol = H, nrow = M, byrow = TRUE)[, 1:p]
        }else{ # No need to filter the grid
          M <- nrow(mesh.cartesian)
          mesh.log.image <- log(eval.target)
          xistar.mat <- matrix(unlist(lapply(
            eval.mesh.cartesian, '[[', 2)),
            ncol = H, nrow = M, byrow = TRUE)
          var.regcoeffs <- matrix(unlist(lapply(lapply(
            eval.mesh.cartesian,'[[',3), diag)),
            ncol = H, nrow = M, byrow = TRUE)[, 1:p]
        }
        omega.m <- exp(mesh.log.image)/sum(exp(mesh.log.image))

        list.out <- list(nquad = M, # Number of quadrature points
                        xistar.mat = xistar.mat, # Matrix of latent field values
                        var.regcoeffs = var.regcoeffs, # Variance of reg. coeffs)
                        omega.m = omega.m,     # Mixture weights
                        v.sample = v.sample,   # Matrix of samples for v
                        pendist.params = pendist.params, # Parameters of sn fit
                        mode.check = mode.check) # Mode check
        return(list.out)
      } else{ # avoid grid
        return(list(mode.check = mode.check))
      }
    }

    mixture.comp <- mix.components()

    if(mixture.comp$mode.check == 1){
      pen.family <- "skew-normal"
      pendist.params <- mixture.comp$pendist.params
      v.sample <- mixture.comp$v.sample
      M <- mixture.comp$nquad
      xistar.mat <- mixture.comp$xistar.mat
      omega.m <- mixture.comp$omega.m
      var.regcoeffs <- as.matrix(mixture.comp$var.regcoeffs)
      latent.hat <- colSums(omega.m * xistar.mat)
      betas.hat <- latent.hat[1:p]
      sdbeta.post <- sqrt(colSums(omega.m * (var.regcoeffs + (xistar.mat[, 1:p] -
                                  matrix(rep(betas.hat, M), nrow = M, ncol = p,
                                  byrow = TRUE)) ^ 2)))
      zscore <- betas.hat / sdbeta.post # Bayesian z-score

      # Credible intervals for regression coefficients
      betas.matCI <- matrix(0, nrow = p, ncol = 2)

      for(j in 1:p){
        support.j <- seq(betas.hat[j] - 4.5 * sdbeta.post[j],
                         betas.hat[j] + 4.5 * sdbeta.post[j],
                         length = 400)
        support.mat <- matrix(0, ncol = length(support.j), nrow = M)
        for(i in 1:length(support.j)){
          for(m in 1:M){
            support.mat[m, i] <- stats::dnorm(support.j[i], mean = xistar.mat[m,j],
                                       sd = sqrt(var.regcoeffs[m, j]))}
        }
        mix.image <- colSums(omega.m * support.mat)
        cumsum.mix <- cumsum(mix.image * diff(support.j)[1])
        lb.j <- support.j[sum(cumsum.mix < ((1 - cred.int) *.5))]
        ub.j <- support.j[sum(cumsum.mix < (1 - (1 - cred.int) * .5))]
        betas.matCI[j, ] <- c(lb.j, ub.j)
      }
    }
  }

  #############################################################################
  ###   When q >=5 inference is based on maximum a posteriori of v          ###
  #############################################################################

  if(q >= 5 || (q < 5 && mixture.comp$mode.check == 0)){
    pendist.params <- v.max
    pen.family <- "gaussian"
    latent.hat <- latmaximum
    betas.hat <- latent.hat[1:p]
    sdbeta.post <- sqrt(diag(Covmaximum)[1:p])
    zscore <- betas.hat / sdbeta.post # Bayesian z-score
    v.sample <- MASS::mvrnorm(n = nv.sample, mu = v.max, Sigma = Cov.vmax)

    # Credible intervals for regression coefficients
    betas.matCI <- matrix(0, nrow = p, ncol = 2)

    for(j in 1:p){

      betas.matCI[j, 1] <- betas.hat[j] - stats::qnorm(1-(1-cred.int) * .5) *
        sdbeta.post[j]
      betas.matCI[j, 2] <- betas.hat[j] + stats::qnorm(1-(1-cred.int) * .5) *
        sdbeta.post[j]
    }
  }

  #############################################################################
  ###                        Output results                                 ###
  #############################################################################

  # 95% HPD interval for effective dimension of smooth terms
  ED.functions.sample <- matrix(0, nrow = nv.sample, ncol = q)
  matdiagH <- matrix(0, nrow = nv.sample, ncol = H)

  for (s in 1:nv.sample) {
    if(class(try(diagHmatrix(v.sample[s, ]), silent = TRUE)) == "numeric") {
      matdiagH[s, ] <- diagHmatrix(v.sample[s, ])
    } else {
      matdiagH[s, ] <- diagHmatrix(v.max)
    }
    edf.v <- matdiagH[s, ]
    for (j in 1:q) {
      edf.vj <- edf.v[((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1))]
      edf.vjsum <- sum(edf.vj)
      if(edf.vjsum < (penorder - 1)) edf.vjsum <- sum(edf.vj[edf.vj > 0])
      if(edf.vjsum < (penorder - 1)) edf.vjsum <- penorder - 1
      ED.functions.sample[s, j] <- edf.vjsum
    }
  }
  ED.functions.sample <- coda::as.mcmc(ED.functions.sample)
  ED.functions.HPD95 <- coda::HPDinterval(ED.functions.sample)[1:q, ]
  if (q == 1) ED.functions.HPD95 <- matrix(ED.functions.HPD95, ncol = 2)
  colnames(ED.functions.HPD95) <- c("lower .95", "upper .95")

  # Fitted values of gam
  fitted.values <- s.gam.deriv(as.numeric(B %*% latent.hat))

  # Response residuals
  residuals <- y - fitted.values

  # Adjusted R-squared
  r2.adj <- 1 - ((sum(residuals ^ 2) / (n - p)) / stats::var(y))

  Fmat <- solve(BWB + Qv(v.max)) %*% BWB
  edf2 <- diag(2 * Fmat - (Fmat %*% Fmat))

  ED.functions2 <- c()
  for (j in 1:q) {
    ED.functions2[j] <- sum(edf2[((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1))])
  }

  # Approximate significance of smooth terms: test H0: fj=0
  approx.signif <- function(j){

    Vthetaj <- Covmaximum[((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1)),
                          ((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1))]

    Vfj <- B.list.trim[[j]] %*% Vthetaj %*% t(B.list.trim[[j]])
    fhatj <- as.numeric(B.list.trim[[j]] %*% thetahat.list[[j]])

    k <- floor(ED.functions2[j]) + 1
    nnu <- ED.functions2[j] - k + 1
    rrho <- sqrt(nnu *  (1 - nnu) * .5)

    eigdecomp <- RSpectra::eigs(Vfj, k = k, which = "LM")
    eigk <- eigdecomp$values
    U <- eigdecomp$vectors

    Lamb.tilde <- diag(utils::tail(eigk, 2) ^ (-.5))
    Brhonu <- matrix(c(1, rrho, rrho, nnu), ncol = 2, byrow = TRUE)
    Bpseudo <- Lamb.tilde %*% Brhonu %*% t(Lamb.tilde)

    if(k > 2) {
      Lambfull.tilde <- diag(0, k)
      Lambfull.tilde[1:(k - 2), 1:(k - 2)] <- diag(eigk[1:(k - 2)] ^ (-1), k - 2)
      Lambfull.tilde[(k - 1):k, (k - 1):k] <- Bpseudo
      Vfjinv <- U %*% Lambfull.tilde %*% t(U)
      Tr <- as.numeric(t(fhatj) %*% Vfjinv %*% fhatj)
    } else {
      Lambfull.tilde <- Bpseudo
      Vfjinv <- U %*% Lambfull.tilde %*% t(U)
      Tr <- as.numeric(t(fhatj) %*% Vfjinv %*% fhatj)
    }
    pvalue <- stats::pgamma(Tr, shape = ED.functions2[j] * .5, scale= 2,
                            lower.tail = FALSE)
    return(list(Tr = Tr, pval = pvalue, dim = j))
  }


  # List of estimated B-spline coefficients
  thetahat.list <- as.list(as.data.frame(t(matrix(latent.hat[(p + 1):H],
                                                  ncol = (K - 1),
                                                  nrow = q, byrow = TRUE))))
  names(thetahat.list) <- paste0(rep("theta", q), seq_len(q))

  Approx.signif <- matrix(unlist(lapply(seq_len(q), approx.signif)),
                          ncol = 3, byrow = TRUE)[, 1:2]

  # The data frame
  if(missing(data)){
    data = NULL
  } else{
    data = data
  }

  # Posterior estimates for linear regression coefficients
  betas_estim <- matrix(0, nrow = p, ncol = 5)
  betas_estim[, 1] <- betas.hat
  betas_estim[, 2] <- sdbeta.post
  betas_estim[, 3] <- zscore
  betas_estim[, 4] <- betas.matCI[, 1]
  betas_estim[, 5] <- betas.matCI[, 2]
  if(p > 1){
    rownames(betas_estim) <- colnames(X[, 1:p])
  }else{rownames(betas_estim) <- "(Intercept)"}
  colnames(betas_estim) <- c("Estimate", "sd.post", "z-score",
                             paste0("lower.", round(cred.int  * 100, 2)),
                             paste0("upper.", round(cred.int  * 100, 2)))

  # List of objects to return
  outlist <- list(formula = formula,          # Model formula
                  family = family,            # Gam family
                  link = link,                # Link used
                  n = n,                      # Sample size
                  q = q,                      # Number of smooth terms
                  K = K,                      # Number of B-splines
                  penalty.order = penorder,   # Chosen penalty order
                  latfield.dim = H,           # Latent field dimension
                  linear.coeff = betas_estim, # Estimated linear coeff.
                  spline.estim = thetahat.list,  # List of B-spline coeff.
                  edf = edf, # Degrees of freedom for each latent field variable
                  Approx.signif = Approx.signif, # Approximate significance smooth
                  EDf = ED.functions, # Effective degrees of freedom of smooths
                  EDfHPD.95 = ED.functions.HPD95, # 95% HPD for EDf
                  ED = ED.global,     # Effective degrees of freedom of model
                  vmap = v.max, # Maximum a posteriori of log(penalty) vector
                  Cov.vmap = Cov.vmax, # Covariance of log(penalty) vector at vmap
                  pen.family = pen.family, # Posterior distribution for v
                  pendist.params = pendist.params, # Parameters of posterior dist of v
                  Covmaximum = Covmaximum, # Covariance of lat field at v.max
                  latmaximum = latmaximum, # latent field at v.max
                  fitted.values = fitted.values, # Gam fitted values
                  residuals = residuals,         # Response residuals
                  r2.adj = r2.adj,               # Adjusted R-squared
                  data = data)                   # The data frame

  attr(outlist, "class") <- "gamlps"
  outlist

}








