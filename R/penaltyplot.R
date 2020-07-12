#' Plot the approximate posterior distribution of the penalty vector.
#'
#' @description The routine gives a graphical representation of the univariate
#' approximate posterior distribution of the (log-)penalty parameters for
#' objects of class \emph{coxlps}, \emph{curelps}, \emph{amlps} and
#' \emph{gamlps}.
#'
#' @param object An object of class \code{coxlps}, \code{curelps}, \code{amlps}
#'  or \code{gamlps}.
#' @param dimension For objects of class \code{amlps} and \code{gamlps}, the
#'  penalty vector can have a dimension larger than one, i.e. more than a
#'  single smooth term is present in the considered additive model. In that case,
#'  \code{dimension} is the penalty dimension to be plotted
#'  corresponding either to a scalar indicating the desired dimension or to a
#'  vector indicating more than one dimension. For instance, dimension = c(1,3)
#'  displays two separate plots of the (approximate) posterior distribution of
#'  the (log-)penalty parameter associated to the first and the third smooth
#'   function respectively.
#' @param ... Further arguments to be passed to the routine.
#'
#' @details When q, the number of smooth term in a (generalized) additive model is
#'  smaller than five, the exploration of the posterior penalty space is based
#'  on a grid strategy. In particular, the multivariate grid of dimension q is
#'  constructed by taking the Cartesian product of univariate grids in each
#'  dimension j = 1,...q. These univariate grids are obtained from a skew-normal
#'  fit to the conditional posterior p(vj|vmap[-j]),D), where vj is the
#'  (log-)penalty value associated to the jth smooth function and vmap[-j] is
#'  the posterior maximum of the (log-)penalty vector omitting the jth
#'  dimension. The routine displays the latter skew-normal distributions. When
#'  q>=5, inference is based on vmap and the grid is omitted to avoid
#'  computational overflow. In that case, the posterior distribution of the
#'  (log-)posterior penalty vector v is approximated by a multivariate Gaussian
#'   and the routine shows the marginal distributions.
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @examples
#' ### Classic simulated data example (with simgamdata)
#'
#' set.seed(123)
#' sim.data <- simgamdata(setting = 2, n = 250, dist = "gaussian", scale = 0.25)
#' plot(sim.data)         # Scatter plot of response
#' data <- sim.data$data  # Simulated data frame
#' # Fit model
#' fit <- amlps(y ~ z1 + z2 + sm(x1) + sm(x2), data = data, K = 15)
#' fit
#'
#' # Penalty plot
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 2))
#' penaltyplot(fit, dimension = c(1, 2))
#' par(opar)
#'
#' @export

penaltyplot <- function(object, dimension, ...){

  # Print penalty plot for a coxlps object
  if(class(object) == "coxlps"){
    # Extract required objects from the fit
    n <- object$n
    tlow <- 0
    tup <- object$tup
    time <- object$event.times
    event <- object$event.indicators
    K <- object$K
    X <- scale(object$X, center = TRUE, scale = FALSE)
    penorder <- object$penalty.order
    H <- object$latfield.dim
    p <- object$p
    vquad <- log(object$penalty.vector)
    latent0 <- as.numeric(c(object$regcoeff[, 1], object$spline.estim))
    vmap <- object$vmap

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
    for (k in 1:penorder) D <- diff(D)
    P <- t(D) %*% D           # Penalty matrix of dimension c(K,K)
    P <- P + diag(1e-06, K)   # Perturbation to ensure P full rank

    # Robust prior following Jullion and Lambert (2007)
    a.delta <- 1e-04
    b.delta <- 1e-04
    nu <- 3
    prec.betas <- 1e-05 # Precision for regression coefficients

    # Log-likelihood function
    loglik <- function(latent) {
      val <- sum(event * (Bdesign %*% latent)) - sum(cumsum(exp((
        Bmiddle %*% latent[1:K])) * width)[upto] * exp(X %*% latent[(K + 1):H]))
      return(val)
    }

    # Latent field prior precision matrix as a function of v = log(lambda)
    Q.v <- function(v) {
      Qmatrix <- matrix(0, nrow = H, ncol = H)
      Qmatrix[1:K, 1:K] <- exp(v) * P
      Qmatrix[(K + 1):H, (K + 1):H] <- diag(prec.betas, p)
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
                                             apply(prodwidth, 2, cumsum)[upto,])
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
      upleft.block <-
        (-1) * (t(Bmiddle) %*% (expBmidtheta * colvecsum * Bmiddle))

      # Upper right block
      matprod.upright <- expBmidtheta * Bmiddle * width
      upright.block <- (-1) * (t(apply(matprod.upright, 2, cumsum)[upto, ] *
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

      y <- loglik(latentstar) - (.5 * sum((Qv * latentstar) %*% latentstar)) +
        (.5 * (K + nu) * v) - sum(log(diag(L))) -
        (a.delta + .5 * nu) * log(b.delta + .5 * nu * exp(v))
      Covmatrix <- solve(Qv - Htilde)
      return(y)
    }

    # Plot of approximate posterior log-penalty distribution
    np <- 100
    vlim <- c(vquad[1], vquad[length(vquad)])
    logpvmax <- logpost.v(vmap, latent0)
    vgrid <- seq(vlim[1], vlim[2], length = np)
    dvgrid <- vgrid[2] - vgrid[1]
    y <- sapply(vgrid, logpost.v, latent0 = latent0)
    y <- exp(y - logpvmax)
    cnorm <- 1/ sum(y * dvgrid)
    y <- cnorm * y
    graphics::plot(vgrid, y, type = "l", xlab = paste0("v"), ylab = "p(v|D)",...)

  }

  # Print penalty plot for a curelps object
  else if(class(object) == "curelps"){
    # Extract required objects from the fit
    n <- object$n
    tlow <- 0
    tup <- object$tup
    time <- object$event.times
    event <- object$event.indicators
    K <- object$K
    X <- scale(object$X, center = TRUE, scale = FALSE)
    X[, 1] <- rep(1, n)
    Z <- scale(object$Z, center = TRUE, scale = FALSE)
    penorder <- object$penalty.order
    H <- object$latfield.dim
    p <- object$p
    latent0 <- as.numeric(c(object$spline.estim, object$coeff.probacure[, 1],
                            object$coeff.cox[, 1]))
    nbeta <- ncol(X)
    ngamma <- ncol(Z)
    vmap <- object$vmap
    vquad <- object$vquad
    constr <- 6

    # Bins and midpoints to estimate baseline survival
    bins <- 300                                      # Total number of bins
    partition <- seq(tlow, tup, length = bins + 1)   # Partition domain
    width <- partition[2] - partition[1]             # Bin width
    middleBins <- partition[1:bins] + (width / 2)    # Middle points bins
    upto <- as.integer(time / width) + 1  # Indicates in which bin time falls
    upto[which(upto == bins + 1)] <- bins # If time == tup (bin 301 to 300)

    # Functions used to construct Hessian blocks
    fill.ones <- function(vecidx) c(rep(1, vecidx[1]), rep(0, vecidx[2]))
    m.idx <- matrix(cbind(upto, (bins - upto)), ncol = 2)
    matones <- matrix(apply(m.idx, 1, fill.ones), nrow = n, ncol = bins,
                      byrow = TRUE)

    # Function to compute outer product between rows of two matrices
    outerprod <- function(i, X, Y) outer(X[i,], Y[i,])
    matxx <- lapply(seq_len(n), FUN = outerprod, X, X) # For Hessian blocks
    matxz <- lapply(seq_len(n), FUN = outerprod, X, Z) # For Hessian blocks
    matzz <- lapply(seq_len(n), FUN = outerprod, Z, Z) # For Hessian blocks
    matxx.fun <- function(factor) {
      mat <- matrix(0, nrow = nbeta, ncol = nbeta)
      for (i in 1:n) {
        mat <- mat + factor[i] * matxx[[i]]
      }
      return(unname(mat))
    }
    matxz.fun <- function(factor) {
      mat <- matrix(0, nrow = nbeta, ncol = ngamma)
      for (i in 1:n) {
        mat <- mat + factor[i] * matxz[[i]]
      }
      return(unname(mat))
    }
    matzz.fun <- function(factor) {
      mat <- matrix(0, nrow = ngamma, ncol = ngamma)
      for (i in 1:n) {
        mat <- mat + factor[i] * matzz[[i]]
      }
      return(unname(mat))
    }

    # B-spline basis
    Bobs <- cubicbs(time, tlow, tup, K)$Bmatrix
    Bmiddle <- cubicbs(middleBins, tlow, tup, K)$Bmatrix

    # Penalty matrix
    D <- diag(K)
    for (k in 1:penorder) D <- diff(D)
    P <- t(D) %*% D            # Penalty matrix
    P <- P + diag(1e-06, K)    # Perturbation to ensure P full rank

    # Robust prior following Jullion and Lambert (2007)
    a.delta <- 1e-04
    b.delta <- 1e-04
    nu <- 3
    prec.regcoeff <- 1e-05 # Precision for regression coefficients

    # Baseline survival based on observed survival times
    S0.time <- function(theta) {
      S0 <- exp(-cumsum(as.numeric(exp(Bmiddle %*% theta) * width))[upto])
      S0[which(S0 < 1e-10)] <- 1e-10
      return(S0)
    }

    # Baseline hazard based on observed survival times
    h0.time <- function(theta) {
      h0 <- as.numeric(exp(Bobs %*% theta))
      return(h0)
    }

    #Prior precision matrix on the latent field (as a function of v = log(lambda))
    Qv <- function(v) {
      Q <- matrix(0, nrow = H, ncol = H)
      Q[1:K, 1:K] <- exp(v) * P
      Q[(K + 1):H, (K + 1):H] <- diag(prec.regcoeff, (nbeta + ngamma))
      return(Q)
    }

    # Log-likelihood of promotion time cure model
    loglik <- function(latent) {
      # Split latent field
      theta <- latent[1:K]
      beta  <- latent[(K + 1):(K + nbeta)]
      gamma <- latent[(K + nbeta + 1):H]
      # Compute log-likelihood
      Xbeta <- as.numeric(X %*% beta)
      Zgamma <- as.numeric(Z %*% gamma)
      S0time <- S0.time(theta)

      loglikval <- sum(event * (Xbeta + Zgamma + log(h0.time(theta)) -
                                  exp(Zgamma) * (-log(S0time))) -
                         exp(Xbeta) * (1 - (S0time) ^ (exp(Zgamma))))
      return(loglikval)
    }

    # Gradient of the log-likelihood
    grad.loglik <- function(latent) {
      # Split latent field
      theta <- latent[1:K]
      beta  <- latent[(K + 1):(K + nbeta)]
      gamma <- latent[(K + nbeta + 1):H]
      Xbeta <- as.numeric(X %*% beta)
      Zgamma <- as.numeric(Z %*% gamma)
      Bmidtheta <- as.numeric(Bmiddle %*% theta)
      S0time <- S0.time(theta)
      omega.0k <- apply((exp(Bmidtheta) * Bmiddle * width), 2, cumsum)[upto,]

      # Derivative w.r.t the spline coefficients
      grad.theta <- colSums(event * (Bobs - exp(Zgamma) * omega.0k) -
                              exp(Xbeta + Zgamma) * (S0time ^ exp(Zgamma)) * omega.0k)

      # Derivative w.r.t the beta coefficients
      grad.beta <- colSums((event - exp(Xbeta) *
                              (1 - (S0time ^ exp(Zgamma)))) * X)

      # Derivative w.r.t the gamma coefficients
      grad.gamma <- colSums(event * Z * (1 - exp(Zgamma) * (-log(S0time))) -
                              exp(Xbeta + Zgamma) * (S0time ^ exp(Zgamma)) *  (-log(S0time)) * Z)

      grad.val <- as.numeric(c(grad.theta, grad.beta, grad.gamma))
      return(grad.val)
    }

    # Hessian of the log-likelihood
    Hessian.loglik <- function(latent) {
      # Split latent field
      theta <- latent[1:K]
      beta  <- latent[(K + 1):(K + nbeta)]
      gamma <- latent[(K + nbeta + 1):H]
      Xbeta <- as.numeric(X %*% beta)
      Zgamma <- as.numeric(Z %*% gamma)
      Bmidtheta <- as.numeric(Bmiddle %*% theta)

      omega.0 <- (-1) * log(S0.time(theta))
      omega.0k <- apply((exp(Bmidtheta) * Bmiddle * width), 2, cumsum)[upto,]
      omega.0kl <- function(vecn) {
        val <- t(Bmiddle) %*% diag(colSums(matones * vecn) *
                                     exp(Bmidtheta) * width) %*% Bmiddle
        return(val)
      }
      expXZ <- exp(Xbeta + Zgamma)
      expexp <- (exp(-omega.0) ^ (exp(Zgamma)))

      ### Block11
      factor <- expXZ * expexp
      Block11 <- omega.0kl((-1) * event * exp(Zgamma) - factor) +
        t(omega.0k * factor * exp(Zgamma)) %*% omega.0k

      ### Block12 and Block21
      Block12 <-  t((-1) * expXZ * expexp * omega.0k) %*% X
      Block21 <- t(Block12)

      ### Block13 and Block31
      Block13 <- t((-1) * event * exp(Zgamma) * omega.0k) %*% Z -
        t(expXZ * expexp * omega.0k) %*% (Z - omega.0 * exp(Zgamma) * Z)
      Block31 <- t(Block13)

      ### Block22
      factor <- (-1) * exp(Xbeta) * (1 - expexp)
      Block22 <- matxx.fun(factor)

      ### Block23 and Block32
      factor <- (-1) * expXZ * expexp * omega.0
      Block23 <- matxz.fun(factor)
      Block32 <- t(Block23)

      ### Block33
      factor1 <- omega.0 * expXZ * expexp
      factor2 <- (-1) * (event * exp(Zgamma) * omega.0 +
                           factor1 - factor1 * omega.0 * exp(Zgamma))
      Block33 <- matzz.fun(factor2)

      ### Inserting the blocks in Hessian matrix
      Hessianval <- matrix(0, nrow = H, ncol = H)
      Hessianval[1:K, 1:K] <- Block11
      Hessianval[1:K, ((K + 1):(K + nbeta))] <- Block12
      Hessianval[1:K, ((K + nbeta + 1):H)] <- Block13
      Hessianval[((K + 1):(K + nbeta)), 1:K] <- Block21
      Hessianval[((K + 1):(K + nbeta)), ((K + 1):(K + nbeta))] <- Block22
      Hessianval[((K + 1):(K + nbeta)), ((K + nbeta + 1):(H))] <- Block23
      Hessianval[((K + nbeta + 1):H), 1:K] <- Block31
      Hessianval[((K + nbeta + 1):H), ((K + 1):(K + nbeta))] <- Block32
      Hessianval[((K + nbeta + 1):H), ((K + 1 + nbeta):H)] <- Block33
      return(Hessianval)
    }

    # Log conditional posterior of the latent field given v
    logplat <- function(latent, v) {
      Q <- Qv(v)
      y <- as.numeric(loglik(latent) - .5 * sum((latent * Q) %*% latent))
      return(y)
    }

    # Gradient of the log conditional posterior of the latent field
    grad.logplat <- function(latent, v) {
      Q <- Qv(v)
      val <- grad.loglik(latent) - as.numeric(Q %*% latent)
      return(val)
    }

    # Hessian of the log conditional posterior of the latent field
    Hess.logplat <- function(latent, v) {
      Q <- Qv(v)
      Hess.mat <- Hessian.loglik(latent) - Q
      hessmatrix <- adjustPD(-Hess.mat)$PD
      hess.mat <- (-1) * hessmatrix
      return(hess.mat)
    }

    # Function that returns the logposterior of the log penalty value v and
    # other quantities related to the Laplace approximation
    logpost.v <- function(v, latent0) {
      # Laplace approximation to the conditional posterior of the latent field
      Laplace.approx <- function(latent0, v) {
        # Compute posterior mode via Newton-Raphson algorithm
        tol <- 1e-3            # Tolerance
        dist <- 3              # Initial distance
        iter <- 0              # Iteration counter
        maxiter <- 100         # maximum permitted iterations

        for (k in 1:maxiter) {
          gradf <- grad.logplat(latent0, v)         # Gradient of logplat at v
          Hessf <- Hess.logplat(latent0, v)         # Hessian of logplat at v
          dlat <- (-1) * solve(Hessf, gradf)        # Ascent direction
          dlat <- (dlat / sqrt(sum(dlat ^ 2))) * 4  # Step damping
          lat.new <- latent0 + dlat
          step <- 1
          iter.halving <- 1
          logplat.current <- logplat(latent0, v)
          logplat.new <- logplat(lat.new, v)
          while (logplat.new <= logplat.current) {
            step <- step * .5
            lat.new <- latent0 + (step * dlat)
            logplat.new <- logplat(lat.new, v)
            iter.halving <- iter.halving + 1
            if (iter.halving > 30) {
              break
            }
          }
          dist <- sqrt(sum((lat.new - latent0) ^ 2))
          iter <- iter + 1
          latent0 <- lat.new
          if (dist < tol) {
            break
          }
        }

        latentstar <- latent0   # Posterior mode of Laplace approximation given v
        Sigmastar <- solve((-1) * Hess.logplat(latentstar, v)) # Covariance matrix

        Sig11 <- Sigmastar[K, K]
        Sig21 <- Sigmastar[K,][-K]
        Sig12 <- t(Sig21)
        Sig22 <- Sigmastar[-K,-K]

        latentstar.c <- latentstar[-K] + Sig21 * (1 / Sig11) *
          (constr - latentstar[K])
        latentstar.cc <- c(latentstar.c[1:(K - 1)], constr,
                           latentstar.c[K:(H - 1)])
        latentstar.cc[(K + 1):H] <- latentstar[(K + 1):H]
        Sigmastar.c <- Sig22 - (1 / Sig11) * (Sig21 %*% Sig12)

        listout <- list(latstar = latentstar,
                        latstar.c = latentstar.c,
                        latstar.cc = latentstar.cc,
                        Sigmastar = Sigmastar,
                        Sigstar.c = Sigmastar.c)

        return(listout)
      }

      # Computation of the Laplace approximation
      Laplace <- Laplace.approx(latent0, v)

      # Output of Laplace approximation
      Sigmastar <- Laplace$Sigmastar # Unconstrained latent  covariance matrix
      Sigstar.c <- Laplace$Sigstar.c # Constrained covariance matrix
      latstar <- Laplace$latstar     # Unconstrained mode of Laplace approx.
      latstar.c <- Laplace$latstar.c # Constrained mode of dimension H-1
      latstar.cc <- Laplace$latstar.cc # Constrained mode of dimension H

      # Computation of the log-posterior of v
      Q.v <- Qv(v)
      LQ <- t(chol(Q.v))
      LSig <- t(chol(Sigstar.c))

      logpv <- loglik(latstar.cc) - .5 *  sum((latstar.cc * Q.v) %*% latstar.cc) +
        sum(log(diag(LQ))) + sum(log(diag(LSig))) + 5 * nu * v -
        (.5 * nu + a.delta) * log(b.delta + .5 * nu * exp(v))

      return(logpv)
    }

    # Plot of approximate posterior log-penalty distribution
    np <- 100
    vlim <- c(vquad[1], vquad[length(vquad)])
    logpvmax <- logpost.v(vmap, latent0)
    vgrid <- seq(vlim[1], vlim[2], length = np)
    dvgrid <- vgrid[2] - vgrid[1]
    y <- sapply(vgrid, logpost.v, latent0 = latent0)
    y <- exp(y - logpvmax)
    cnorm <- 1/ sum(y * dvgrid)
    y <- cnorm * y
    graphics::plot(vgrid, y, type = "l", xlab = paste0("v"), ylab = "p(v|D)",...)
  }

  # Print penalty plot for a amlps or gamlps object
  else if(class(object) == "amlps" || class(object) == "gamlps"){

    if (missing(dimension))
      stop("A dimension vector needs to be specified")
    if (!is.vector(dimension, mode = "numeric"))
      stop("Dimension must be a numeric vector")
    if (any(dimension < 1) || any(dimension > object$q))
      stop("Wrong dimension. The dimension vector should contain values
           between 1 and the total number of smooth terms in the model")

    if (object$pen.family == "skew-normal") {
      for (j in 1:length(dimension)) {
        lb <- sn::qsn(0.001, dp = object$pendist.params[dimension[j],])
        dsnlb <- sn::dsn(lb, dp = object$pendist.params[dimension[j],])
        while (dsnlb > 1e-4) {
          lb <- lb - .3
          dsnlb <- sn::dsn(lb, dp = object$pendist.params[dimension[j], ])
        }
        ub <- sn::qsn(0.999, dp = object$pendist.params[dimension[j],])
        dsnub <- sn::dsn(ub, dp = object$pendist.params[dimension[j],])
        while (dsnub > 1e-4) {
          ub <- ub + .3
          dsnub <- sn::dsn(ub, dp = object$pendist.params[dimension[j], ])
        }
        vdom <- seq(lb, ub, length = 500)
        yylab <- paste0("p(v", dimension[j], "|D)")


        graphics::plot(vdom, sn::dsn(vdom, dp = object$pendist.params[dimension[j],]),
             type = "l", xlab = paste0("v", dimension[j]), ylab = yylab, lwd = 2)
        location <- format(round(object$pendist.params[dimension[j], 1], 2),
                           nsmall = 2)
        scale <- format(round(object$pendist.params[dimension[j], 2], 2),
                        nsmall = 2)
        shape <- format(round(object$pendist.params[dimension[j], 3], 2),
                        nsmall = 2)
        graphics::legend("topleft", paste0("SN(", location, ",", scale, ",", shape, ")"),
               lty = 1, lwd = 2, bty = "n", cex = 0.8)
      }
    } else if (object$pen.family == "gaussian") {
      for (j in 1:length(dimension)) {
        mean.norm <- object$vmap[dimension[j]]
        sd.norm <- sqrt(diag(object$Cov.vmap)[dimension[j]])
        lb <- mean.norm - 4 * sd.norm
        ub <- mean.norm + 4 * sd.norm
        vdom <- seq(lb, ub, length = 500)

        graphics::plot(vdom, stats::dnorm(vdom, mean = mean.norm, sd = sd.norm), type = "l",
             xlab = paste0("v", dimension[j]),
             ylab = paste0("p(v", dimension[j], "|D)"),lwd = 2)
        mean.norm <- format(round(object$vmap[dimension[j]], 2), nsmall = 2)
        var.norm <- format(round(diag(object$Cov.vmap)[dimension[j]], 2),
                           nsmall = 2)
        graphics::legend("topleft", paste0("N(", mean.norm, ",", var.norm, ")"), lty = 1,
               lwd = 2, bty = "n", cex = 0.8)
      }
    }

  }

  else stop("Specify a valid object class")

}
