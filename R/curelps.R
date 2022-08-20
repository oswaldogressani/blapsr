#' Promotion time cure model with Laplace P-splines.
#'
#' Fits a promotion time cure model with the Laplace-P-spline methodology. The
#' routine can be applied to survival data for which a plateau is observed in
#' the Kaplan-Meier curve. In this case, the follow-up period is considered to
#' be sufficiently long to intrinsically account for long-term survivors and
#' hence a cured fraction. The user can separately specify the model covariates
#' influencing the cure probability (long-term survival) and the population
#' hazard dynamics (short-term survival).
#'
#' @param formula A formula object where the ~ operator separates the response
#'   from the covariates. In a promotion time cure model, it takes the form
#'   \emph{response ~ covariates}, where \emph{response} is a survival object
#'   returned by the \code{\link{Surv}} function of the \emph{survival} package.
#'   The model covariates influencing the long-term survival can be specified
#'   in the function \code{lt(.)} separated by '+', while the
#'   covariates affecting the short-term survival can be specified
#'   in \code{st(.)}. For instance, a promotion time cure model with
#'   covariates specified as \code{lt(x1+x2)+st(x1)}, means that \code{x1} will
#'   jointly influence the long- and short-term survival, while \code{x2} will
#'   only influence the long-term survival.
#' @param data Optional. A data frame to match the variable names provided in
#'  \code{formula}.
#' @param K A positive integer specifying the number of cubic B-spline
#'   functions in the basis. Default is \code{K = 30} and allowed values
#'   are \code{10 <= K <= 60}.
#' @param penorder The penalty order used on finite differences of the
#'  coefficients of contiguous B-splines. Can be either 2 for a second-order
#'  penalty (the default) or 3 for a third-order penalty.
#' @param tmax A user-specified value for the upper bound of the B-spline
#'  basis. The default is NULL, so that the B-spline basis is specified
#'  in the interval \emph{[0, tup]}, where \emph{tup} is the upper bound of
#'  the follow-up times. It is required that \emph{tmax} > \emph{tup}.
#' @param constr Constraint imposed on last B-spline coefficient
#'  (default is 6).
#'
#' @details The log-baseline hazard is modeled as a linear combination of
#'   \code{K} cubic B-splines as obtained from \code{\link{cubicbs}}. A
#'   robust Gamma prior is imposed on the roughness penalty parameter.
#'   A grid-based approach is used to explore the posterior penalty space and
#'   the resulting quadrature points serve to compute the approximate (joint)
#'   posterior of the latent field vector. Point and set estimates of latent
#'   field elements are obtained from a finite mixture of Gaussian densities.
#'   The routine centers the columns of the covariate matrix around their mean
#'   value for numerical stability. See \code{\link{print.curelps}} for a
#'   detailed explanation on the output printed by the curelps
#'   function.
#'
#' @return An object of class \code{curelps} containing various components
#'  from the promotion time cure model fit. Details can be found in
#'  \code{\link{curelps.object}}. Estimates on the baseline survival,
#'  population survival (for a chosen covariate profile) and cure probability
#'  can be obtained with the \code{\link{plot.curelps}} and
#'  \code{\link{curelps.extract}} routines.
#'
#' @examples
#'
#' ## Fit a promotion time cure model on malignant melanoma data
#'
#' data(melanoma)
#' medthick <- median(melanoma$thickness)
#'
#' # Kaplan-Meier estimate to check the existence of a plateau
#' KapMeier <- survfit(Surv(time,status) ~ 1, data = melanoma)
#' plot(KapMeier, mark.time = TRUE, mark = 4, xlab = "Time (in years)")
#'
#' # Fit with curelps
#' fit <- curelps(Surv(time , status) ~ lt(thickness + ulcer) +
#'                    st(thickness + ulcer), data = melanoma, K = 40)
#' fit
#'
#' # Cure prediction for median thickness and absence of ulceration
#' curelps.extract(fit, time = c(2, 4 ,6, 8), curvetype = "probacure",
#'                 cred.int = 0.90, covar.profile = c(medthick, 0, medthick, 0))
#'
#' # Plot of baseline and population survival functions
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 2))
#' # Baseline survival
#' plot(fit, curvetype = "baseline", plot.cred = FALSE, ylim = c(0,1))
#' # Population survival
#' plot(fit, curvetype = "population", covar.profile = c(medthick, 0, medthick, 0),
#' plot.cred = FALSE, ylim = c(0,1))
#' par(opar)
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @seealso \code{\link{curelps.object}}, \code{\link{curelps.extract}},
#'  \code{\link{plot.curelps}}, \code{\link{print.curelps}},
#'  \code{\link{Surv}}.
#'
#' @references Cox, D.R. (1972). Regression models and life-tables.
#'   \emph{Journal of the Royal Statistical Society: Series B (Methodological)}
#'   \strong{34}(2): 187-202. \url{https://doi.org/10.1111/j.2517-6161.1972.tb00899.x}
#' @references Bremhorst, V. and Lambert, P. (2016). Flexible estimation in
#'   cure survival models using Bayesian P-splines.
#'   \emph{Computational Statistical & Data Analysis} \strong{93}: 270-284.
#'   \url{https://doi.org/10.1016/j.csda.2014.05.009}
#' @references Gressani, O. and Lambert, P. (2018). Fast Bayesian inference
#'   using Laplace approximations in a flexible promotion time cure model based
#'   on P-splines. \emph{Computational Statistical & Data Analysis} \strong{124}:
#'   151-167. \url{https://doi.org/10.1016/j.csda.2018.02.007}
#' @references Lambert, P. and Bremhorst, V. (2019). Estimation and
#'   identification issues in the promotion time cure model when the same
#'   covariates influence long- and short-term survival. \emph{Biometrical
#'   Journal} \strong{61}(2): 275-289. \url{https://doi.org/10.1002/bimj.201700250}
#' @export

curelps <- function(formula, data, K = 30, penorder = 2, tmax = NULL,
                    constr = 6){

  if (!inherits(formula, "formula"))
    stop("Incorrect model formula")
  formula.inp <- formula
  if (as.character(stats::terms(formula)[[3]][[2]])[1] == "st" &&
      as.character(stats::terms(formula)[[3]][[3]])[1] == "lt") {
    formula <- stats::as.formula(paste(
      paste0("Surv(", all.vars(formula)[1], ",", all.vars(formula)[2], ")~"),
      attr(stats::terms(formula), "term.labels")[2],
      "+",
      attr(stats::terms(formula), "term.labels")[1], sep = ""))
  }
  varnames <- all.vars(formula, unique = FALSE)
  collect <- parse(text = paste0("list(", paste(varnames, collapse = ","), ")"),
                   keep.source = FALSE)
  if (missing(data)) {
    mff <- data.frame(matrix(unlist(eval(collect)),
                             ncol = length(varnames), byrow = FALSE))
    colnames(mff) <- varnames
    var.type <- unlist(lapply(eval(collect),class))
    which.col <- which(var.type == "factor")
    if (length(which.col) > 0) {
      for (j in 1:length(which.col)) {
        col.idx <- which.col[[j]]
        mff[, col.idx] <- as.factor(mff[, col.idx])
        levels(mff[, col.idx]) <- levels(eval(collect)[[col.idx]])
      }
    }
    list.mff <- list()
    list.mff[[1]] <- mff[, 1:2]
    lt.vars <- as.character(stats::terms(formula)[[3]][[2]])[2] # Variables in "lt"
    st.vars <- as.character(stats::terms(formula)[[3]][[3]])[2] # Variables in "st"
    ncovar.lt <- 0
    ncovar.st <- 0

    for (j in 3:ncol(mff)) {
      col.pos <- match(1, as.numeric(colnames(mff[j]) == colnames(mff)))

      if (is.factor(class(mff[, col.pos]))) {
        if (grepl(colnames(mff)[col.pos], lt.vars, fixed = TRUE)) {
          ncovar.lt <- ncovar.lt + (nlevels(mff[, col.pos]) - 1)
          lt.vars <- gsub(colnames(mff)[col.pos], replacement = "", lt.vars)
        }
        if (grepl(colnames(mff)[col.pos], st.vars, fixed = TRUE)) {
          ncovar.st <- ncovar.st +  (nlevels(mff[, col.pos]) - 1)
          st.vars <- gsub(colnames(mff)[col.pos], replacement = "", st.vars)
        }
        list.inpj <- as.data.frame(stats::model.matrix( ~ 0 + mff[, j])[,-1])
        colnames(list.inpj) <- levels(mff[, j])[-1]
        list.mff[[j - 1]] <- list.inpj
      } else {
        list.inpj <- as.data.frame(mff[, j])
        colnames(list.inpj) <- colnames(mff[j])
        list.mff[[j - 1]] <- list.inpj
        if (grepl(colnames(mff)[col.pos], lt.vars, fixed = TRUE)) {
          ncovar.lt <- ncovar.lt + 1
          lt.vars <- gsub(colnames(mff)[col.pos], replacement = "", lt.vars)
        }
        if (grepl(colnames(mff)[col.pos], st.vars, fixed = TRUE)) {
          ncovar.st <- ncovar.st + 1
          st.vars <- gsub(colnames(mff)[col.pos], replacement = "", st.vars)
        }
      }

    }
    mff <- do.call(cbind.data.frame, list.mff)
  } else {
    mff <- data.frame(matrix(unlist(eval(collect, envir = data)),
                             ncol = length(varnames), byrow = FALSE))
    colnames(mff) <- varnames
    list.mff <- list()
    list.mff[[1]] <- mff[, 1:2]
    lt.vars <- as.character(stats::terms(formula)[[3]][[2]])[2] # Variables in "lt"
    st.vars <- as.character(stats::terms(formula)[[3]][[3]])[2] # Variables in "st"
    ncovar.lt <- 0
    ncovar.st <- 0

    for (j in 3:ncol(mff)) {
      col.pos <- match(1, as.numeric(colnames(mff[j]) == colnames(data)))

      if (is.factor(class(data[, col.pos]))) {
        if (grepl(colnames(data)[col.pos], lt.vars, fixed = TRUE)) {
          ncovar.lt <- ncovar.lt + (nlevels(data[, col.pos]) - 1)
          lt.vars <-
            gsub(colnames(data)[col.pos], replacement = "", lt.vars)
        }
        if (grepl(colnames(data)[col.pos], st.vars, fixed = TRUE)) {
          ncovar.st <- ncovar.st +  (nlevels(data[, col.pos]) - 1)
          st.vars <-
            gsub(colnames(data)[col.pos], replacement = "", st.vars)
        }
        mff[, j] <- as.factor(mff[, j])
        levels(mff[, j]) <- levels(data[, col.pos])
        list.inpj <- as.data.frame(stats::model.matrix( ~ 0 + mff[, j])[,-1])
        colnames(list.inpj) <- levels(mff[, j])[-1]
        list.mff[[j - 1]] <- list.inpj
      } else {
        list.inpj <- as.data.frame(mff[, j])
        colnames(list.inpj) <- colnames(mff[j])
        list.mff[[j - 1]] <- list.inpj
        if (grepl(colnames(data)[col.pos], lt.vars, fixed = TRUE)) {
          ncovar.lt <- ncovar.lt + 1
          lt.vars <- gsub(colnames(data)[col.pos], replacement = "", lt.vars)
        }
        if (grepl(colnames(data)[col.pos], st.vars, fixed = TRUE)) {
          ncovar.st <- ncovar.st + 1
          st.vars <- gsub(colnames(data)[col.pos], replacement = "", st.vars)
        }
      }

    }
    mff <- do.call(cbind.data.frame, list.mff)
  }

  time <- mff[, 1]    # Event times
  event <- mff[, 2]   # Event indicator
  n <- nrow(mff)
  X <- as.matrix(cbind(rep(1, n), mff[, (3:(2 + ncovar.lt))]))
  colnames(X) <- c("(Intercept)", colnames(mff[(3:(2 + ncovar.lt))]))
  if(any(is.infinite(X)))
    stop("Data contains infinite covariates")
  if(ncol(X) == 0)
    stop("The model has no covariates for the long-term survival part")
  Xdata <- X # Uncentered covariate matrix of long-term survival
  X <- scale(X, center = TRUE, scale = FALSE) # Center X matrix
  X[, 1] <- rep(1, n) # Fill column of ones for intercept
  Z <- as.matrix(mff[, (2 + ncovar.lt + 1) : ncol(mff)], ncol = ncovar.st)
  if(any(is.infinite(Z)))
    stop("Data contains infinite covariates")
  if(ncol(Z) == 0)
    stop("The model has no covariates for the short-term survival part")
  colnames(Z) <- colnames(mff[(2 + ncovar.lt + 1) : ncol(mff)])
  Zdata <- Z # Uncentered covariate matrix of short-term survival
  Z <- scale(Z, center = TRUE, scale = FALSE) # Center Z matrix
  if (!is.vector(K, mode = "numeric") || length(K) > 1 || is.na(K))
    stop("K must be a numeric scalar")
  if (K < 10 || K > 60)
    stop("K must be between 10 and 60")
  tlow <- 0                  # Lower bound of follow-up
  tup <- max(time)           # Upper bound of follow-up
  if(!is.null(tmax)){
    if(tmax <= 0) stop("tmax must be positive")
    else if(tmax < tup) stop("tmax must be larger than the largest observed
                             event time")
    tup <- tmax
  }
  nbeta <- ncol(X)           # Number of regression coeff. in proba to be cured
  ngamma <- ncol(Z)          # Number of regression coeff. in Cox PH part
  n <- length(time)          # Sample size
  H <- K + nbeta + ngamma    # Latent field dimension
  if(n < H)
    warning("Number of coefficients to be estimated is larger than sample size")
  constr <- constr                # Constraint on last B-spline coefficient
  penorder <- floor(penorder)
  if(penorder < 2 || penorder > 3)
    stop("Penalty order must be either 2 or 3")


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
    Sigmastar <- Laplace$Sigmastar   # Unconstrained latent  covariance matrix
    Sigstar.c <- Laplace$Sigstar.c   # Constrained covariance matrix
    latstar <- Laplace$latstar       # Unconstrained mode of Laplace approx.
    latstar.c <- Laplace$latstar.c   # Constrained mode of dimension H-1
    latstar.cc <- Laplace$latstar.cc # Constrained mode of dimension H

    # Computation of the log-posterior of v
    Q.v <- Qv(v)
    LQ <- t(chol(Q.v))
    LSig <- t(chol(Sigstar.c))

    logpv <- loglik(latstar.cc) - .5 *  sum((latstar.cc * Q.v) %*% latstar.cc) +
      sum(log(diag(LQ))) + sum(log(diag(LSig))) + 5 * nu * v -
      (.5 * nu + a.delta) * log(b.delta + .5 * nu * exp(v))
    varbetas <- diag(Sigmastar)[(K + 1):(K + nbeta)]
    vargammas <- diag(Sigmastar)[(K + nbeta + 1):H]

    list.out <- list(y = logpv, # Log- posterior of v
                     latstar = latstar, # Posterior mode of Laplace approx.
                     latstar.c = latstar.c, # Constrained posterior mode
                     latstar.cc = latstar.cc, # Constrained posterior mode
                     Sigstar.c = Sigstar.c, # Constrained covariance matrix
                     varbetas = varbetas, # Variance of beta coeffs.
                     vargammas = vargammas) # Variance of gamma coeffs.
    return(list.out)
  }

  #Grid construction
  grid.construct <- function() {
    v.grid <- c()
    logpv.v <- c()
    iter <- 1
    jump <- 3
    v.0 <- 17
    penlow <- 0
    latent0 <-
      rep(sum(event) / sum(time), H) # Initial latent field vector
    logpv <- logpost.v(v.0, latent0)
    logpv.v[iter] <- logpv$y
    latent0 <- logpv$latstar
    v.grid[iter] <- v.0

    # First nesting of vmap
    while (jump > 0) {
      iter <- iter + 1
      v.0 <- v.0 - 1
      v.grid[iter] <- v.0
      logpv <- logpost.v(v.0, latent0)
      logpv.v[iter] <- logpv$y
      latent0 <- logpv$latstar
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
      omega.m <-
        exp(logpvquad - logpvmap) / sum(exp(logpvquad - logpvmap))
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
      latent0 <- logpv$latstar

      while (jump > 0) {
        iter <- iter + 1
        v.0 <- v.0 - 0.1
        v.grid[iter] <- v.0
        logpv <- logpost.v(v.0, latent0)
        logpv.v[iter] <- logpv$y
        latent0 <- logpv$latstar
        jump <- logpv.v[iter] -  logpv.v[iter - 1]
      }

      v.grid <- rev(v.grid)
      logpv.v <- rev(logpv.v)

      if (length(v.grid) >= 3) {
        vmap <- (v.grid[1] + v.grid[3]) / 2
      } else{
        vmap <- (v.grid[1] + v.grid[2]) / 2
      }
      logpvmap <- logpost.v(vmap, latent0)$y
      nquad <- 8
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
      omega.m <-
        exp(logpvquad - logpvmap) / sum(exp(logpvquad - logpvmap))
    }

    listout <- list(
      vmap = vmap,
      logpvmap = logpvmap,
      nquad = nquad,
      vquadrature = vquadrature,
      logpvquad = logpvquad,
      omega.m = omega.m,
      list.star = list.star,
      latent.0 = latent0
    )
    return(invisible(listout))
  }
  gridpoints <- grid.construct()
  vmap <- gridpoints$vmap
  vquadrature <- gridpoints$vquadrature
  nquad <- gridpoints$nquad
  logpvquad <- gridpoints$logpvquad
  omega.m <- gridpoints$omega.m
  latent0 <- gridpoints$latent.0
  logpvmax <- logpost.v(vmap, latent0)$y
  list.star <- gridpoints$list.star

  # Posterior inference
  latent.matrix <- matrix(unlist(lapply(list.star, "[[", 4)),
                          ncol = H, byrow = TRUE)
  varbetas <- matrix(unlist(lapply(list.star, "[[", 6)),
                     ncol = nbeta, byrow = TRUE)
  vargammas <- matrix(unlist(lapply(list.star, "[[", 7)),
                      ncol = ngamma, byrow = TRUE)
  Covlatc.map <- logpost.v(vmap, gridpoints$latent.0)$Sigstar.c
  Covtheta.map <- Covlatc.map[1:(K - 1), 1:(K - 1)]

  # Posterior point estimate of latent field vector
  latent.hat <- colSums(omega.m * latent.matrix)
  theta.hat <- latent.hat[1:K]
  beta.hat <- latent.hat[(K + 1):(K + nbeta)]
  gamma.hat <- latent.hat[(K + nbeta + 1):H]

  # Posterior standard deviation of regression coefficients
  sdbeta.post <- sqrt(colSums(omega.m * varbetas) +
                        colSums(omega.m * (
                          latent.matrix[, (K + 1):(K + nbeta)] -
                            matrix(rep(beta.hat, nquad), nrow = nquad,
                                   byrow = TRUE)) ^ 2))

  sdgamma.post <- sqrt(colSums(omega.m * vargammas) +
                         colSums(omega.m * (
                           latent.matrix[, (K + nbeta + 1):H] -
                             matrix(rep(gamma.hat, nquad), nrow = nquad,
                                    byrow = TRUE)) ^ 2))

  # Wald test statistic
  z.Waldcure <- beta.hat / sdbeta.post # Wald statistic for cure probability
  z.Waldcox <- gamma.hat / sdgamma.post # Wald statistic for Cox model

  # Posterior credible intervals for beta coefficients
  postbetas.CI95 <- matrix(0, nrow = nbeta, ncol = 2)

  for (j in 1:nbeta) {
    means.betaj <- latent.matrix[, (K + j)]
    variances.betaj <- varbetas[, j]
    betaj.dom <- seq(beta.hat[j] - 4 * sdbeta.post[j],
                     beta.hat[j] + 4 * sdbeta.post[j],
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
    postbetas.CI95[j,] <- c(lb95.j, ub95.j)
  }

  # Posterior credible intervals for gamma coefficients
  postgammas.CI95 <- matrix(0, nrow = ngamma, ncol = 2)

  for (j in 1:ngamma) {
    means.gammaj <- latent.matrix[, (K + nbeta + j)]
    variances.gammaj <- vargammas[, j]
    gammaj.dom <- seq(gamma.hat[j] - 4 * sdgamma.post[j],
                      gamma.hat[j] + 4 * sdgamma.post[j],
                      length = 500)
    dj <- gammaj.dom[2] - gammaj.dom[1]

    post.gammaj <- function(gammaj) {
      mixture <- sum(omega.m * stats::dnorm(gammaj, mean = means.gammaj,
                                     sd = sqrt(variances.gammaj)))
      return(mixture)
    }
    mix.image <- sapply(gammaj.dom, post.gammaj)
    cumsum.image <- cumsum(mix.image * dj)
    lb95.j <- gammaj.dom[sum(cumsum.image < 0.025)]
    ub95.j <- gammaj.dom[sum(cumsum.image < 0.975)]
    postgammas.CI95[j,] <- c(lb95.j, ub95.j)
  }

  # Output for beta coefficients (related to proba to be cured or long-term survival)
  coeff.probacure <- cbind(beta.hat, sdbeta.post, z.Waldcure,
                           postbetas.CI95[, 1], postbetas.CI95[, 2])
  colnames(coeff.probacure) <- c("coef", "sd.post", "z",
                                 "lower.95", "upper.95")
  rownames(coeff.probacure) <- colnames(X)

  # Output for gamma coefficients (related to cox model or short-term survival)
  coeff.cox <- cbind(gamma.hat, exp(gamma.hat), sdgamma.post, z.Waldcox,
                     exp(gamma.hat), exp(-gamma.hat), exp(postgammas.CI95[, 1]),
                     exp(postgammas.CI95[, 2]))
  colnames(coeff.cox) <- c("coef", "exp(coef)", "sd.post", "z", "exp(coef)",
                           "exp(-coef)", "lower.95", "upper.95")
  rownames(coeff.cox) <- colnames(Z)

  # Effective model dimension
  I.loglik <- (-1) * Hessian.loglik(latent.hat)
  EDmat <- solve(I.loglik + Qv(vmap)) %*% I.loglik
  edf <- diag(EDmat)
  names(edf) <- c(paste0(rep("spline.coeff"), seq_len(K)), colnames(X),
                  colnames(Z))
  ED <- sum(edf)

  # AIC and BIC
  p <- nbeta + ngamma
  loglik.value <- loglik(latent.hat)
  AIC.p <-  (-2 * loglik.value) + 2 * p
  AIC.ED <- (-2 * loglik.value) + 2 * ED
  BIC.p <-  (-2 * loglik.value) + p * log(sum(event))
  BIC.ED <- (-2 * loglik.value) + ED * log(sum(event))

  # Output results
  listout <- list(formula = formula.inp,       # Model formula
                  K = K,                       # Number of B-splines in basis
                  penalty.order = penorder,    # Penalty order
                  latfield.dim = H,            # Latent field dimension
                  event.times = time,          # Event times
                  n = n,                       # Sample size
                  num.events = sum(event),     # Number of events
                  tup = tup,                   # Upper bound of follow-up
                  event.indicators = event,    # Event indicators
                  coeff.probacure = coeff.probacure, # Coeffs. of proba cure
                  coeff.cox = coeff.cox,       # Cox model part
                  vmap = vmap,                 # Mamixmum a posteriori for v
                  vquad = vquadrature,         # Quadrature points for mixture
                  spline.estim = theta.hat,    # Estimated B-spline amplitudes
                  edf = round(edf, 4),         # Degrees of freedom of latent elements
                  ED = ED,                     # Effective model dimension
                  Covtheta.map = Covtheta.map, # Covariance of spline coeffs
                  Covlatc.map = Covlatc.map,   # Constrained latent covariance
                  X = Xdata,       # Covariate matrix of long-term survival
                  Z = Zdata,       # Covariate matrix of short-term survival
                  loglik = loglik.value, # Log-likelihood at estimated parameters
                  p = p,           # Number of parametric terms
                  AIC.p = AIC.p,   # AIC with penalty term p
                  AIC.ED = AIC.ED, # AIC with penalty term ED
                  BIC.p = BIC.p,   # BIC with penalty term p
                  BIC.ED = BIC.ED) # BIC with penalty term ED
  attr(listout, "class") <- "curelps"
  listout
}

























