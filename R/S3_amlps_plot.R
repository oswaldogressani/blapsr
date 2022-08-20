#' Plot smooth functions of an additive model object.
#'
#' @description Displays a plot of the fitted additive smooth components of an
#'  \code{amlps} object. The routine can also be used to
#'  print a table of point and set estimates of an additive smooth term for a
#'  user-specified grid of values.
#'
#' @usage
#' \method{plot}{amlps}(x, xp, smoo.index, cred.int = 0.95, plot.cred = TRUE,
#'      np = 100, fit.col = "blue", shade.col = "gray75", show.plot = TRUE,
#'      show.info = TRUE, ...)
#'
#' @param x An object of class \code{amlps}.
#' @param xp A numeric vector of grid values on which to compute a point estimate
#'  and pointwise credible interval for the smooth function specified in
#'  \code{smoo.index}. The components of \code{xp} must be within the range
#'  of the observed covariate values for the corresponding smooth function.
#'  Results will be displayed in a table.
#' @param smoo.index The index of the smooth function. For instance
#'  \code{smoo.index = 2} refers to the second smooth function specified in
#'  the \code{formula} of the \code{amlps} routine.
#' @param cred.int The level of the pointwise credible interval to be
#'  computed for the smooth additive term. Default is \code{0.95}.
#' @param plot.cred Logical. Should the credible intervals be plotted?
#'  Default is \code{TRUE}.
#' @param np The number of points used to construct the plot of the smooth
#'  additive function. Default is 100 and allowed values are between 20 and
#'  200.
#' @param fit.col The color of the fitted curve.
#' @param shade.col The shading color for the credible intervals.
#' @param show.plot Logical. Should the plot be displayed? Default is
#'  \code{TRUE}.
#' @param show.info Logical. Should the table of point and set estimates of
#'  the smooth function on the specified \code{xp} values be displayed? Default
#'  is \code{TRUE}.
#' @param ... Further arguments to be passed to \code{plot}.
#'
#' @details Produces a plot of a smooth additive term fitted with the
#'  \code{\link{amlps}} function. On the y-axis, the estimated effective
#'  dimension of the smooth term is also displayed. At the bottom of each
#'  plot, vertical ticks indicate the location of the covariate values. The
#'  labels on the x-axis correspond to the covariate name associated to the
#'  smooth term.
#'
#' @return If \code{xp} is unspecified (the default), the routine will only
#' return a plot of the estimated smooth curve. Otherwise, it provides a
#' list with the following components:
#'
#' \item{xp}{The chosen points on which to compute the smooth fit.}
#'
#' \item{sm.xp}{The estimated smooth fit at points specified in \code{xp}.}
#'
#' \item{sm.low}{The lower bound of the pointwise credible interval for the
#'  smooth additive function at points specified in \code{xp}.}
#'
#' \item{sm.up}{The upper bound of the pointwise credible interval for the
#'  smooth additive function at points specified in \code{xp}.}
#'
#'  \item{cred.int}{The chosen level to compute credible intervals.}
#'
#'  \item{smoo.index}{The index of the smooth function.}
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @seealso \code{\link{amlps}}, \code{\link{amlps.object}},
#'  \code{\link{print.amlps}}
#'
#' @examples
#'
#' ### Classic simulated data example
#'
#' set.seed(3)
#'
#' sim.data <- simgamdata(setting = 2, n = 200, dist = "gaussian", scale = 0.3)
#' plot(sim.data)         # Scatter plot of response
#' data <- sim.data$data  # Simulated data frame
#'
#' # Fit model
#' fit <- amlps(y ~ z1 + z2 + sm(x1) + sm(x2), data = data, K = 20)
#' fit
#'
#' # Plot fit of second function and results for a specific grid x
#' plot(fit, xp = c(-0.8, -0.4, 0, 0.4, 0.8), smoo.index = 2, ylim=c(-3, 3))
#'
#' @export

plot.amlps <- function(x, xp, smoo.index, cred.int = 0.95, plot.cred = TRUE,
                       np = 100, fit.col = "blue", shade.col = "gray75",
                       show.plot = TRUE, show.info = TRUE, ...){
  if (smoo.index < 1 || smoo.index > x$q)
    stop("smoo.index wrongly specified")
  smoo.index <- as.integer(smoo.index)
  if (!is.vector(cred.int, mode = "numeric") || length(cred.int) > 1 ||
      is.na(cred.int) || is.infinite(cred.int) || cred.int <= 0 ||
      cred.int >= 1)
    stop("cred.int must be between 0 and 1")
  if (!is.logical(plot.cred))
    stop("plot.cred must be either TRUE or FALSE")
  if(np < 20 || np > 200)
    stop("choose np between 20 and 200")

  ### Pointwise estimation of the smooth additive terms
  if(is.null(x$data)) {
    mf <- stats::model.frame(x$formula) # Extract model frame from formula
    X  <- stats::model.matrix(mf)     # Full design matrix
    colXnames <- colnames(X)
    smterms <- grepl("sm(", colnames(X), fixed = TRUE)
    X <- cbind(X[, as.logical(1 - smterms)], X[, smterms])
    colnames(X) <- colXnames
  } else{
    mf <- stats::model.frame(x$formula, data = x$data)
    X <- stats::model.matrix(mf, data = x$data)
    colXnames <- colnames(X)
    smterms <- grepl("sm(", colnames(X), fixed = TRUE)
    X <- cbind(X[, as.logical(1 - smterms)], X[, smterms])
    colnames(X) <- colXnames
  }
  q  <- x$q      # Number of smooth terms in model
  p  <- ncol(X) - q   # Number of regression coefficients in linear part
  n  <- x$n      # Sample size
  K <- x$K       # Number of B-splines in basis
  splines <- x$spline.estim # Estimated spline coefficients

  # Estimate smooth functions
  j <- smoo.index
  xj <- as.numeric(X[, p + j])
  min.xgrid <- min(xj) # lower bound for B-spline basis
  max.xgrid <- max(xj) # upper bound for B-spline basis
  if(!missing(xp)){
    if (any(xp < min.xgrid) || any(xp > max.xgrid))
      stop("values in xp not in the range of observed covariates")
  }
  xgrid <- seq(min.xgrid, max.xgrid, length = np) # fine grid
  xj.fine <- seq(min.xgrid, max.xgrid, length = 1000) # for centering
  Bj.fine <- cubicbs(xj.fine, lower = min.xgrid, upper = max.xgrid,
                     K = K)$Bmatrix
  Bj.fine.mean <- colMeans(Bj.fine)
  Bxgrid <- cubicbs(xgrid, lower = min.xgrid,
                    upper = max.xgrid, K = K)$Bmatrix
  Bxgrid.centered <- Bxgrid - matrix(rep(Bj.fine.mean, np), nrow = np,
                                     byrow = TRUE)
  Bx <- Bxgrid.centered[, -K]
  fhat  <- as.numeric(Bx %*% splines[[j]])

  ### Pointwise credible intervals for functions fj on xgrid
  alpha <- 1 - cred.int
  latmaximum <- x$latmaximum
  Covmaximum <- x$Covmaximum
  thetaj.max <- latmaximum[((p + 1) + (j - 1) * (K - 1)) : (p + j * (K - 1))]
  Sigj.max <- Covmaximum[((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1)),
                         ((p + 1) + (j - 1) * (K - 1)):(p + j * (K - 1))]
  postfj.mean <- as.numeric(Bx %*% thetaj.max)
  postfj.sd <- sqrt(diag(Bx %*% Sigj.max %*% t(Bx)))
  fj.lb <- fhat - stats::qnorm((1 - (alpha * .5))) * postfj.sd
  fj.ub <- fhat + stats::qnorm((1 - (alpha * .5))) * postfj.sd


  # Plot the desired smooth additive term
  minf <- min(fj.lb)
  maxf <- max(fj.ub)
  shift.frame <- 0.3 * (maxf - minf)
  covariate.name <- colnames(X)[(p + 1):(p + q)][smoo.index]
  nchar.covariate <- nchar(covariate.name)
  covariate.name <- substr(covariate.name, 4, nchar.covariate - 1)

  if(show.plot == TRUE){
    graphics::plot(xgrid, fhat, type = "l", col = fit.col, xlab =  covariate.name,
         ylab = paste0("sm(", covariate.name, ",",
                       format(round(x$EDf[smoo.index], 2), nsmall = 2),
                       ")"), ...)
    if(plot.cred == TRUE) {
      graphics::polygon(x = c(xgrid, rev(xgrid)), y = c(fj.lb, rev(fj.ub)),
              col = shade.col, border = NA)
    }
    graphics::lines(xgrid, fhat, type="l", lwd=2, col = fit.col)
    graphics::rug(X[, p + smoo.index])
  }

  ### Pointwise credible intervals for values in x
  if (!missing(xp)) {
    xxgrid <- xp
    xlen <- length(xp)
    Bxxgrid <- cubicbs(xp, lower = min.xgrid, upper = max.xgrid, K = K)$Bmatrix
    Bxxgrid.centered <- Bxxgrid - matrix(rep(Bj.fine.mean, xlen), nrow = xlen,
                                         byrow = TRUE)
    Bxx <- Bxxgrid.centered[, -K]
    fhatxx  <- as.numeric(Bxx %*% splines[[j]])
    postfjxx.mean <- as.numeric(Bxx %*% thetaj.max)
    postfjxx.sd <- sqrt(diag(Bxx %*% Sigj.max %*% t(Bxx)))
    fjxx.lb <- fhatxx - stats::qnorm((1 - (alpha * .5))) * postfjxx.sd
    fjxx.ub <- fhatxx + stats::qnorm((1 - (alpha * .5))) * postfjxx.sd

    ftable <- matrix(0, nrow = xlen, ncol = 4)
    colnames(ftable) <- c("xp", "sm", "sm.low", "sm.up")
    ftable[, 1] <- xp
    ftable[, 2] <- fhatxx
    ftable[, 3] <- fjxx.lb
    ftable[, 4] <- fjxx.ub
    ftable <- round(ftable, 4)

    if(show.info == TRUE){
      cat("Estimated smooth function", paste0("sm(", covariate.name,")"),
          "at specified grid points (*): \n")
      cat("\n")
      print.table(format(ftable, nsmall = 4), right = TRUE)
      cat("--- \n")
      cat("* Bounds correspond to a", paste(format(round(cred.int  * 100, 2),
                  nsmall = 2), "%", sep = ""), "credible interval. \n")
    }

    listout <- list(xp = xp,
                    sm.xp = fhatxx,
                    sm.low = fjxx.lb,
                    sm.up = fjxx.ub,
                    cred.int = cred.int,
                    smoo.index = smoo.index)

    return(invisible(listout))
  }
}
