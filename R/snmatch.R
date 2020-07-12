#' Fit a skew-normal distribution to a target density.
#'
#' The routine fits a skew-normal univariate distribution to a target density.
#' Parameters of the resulting skew-normal fit are estimated by the method
#' of moments.
#'
#'
#' @param x A numeric vector on the domain of the target density.
#' @param y The y-coordinates of the target density on grid \code{x}.
#' @param p Vector of probabilities at which to compute quantiles of the
#' skew-normal fit.
#'
#' @details The skew-normal density is parameterized by a location parameter
#' \eqn{\mu}, a scale parameter \eqn{\omega} and a shape parameter
#' \eqn{\rho} that regulates skewness. The probability density function at any
#' x on the real line is:
#'
#' \deqn{p(x) = (2/\omega) \phi((x-\mu)/\omega) \psi(\rho (x-\mu)/\omega),}
#'
#' where \eqn{\phi()} and \eqn{\psi()} denote the standard Gaussian density and
#' cumulative distribution function respectively (see Azzalini 2018).
#' The first moment and second and third central moments of the target density
#' are computed based on the \code{x, y} coordinates using the trapezoidal rule
#' and matched against the theoretical moments of a skew-normal distribution.
#' The solution to this system of equations is the method of moment estimate of
#' the location, scale and shape parameters of a skew-normal density.
#'
#' @return A list with the following components:
#'
#' \item{location \verb{ }}{Estimated location parameter.}
#'
#' \item{scale}{Estimated scale parameter.}
#'
#' \item{shape}{Estimated shape parameter.}
#'
#' \item{snfit}{Fitted values of the skew-normal density
#'   computed on an equally spaced grid between \code{min(x)}
#'   and \code{max(x)}.}
#'
#' \item{quant}{Vector of quantiles of the skew-normal fit
#'   computed on the input vector of probabilities \code{p}.}
#'
#' \item{xgrid}{Equidistant grid on which the skew-normal fitted
#'   density is computed.}
#'
#' @examples
#' # Pdf of skew-normal density
#' sn.target <- function(x, location, scale, shape){
#'                val <- 2 * stats::dnorm(x, mean = location, sd = scale) *
#'                pnorm(shape * (x - location) / scale)
#'                return(val)
#'               }
#'
#' # Extract x and y coordinates from target
#' x.grid <- seq(-2, 6, length = 200)
#' y.grid <- sapply(x.grid, sn.target, location = 0, scale = 2, shape = 3)
#'
#' # Computation of the fit and graphical illustration
#' fit <- snmatch(x.grid, y.grid)
#' domx <- seq(-2, 6, length = 1000)
#' plot(domx, sapply(domx, sn.target, location = 0, scale = 2, shape = 3),
#'      type = "l", ylab = "f(x)", xlab = "x", lwd= 2)
#' lines(fit$xgrid, fit$snfit, type="l", col = "red", lwd = 2, lty = 2)
#' legend("topright", lty = c(1,2), col = c("black", "red"), lwd = c(2, 2),
#'        c("Target","SN fit"), bty="n")
#'
#' # Extract estimated parameters
#' fit$location # Estimated location parameter
#' fit$scale    # Estimated scale parameter
#' fit$shape    # Estimated shape parameter
#'
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' @references Azzalini, A. (2018). The Skew-Normal and Related families.
#' \emph{Cambridge University Press}.
#' @export

snmatch <- function(x, y, p = c(0.025, 0.5, 0.975)){

  if(!is.vector(x, mode = "numeric") || !is.vector(y, mode = "numeric") ||
     !is.vector(p, mode = "numeric"))
    stop("x, y, p must be numeric vectors")
  if(length(x) != length(y))
    stop("x and y must have the same length")
  if(any(is.infinite(x)) | any(is.infinite(y)))
    stop("x, y must contain finite values")
  if(anyNA(x) | anyNA(y) | anyNA(p))
    stop("x or y or p contain NA or NaN")
  if(any(p < 0)  | any(p > 1))
    stop("probabilities must be between 0 and 1")

  # sort of x and y
  x <- sort(x)
  y <- y[order(x)]

  # Cubic root function
  cubroot <- function(x) sign(x) * (abs(x)) ^ (1 / 3)

  # Empirical moments computed using trapezoidal rule
  n <- length(x) - 1
  dx <- diff(x)
  gx <- x * y
  trapezoid.vec <- rowSums(cbind(gx[1:n], gx[-1])) * .5
  mom1 <- sum(trapezoid.vec * dx)
  gx2 <- ((x - mom1) ^ 2) * y
  trapezoid.vec <- rowSums(cbind(gx2[1:n], gx2[-1])) * .5
  mom2 <- sum(trapezoid.vec * dx)
  gx3 <- ((x - mom1) ^ 3) * y
  trapezoid.vec <- rowSums(cbind(gx3[1:n], gx3[-1])) * .5
  mom3 <- sum(trapezoid.vec * dx)

  # Matching skew-normal parameters
  const.kappa <- (cubroot(mom3) * sqrt(pi)) / (cubroot(4 - pi) * (2 ^ (1 / 6)) *
                                                 sqrt(mom2))
  psi.star <- sign(mom3) * sqrt(4 * (const.kappa ^ 2 +
              (2 * (const.kappa ^ 4) / pi))) / (2 + (4 * const.kappa ^ 2) / pi)
  if (psi.star < -1) {
    psi.star <- (-1 + 1e-3)
  }
  if (psi.star > 1) {
    psi.star <- (1 - 1e-3)
  }
  rho.star <- psi.star / (sqrt(1 - psi.star ^ 2))           # shape
  zeta.star <- sqrt(mom2 / (1 - (2 / pi) * (psi.star ^ 2))) # scale
  mu.star <- mom1 - zeta.star * sqrt(2 / pi) * psi.star     # location
  xfine <- seq(min(x), max(x), length = 500)
  snfit <- 2 * stats::dnorm(xfine, mean = mu.star, sd = zeta.star) *
    stats::pnorm(rho.star * (xfine - mu.star) / zeta.star)
  quantiles <- sapply(p, sn::qsn, xi = mu.star, omega = zeta.star,
                      alpha = rho.star)
  list.out <- list(location = mu.star,
                   scale = zeta.star,
                   shape = rho.star,
                   snfit = snfit,
                   quant = quantiles,
                   xgrid = xfine)
  return(list.out)
}

