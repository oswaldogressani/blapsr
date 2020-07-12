#' Construct a cubic B-spline basis.
#'
#' Computation of a cubic B-spline basis matrix.
#'
#' @param x A numeric vector containing the values on which to evaluate the
#'          B-spline basis.
#'
#' @param lower,upper The lower and upper bounds of the B-spline basis domain.
#'                     Must be finite with \code{lower < upper}.
#' @param K A positive integer specifying the number of B-spline functions in
#'          the basis.
#' @author Gressani Oswaldo \email{oswaldo_gressani@hotmail.fr}.
#'
#' The core algorithm of the \code{cubicbs} function owes much to a code written
#' by Phlilippe Lambert.
#
#' @return An object of class \code{cubicbs} for which
#' \code{print} and \code{plot} methods are available. The \code{cubicbs}
#' class consists of a list with the following components:
#'
#' \item{x}{A numeric vector on which the basis is evaluated.}
#'
#' \item{lower, upper \verb{ }}{The lower and upper bounds of the basis domain.}
#'
#' \item{K}{The number of cubic B-spline functions in the basis.}
#'
#' \item{knots}{The knot sequence to build the basis.}
#'
#' \item{nknots}{Total number of knots.}
#'
#' \item{dimbasis}{The dimension of the B-spline basis matrix.}
#'
#' \item{Bmatrix}{The B-spline basis matrix.}
#'
#' The \code{print} method summarizes the B-spline basis and the \code{plot}
#' method gives a graphical representation of the basis
#' with dashed vertical lines indicating knot placement and blue ticks the
#' coordinates of \code{x}.
#'
#' @examples
#' lb <- 0  # Lower bound
#' ub <- 1  # Upper bound
#' xdom <- runif(100, lb, ub) # Draw uniform values between lb and ub
#' Bsmat <- cubicbs(xdom, lb, ub, 25) # 100 x 25 B-spline matrix
#' Bsmat
#' plot(Bsmat) # Plot the basis
#' @references Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with
#' B-splines and penalties. \emph{Statistical Science}, \strong{11}(2): 89-121.
#' @export

cubicbs <- function(x, lower, upper, K) {
  if (!is.vector(x, mode = "numeric"))
    stop("x must be a numeric vector")
  if (anyNA(x))
    stop("x cannot contain NA or NaN values")
  if (!is.vector(lower, mode = "numeric") ||
      !is.vector(upper, mode = "numeric"))
    stop("Lower bound and/or upper bound is not a numeric vector")
  if (length(lower) > 1 || length(upper) > 1)
    stop("Lower bound and/or upper bound must be of length 1")
  if (is.infinite(lower) || is.infinite(upper))
    stop("Lower bound and/or upper bound must be finite")
  if (lower >= upper)
    stop("Lower bound must be smaller than upper bound")
  if (any(x < lower) || any(x > upper))
    stop("values in x must be between lower and upper")
  if (!is.vector(K, mode = "numeric") || length(K) > 1)
    stop("K must be a numeric vector of length 1")
  if (is.na(K))
    stop("K cannot be NA or NaN")
  if (floor(K) <= 3 || is.infinite(K))
    stop("K must be a finite integer larger than 3")

  nx <- length(x) # length of input vector
  B <- matrix(0, nrow = nx, ncol = K) # dimension of B-spline matrix
  dimB <- dim(B) # dimension of B-spline matrix
  ndx <- K - 3 # number of intervals between lower and upper
  dx <- (upper - lower) / ndx # interval width
  nknots <- ndx + 2 * 3 + 1 # total number of knots in basis
  knots <- seq(lower - 3 * dx, upper + 3 * dx, by = dx) # knots

  for (i in 1:nx) {
    for (j in 1:(nknots - 4)) {
      temp <- 0
      cub  <- x[i] - knots[j]
      if (cub > 0) {
        temp <- temp + cub ^ 3
        cub  <- x[i] - knots[j + 1]
        if (cub > 0) {
          temp <- temp - 4 * cub ^ 3
          cub  <- x[i] - knots[j + 2]
          if (cub > 0) {
            temp <- temp + 6 * cub ^ 3
            cub  <- x[i] - knots[j + 3]
            if (cub > 0) {
              temp <- temp - 4 * cub ^ 3
              cub  <- x[i] - knots[j + 4]
              if (cub > 0) {
                temp <- temp + cub ^ 3
              }
            }
          }
        }
      }
      B[i, j] <- temp / (6 * dx ^ 3)
      if (abs(B[i, j]) < 1e-10)
        B[i, j] <- 0
    }
  }

  listout <- list(x = x,
                  lower = lower,
                  upper = upper,
                  K = K,
                  knots = knots,
                  nknots = nknots,
                  dimbasis = dimB,
                  Bmatrix = B)

  attr(listout, "class") <- "cubicbs"
  listout
}


