#' Test positive definiteness and adjust positive definite matrix.
#'
#' The routine checks if a real symmetric matrix is positive definite and
#' if required adjusts for positive definiteness.
#'
#' @keywords internal
#' @param x A real symmetric matrix.
#' @param eigentol A tolerance level for the eigenvalues below which
#'        the input matrix is not considered to be positive definite.
#' @param correct logical; if TRUE a corrected positive definite matrix
#'        is returned.
#' @details  In statistical models, Newton's method is often used to solve a
#' (multivariate) maximization problem. For this iterative technique to
#' converge, the search direction must be controlled to guarantee an ascent
#' in the objective function at each step. For maximization with Newton-like
#' methods, an iteration yields an ascent provided that the negative of the
#' Hessian matrix is positive definite. When an iteration or initial value
#' choice does note lie in a ``small'' neighborhood of the maximum, positive
#' definiteness of the negative Hessian may not be guaranteed and a correction
#' should be implemented to ensure that the trajectory points in the ascent
#' direction.
#' This routine checks if a symmetric input matrix is positive definite and
#' if required computes a positive definite adjusted matrix following the work
#' of Levenberg (1944), Marquardt (1963) and Goldfeld (1966). The correction
#' consists in perturbing the main diagonal of the original input matrix
#' \code{x} by adding to each diagonal entry the absolute value of the most
#' negative eigenvalue of \code{x} incremented by \emph{10^(-4)}.
#'
#' @return A list with the following components:
#'
#' \code{isPD} \verb{   } logical; if TRUE the input matrix \code{x} is
#'             positive definite.
#'
#' \code{PD}  \verb{     } The corrected positive definite matrix.
#'
#' @examples
#' # Generate a 3 x 3 matrix that fails to be positive definite
#' A <- matrix(c( - 0.478, - 0.013, 0.001, - 0.013, 1.256, 0.001,
#'             0.001, 0.001, 0.024), ncol = 3, byrow = TRUE)
#' adjustPD(A)
#' @references Levenberg, K. (1944). A method for the solution of certain
#' non-linear problems in least squares,
#' \emph{Quarterly of Applied Mathematics} \strong{2}(2): 164-168.
#' @references Marquardt, D. W. (1963). An algorithm for least-squares
#' estimation of nonlinear parameters. \emph{Journal of the society for
#' Industrial and Applied Mathematics} \strong{11}(2): 431-441.
#' @references Goldfeld, S. M., Quandt, R. E., and Trotter, H. F. (1966).
#' Maximization by Quadratic Hill-Climbing
#' \emph{Econometrica} \strong{34}(3): 541-551.
#' @export

adjustPD <- function(x, eigentol = 1e-06, correct = TRUE){
  if(!is.matrix(x))
    stop("Input argument must be a matrix")
  if(!isSymmetric(x))
    stop("Input matrix must be symmetric")
  if(ncol(x)!=nrow(x))
    stop("Input matrix must be square")
  if(eigentol <= 0)
    stop("eigentol must be a positive number")
  if(!is.logical(correct))
    stop("correct must be either T or F")
  eigvals <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  checkPD <- !any(eigvals < eigentol) # Check if matrix is positive definite
  if(correct == TRUE){  # Correction to make matrix positive definite
    if(checkPD == FALSE){
      xcorrect <- x + diag(abs(min(eigvals)) + 1e-04, ncol(x))
    }else{xcorrect <- x}
    list_out <- list(isPD = checkPD, PD = xcorrect)
  }else{list_out <- list(isPD = checkPD)}
  return(list_out)
}
