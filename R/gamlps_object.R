#' Object resulting from the fit of a generalized additive model.
#'
#' An object returned by the \code{\link{gamlps}} function consists in a list
#' with various components related to the fit of a generalized additive model
#' with the Laplace-P-spline approach.
#'
#' @return A \code{gamlps} object has the following elements:
#'
#' \item{formula}{The formula of the generalized additive model.}
#' \item{family}{The chosen exponential family.}
#' \item{link}{The link function used for the fit.}
#' \item{n}{Sample size.}
#' \item{q}{Total number of smooth terms.}
#' \item{K}{Number of B-spline basis functions used for the fit.}
#' \item{penalty.order}{Chosen penalty order.}
#' \item{latfield.dim}{The dimension of the latent field. This is equal
#'  to the sum of the number of B-spline coefficients and the number of
#'  regression parameters related to the covariates in the linear part.}
#' \item{linear.coeff}{Estimated linear regression coefficients. This is a
#'  matrix containing the posterior point estimate, standard deviation, z-score
#'  and lower/upper bounds of the credible interval.}
#' \item{spline.estim}{The estimated B-spline coefficients. This is a list
#'  with \code{q} vectors of size \code{K-1} representing the estimated B-spline
#'  amplitudes for each smooth term.}
#' \item{edf}{Estimated effective degrees of freedom for each latent field
#'  variable.}
#' \item{Approx.signif}{A matrix returning the observed test statistic and
#'  p-value for the approximate significance of smooth terms.}
#' \item{EDf}{The estimated effective degrees of freedom of the smooth terms.}
#' \item{EDfHPD.95}{95\% HPD interval for the degrees of freedom of the smooth
#'  terms.}
#' \item{ED}{The estimated degrees of freedom of the GAM model.}
#' \item{vmap}{The maximum a posteriori of the (log) posterior penalty vector.}
#' \item{Cov.vmap}{Covariance matrix of the (log) posterior penalty vector
#'  evaluated at vmap.}
#' \item{pen.family}{The family of the posterior distribution for v. It is
#'  either "skew-normal" or "gaussian".}
#' \item{pendist.params}{The parameterization for the posterior distribution of
#'  v. If the posterior of v belongs to the skew-normal family, then
#'  \code{pendist.params} is a matrix with as many rows as the number of smooth
#'  terms \code{q}. Each row contains the location, scale and shape parameter
#'  of the skew-normal distribution. If the posterior of v belongs to the
#'  Gaussian family, then \code{pendist.params} is a vector of length \code{q},
#'  corresponding to \code{vmap}.}
#' \item{Covmaximum}{The covariance matrix of the latent field evaluated at the
#'  posterior maximum value of the penalty vector.}
#' \item{latmaximum}{The latent field value evaluated at the posterior
#'  maximum value of the penalty vector.}
#' \item{fitted.values}{The fitted response values.}
#' \item{residuals}{The response residuals.}
#' \item{r2.adj}{The adjusted r-squared of the model indicating the proportion
#'  of the data variance explained by the model fit.}
#' \item{data}{The data frame of the model.}
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @seealso \code{\link{gamlps}}, \code{\link{print.gamlps}},
#' \code{\link{plot.gamlps}}
#'
#' @name gamlps.object
NULL
