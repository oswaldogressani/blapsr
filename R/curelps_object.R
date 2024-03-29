#' Object from a promotion time model fit with Laplace-P-splines.
#'
#' An object returned by the \code{\link{curelps}} function consists in a list
#' with various components related to the fit of a promotion time cure model
#' using the Laplace-P-spline methodology.
#'
#' @return A \code{curelps} object has the following elements:
#'
#' \item{formula}{The formula of the promotion time cure model.}
#' \item{K}{Number of B-spline basis functions used for the fit.}
#' \item{penalty.order}{Chosen penalty order.}
#' \item{latfield.dim}{The dimension of the latent field. This is equal
#'  to the sum of the number of B-spline coefficients and the number of
#'  regression parameters related to the covariates.}
#' \item{event.times}{The observed event times.}
#' \item{n}{Sample size.}
#' \item{num.events}{The number of events that occurred.}
#' \item{tup}{The upper bound of the follow up, i.e. \code{max(event.times)}.}
#' \item{event.indicators}{The event indicators.}
#' \item{coeff.probacure}{Posterior estimates of the regression coefficients
#'  related to the cure probability (or long-term survival).}
#' \item{coeff.cox}{Posterior estimates of the regression coefficients
#'  related to the population hazard dynamics (or short-term survival).}
#' \item{vmap}{The maximum a posteriori of the (log-)posterior penalty parameter.}
#' \item{vquad}{The quadrature points of (log-) posterior penalty parameters
#'  used to compute the Gaussian mixture posterior of the latent field vector.}
#' \item{spline.estim}{The estimated B-spline coefficients.}
#' \item{edf}{Estimated effective degrees of freedom for each latent field
#'  variable.}
#' \item{ED}{The effective model dimension.}
#' \item{Covtheta.map}{The posterior covariance matrix of the B-spline
#'  coefficients for a penalty fixed at its maximum posterior value.}
#' \item{Covlatc.map}{The posterior covariance matrix of the latent field
#'  for a penalty fixed at its maximum posterior value.}
#' \item{X}{The covariate matrix for the long-term survival part.}
#' \item{Z}{The covariate matrix for the short-term survival part.}
#' \item{loglik}{The log-likelihood evaluated at the posterior latent field
#'  estimate.}
#' \item{p}{Number of parametric coefficients in the model.}
#' \item{AIC.p}{The AIC computed with the formula \emph{-2*loglik+2*p},
#'  where \emph{p} is the number of parametric coefficients.}
#' \item{AIC.ED}{The AIC computed with the formula \emph{-2*loglik+2*ED}, where
#'  \emph{ED} is the effective model dimension.}
#' \item{BIC.p}{The BIC computed with the formula \emph{-2*loglik+p*log(ne)},
#'  where \emph{p} is the number of parametric coefficients and \emph{ne} the
#'  number of events.}
#' \item{BIC.ED}{The BIC computed with the formula \emph{-2*loglik+ED*log(ne)},
#'  where \emph{ED} is the effective model dimension and \emph{ne} the
#'  number of events.}
#'
#' @seealso \code{\link{curelps}}
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @name curelps.object
NULL
