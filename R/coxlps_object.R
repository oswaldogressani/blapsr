#' Object from a Cox proportional hazards fit with Laplace-P-splines.
#'
#' An object returned by the \code{\link{coxlps}} function consists in a list
#' with various components related to the fit of a Cox model using the
#' Laplace-P-spline methodology.
#'
#' @return A \code{coxlps} object has the following elements:
#'
#' \item{formula}{The formula of the Cox model.}
#' \item{K}{Number of B-spline basis functions used for the fit.}
#' \item{penalty.order}{Chosen penalty order.}
#' \item{latfield.dim}{The dimension of the latent field. This is equal
#'  to the sum of the number of B-spline coefficients and the number of
#'  regression parameters related to the covariates.}
#' \item{n}{Sample size.}
#' \item{num.events}{The number of events that occurred.}
#' \item{event.times}{The standardized event times, i.e. if \emph{t} denotes
#'  the original time scale, then \code{event.times = t / sd(t)}, where
#'  \code{sd} is the standard deviation.}
#' \item{tup}{The upper bound of the follow-up, i.e. \code{max(event.times)}.}
#' \item{sd.time}{The standard deviation of the event times in original scale.}
#' \item{event.indicators}{The event indicators.}
#' \item{regcoeff}{Posterior estimates of the regression coefficients.
#' \emph{coef} gives the point estimate, \emph{sd.post} gives the posterior
#'  standard deviation, \emph{z} is the Wald test statistic, \emph{lower .95}
#'  and \emph{ upper .95} the posterior approximate 95\% quantile-based credible
#'  interval.}
#' \item{penalty.vector}{The selected grid of penalty values.}
#' \item{vmap}{The maximum a posteriori of the (log) penalty parameter.}
#' \item{spline.estim}{The estimated B-spline coefficients.}
#' \item{edf}{Estimated effective degrees of freedom for each latent field
#'  variable.}
#' \item{ED}{The effective model dimension.}
#' \item{Covthetamix}{The posterior covariance matrix of the B-spline
#'  coefficients.}
#' \item{X}{The matrix of covariate values.}
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
#' @seealso \code{\link{coxlps}}, \code{\link{coxlps.baseline}}
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @name coxlps.object
NULL
