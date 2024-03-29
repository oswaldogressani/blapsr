% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simcuredata.R
\name{simcuredata}
\alias{simcuredata}
\title{Simulation of survival times for the promotion time cure model.}
\usage{
simcuredata(n, censor = c("Uniform", "Weibull"), cure.setting = 1,
            info = TRUE, KapMeier = FALSE)
}
\arguments{
\item{n}{Sample size.}

\item{censor}{The censoring scheme. Either Uniform (the default) or Weibull.}

\item{cure.setting}{A number indicating the desired cure percentage. If
\code{cure.setting = 1} (default) the cure percentage is around 20\%.
With  \code{cure.setting = 2} the cure percentage is around 30\%.}

\item{info}{Should information regarding the simulation setting be printed
to the console? Default is \code{TRUE}.}

\item{KapMeier}{Logical. Should the Kaplan-Meier curve of the generated
data be plotted? Default is \code{FALSE}.}
}
\value{
A list with the following components:

\item{n}{Sample size.}

\item{survdata}{A data frame containing the simulated data.}

\item{beta.coeff}{The regression coefficients pertaining to long-term
 survival.}

\item{gamma.coeff}{The regression coefficients pertaining to short-term
 survival.}

\item{cure.perc}{The cure percentage.}

\item{censor.perc}{The percentage of censoring.}

\item{censor}{The censoring scheme.}

\item{S0}{The baseline survival function under the chosen Weibull
  parameterization.}
}
\description{
Generates right censored time-to-event data with a plateau in the
Kaplan-Meier estimate.
}
\details{
Latent event times are generated following Bender et al. (2005),
 with a baseline distribution chosen to be a Weibull with mean 8 and variance
 17.47. When \code{cure.setting = 1} the regression coefficients of the
 long-term survival part are chosen to yield a cure percentage around 20\%,
 while \code{cure.setting = 2} yields a cure percentage around 30\%.
 Censoring is either governed by a Uniform distribution on the support
 [20, 25] or by a Weibull distribution with shape parameter 3 and
 scale parameter 25.
}
\examples{
set.seed(10)
sim <- simcuredata(n = 300, censor = "Weibull", KapMeier = TRUE)

}
\references{
Bender, R., Augustin, T. and Blettner, M. (2005). Generating
 survival times to simulate Cox proportional hazards models,
 \emph{Statistics in Medicine} \strong{24}(11): 1713-1723.

Bremhorst, V. and Lambert, P. (2016). Flexible estimation in
  cure survival models using Bayesian P-splines. \emph{Computational
  Statistics & Data Analysis} \strong{93}: 270-284.

Gressani, O. and Lambert, P. (2018). Fast Bayesian inference
  using Laplace approximations in a flexible promotion time cure model based
  on P-splines. \emph{Computational Statistics & Data Analysis} \strong{124}:
  151-167.
}
\author{
Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.

This function is based on a routine used to describe a simulation setting in
Bremhorst and Lambert (2016). Special thanks go to Vincent Bremhorst who
shared this routine during his PhD thesis.
}
