#' Simulation of right censored survival times for the Cox model.
#'
#' Generates right censored time-to-event data. Latent event times are drawn
#' from a Weibull distribution, while censoring times are generated from an
#' exponential distribution.
#'
#' @param a,b The shape parameter `a>0` and scale parameter `b>0`
#' of the Weibull.
#' @param n Sample size.
#' @param betas A numeric vector of regression coefficients. Allowed components
#'  of `betas` are in the interval \emph{[-1 ,1]} and the total
#'  number of components cannot exceed 5.
#' @param censperc A numeric value in \emph{[0,100]} corresponding to
#'        the targeted percentage of censoring.
#' @param tmax A maximum upper bound for the generated latent event times.
#'  Especially useful for a simulation study in which the observed event times
#'  are constrained to be generated in a fixed range.
#'
#'
#' @details The Weibull baseline hazard is parameterized as follows (see Hamada
#'   et al. 2008 pp. 408-409) :
#'   \deqn{h_0(t) = (a/(b^a)) t^(a-1), t > 0.}
#'   The i\emph{th} latent event time is denoted by \emph{T_i} and is generated
#'   following Bender et al. (2005) as follows:
#'   \deqn{T_i = b (-log(U_i) exp(-\beta^T x_i))^(1/a),}
#'   where \emph{U_i} is a uniform random variable obtained with `runif(1)`
#'   , \emph{x_i} is the i\emph{th} row of a covariate matrix X of dimension
#'   `c(n, length(betas))` where each component is generated from a
#'   standard Gaussian distribution and \eqn{\beta} is the vector of
#'   regression coefficients given by `betas`.
#'
#' @return An object of class `simsurvdata` which is a list
#' with the following components:
#'
#' \item{sample.size}{Sample size.}
#'
#' \item{censoring}{Censoring scheme. Either \emph{No censoring}
#'   or \emph{Exponential}.}
#'
#' \item{num.events}{Number of events.}
#'
#' \item{censoring.percentage}{The effective censoring percentage.}
#'
#' \item{survdata}{A data frame containing the simulated data.}
#'
#' \item{regcoeffs}{The true regression coefficients used to simulate
#'   the data.}
#'
#' \item{S0}{The baseline survival function under the chosen Weibull
#'   parameterization.}
#'
#' \item{h0}{The baseline hazard function under the chosen Weibull
#'   parameterization.}
#'
#' \item{Weibull.mean}{The mean of the Weibull used to generate latent
#'   event times.}
#'
#'  \item{Weibull.variance}{The variance of the Weibull used to generate latent
#'   event times.}
#'
#' The `print` method summarizes the generated right censored data and
#' the `plot` method produces a graph with time on the x axis and
#' horizontal bars on the y axis corresponding either to an event or a
#' right censored observation. If `n > 25`, only the 25 first observations
#' are plotted.
#'
#' @examples
#' set.seed(10)
#' sim <- simsurvdata(a = 2, b = 1, n = 300, betas = c(0.8, -0.6), censperc = 25)
#' sim
#' plot(sim)
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @references Bender, R., Augustin, T. and Blettner, M. (2005). Generating
#'             survival times to simulate Cox proportional hazards models,
#'             \emph{Statistics in Medicine} \strong{24}(11): 1713-1723.
#' @references Hamada, M. S., Wilson, A., Reese, C. S. and Martz, H. (2008).
#'            \emph{Bayesian Reliability}. Springer Science and Business Media.
#' @export


simsurvdata <- function(a, b, n, betas, censperc, tmax = NULL) {
  if (!is.vector(a, mode = "numeric") || !is.vector(b, mode = "numeric"))
    stop("a, b must be numeric of length 1")
  if (length(a) > 1 || length(b) > 1)
    stop("a, b must be numeric of length 1")
  if (a <= 0 || b <= 0 || is.infinite(a) || is.infinite(b))
    stop("a, b must be postive and finite")
  if (!is.vector(betas, mode = "numeric"))
    stop("The vector of regression coefficients must be numeric")
  if(any(betas < -1) || any(betas > 1))
    stop("Components in betas must be between -1 and 1")
  if(length(betas) > 5)
    stop("betas must be a vector of maximum length 5")
  if (!is.vector(censperc, mode = "numeric"))
    stop("Censoring percentage must be a number between 0 and 100")
  if (length(censperc) > 1)
    stop("Censoring percentage must be of length 1")
  if (censperc < 0 || censperc > 100)
    stop("Censoring percentage must be a number between 0 and 100")

  n <- as.integer(n)
  if(n < 1)
    stop("Sample size must be at least 1")
  Weibull.mean <- b * gamma(1 + (1 / a))
  Weibull.variance <- (b ^ 2) * (gamma(1 + 2 / a) - (gamma(1 + 1 / a) ^ 2))
  p <- length(betas)
  X <- matrix(stats::rnorm((n * p), mean = 0, sd = 1), ncol = p, nrow = n)
  alpha <- censperc / 100 # Censoring fraction
  Tlat <- as.numeric(b * (-log(stats::runif(n)) * exp(-X %*% betas)) ^ (1 / a)) # times
  if (!is.null(tmax)) {
    while (max(Tlat) > tmax) {
      num.replace <- sum(Tlat > tmax)
      index.replace <- which(Tlat > tmax)
      for (l in 1:num.replace) {
        X[index.replace[l], ] <- stats::rnorm(p, mean = 0, sd = 1)
        Tlat[index.replace[l]] <- as.numeric(b * (-log(stats::runif(1)) *
                              exp(-X[index.replace[l], ] %*% betas)) ^ (1 / a))
      }
    }
  }
  tobs <- Tlat
  delta <- rep(1, n)  # Event indicators

  # Absence of censoring
  if (censperc == 0)
    cens.message <- "No censoring"

  # Censoring follows exponential disribution
  if (censperc > 0) {
    foo.v <- function(v) n * (1 - alpha) - sum(exp(-exp(v) * Tlat))
    foo.val <- foo.v(0)
    if (foo.val == 0) {
      lambda.censor <- 1
    } else if (foo.val > 0) {
      v.left <- 0
      while (foo.val > 0) {
        v.left <- v.left - 2 # go left
        foo.val <- foo.v(v.left)
      }
      lambda.censor <- exp(stats::uniroot(foo.v, interval = c(v.left, 0))$root)
    } else{
      v.right <- 0
      while (foo.val < 0) {
        v.right <- v.right + 2 # go right
        foo.val <- foo.v(v.right)
      }
      lambda.censor <- exp(stats::uniroot(foo.v, interval = c(0, v.right))$root)
    }

    C.censor <- stats::rexp(n, rate = lambda.censor) # generate censoring times
    replace.index <- Tlat > C.censor # event times to be replaced
    tobs[replace.index] <- C.censor[replace.index] # final event times
    delta[replace.index] <- rep(0, sum(replace.index)) # final event indicators
    cens.message <- "Exponential"
  }

  S0 <- function(t) exp( - (t / b) ^ a)         # True baseline survival
  h0 <- function(t) (a / (b ^ a)) * t ^ (a - 1) # True baseline hazard
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  survdata <- data.frame(time = tobs, delta = delta, X)

  # Output list
  listout <- list(sample.size = n,
                  censoring = cens.message,
                  num.events = sum(delta),
                  censoring.percentage = paste(round(mean(1 - delta) * 100, 2),
                                               "%", sep = ""),
                  survdata = survdata,
                  regcoeffs = betas,
                  S0 = S0,
                  h0 = h0,
                  Weibull.mean = Weibull.mean,
                  Weibull.variance = Weibull.variance)

  attr(listout, "class") <- "simsurvdata"
  listout
}


























