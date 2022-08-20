#' Simulation of survival times for the promotion time cure model.
#'
#' Generates right censored time-to-event data with a plateau in the
#' Kaplan-Meier estimate.
#'
#' @usage simcuredata(n, censor = c("Uniform", "Weibull"), cure.setting = 1,
#'             info = TRUE, KapMeier = FALSE)
#'
#' @param n Sample size.
#' @param censor The censoring scheme. Either Uniform (the default) or Weibull.
#' @param cure.setting A number indicating the desired cure percentage. If
#'  \code{cure.setting = 1} (default) the cure percentage is around 20\%.
#'  With  \code{cure.setting = 2} the cure percentage is around 30\%.
#' @param info Should information regarding the simulation setting be printed
#'  to the console? Default is \code{TRUE}.
#' @param KapMeier Logical. Should the Kaplan-Meier curve of the generated
#'  data be plotted? Default is \code{FALSE}.
#'
#' @details Latent event times are generated following Bender et al. (2005),
#'  with a baseline distribution chosen to be a Weibull with mean 8 and variance
#'  17.47. When \code{cure.setting = 1} the regression coefficients of the
#'  long-term survival part are chosen to yield a cure percentage around 20\%,
#'  while \code{cure.setting = 2} yields a cure percentage around 30\%.
#'  Censoring is either governed by a Uniform distribution on the support
#'  [20, 25] or by a Weibull distribution with shape parameter 3 and
#'  scale parameter 25.
#'
#' @return A list with the following components:
#'
#' \item{n}{Sample size.}
#'
#' \item{survdata}{A data frame containing the simulated data.}
#'
#' \item{beta.coeff}{The regression coefficients pertaining to long-term
#'  survival.}
#'
#' \item{gamma.coeff}{The regression coefficients pertaining to short-term
#'  survival.}
#'
#' \item{cure.perc}{The cure percentage.}
#'
#' \item{censor.perc}{The percentage of censoring.}
#'
#' \item{censor}{The censoring scheme.}
#'
#' \item{S0}{The baseline survival function under the chosen Weibull
#'   parameterization.}
#'
#' @examples
#' set.seed(10)
#' sim <- simcuredata(n = 300, censor = "Weibull", KapMeier = TRUE)
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' This function is based on a routine used to describe a simulation setting in
#' Bremhorst and Lambert (2016). Special thanks go to Vincent Bremhorst who
#' shared this routine during his PhD thesis.
#'
#' @references Bender, R., Augustin, T. and Blettner, M. (2005). Generating
#'  survival times to simulate Cox proportional hazards models,
#'  \emph{Statistics in Medicine} \strong{24}(11): 1713-1723.
#' @references Bremhorst, V. and Lambert, P. (2016). Flexible estimation in
#'   cure survival models using Bayesian P-splines. \emph{Computational
#'   Statistics & Data Analysis} \strong{93}: 270-284.
#' @references Gressani, O. and Lambert, P. (2018). Fast Bayesian inference
#'   using Laplace approximations in a flexible promotion time cure model based
#'   on P-splines. \emph{Computational Statistics & Data Analysis} \strong{124}:
#'   151-167.
#' @export

simcuredata <- function(n, censor = c("Uniform", "Weibull"),
                  cure.setting = 1, info = TRUE, KapMeier = FALSE){

  if(n <= 1) stop("Sample size n must be larger than one.")
  if(floor(cure.setting) < 1 || floor(cure.setting) > 2){
    stop("Cure setting must be either 1 or 2.")
  }
  rcens <- 23   # Upper bound for observed failure time
  xtrunk <- 25  # Truncation value for the Weibull distribution )
  a <- 2
  b <- (1 / exp(-4.4)) ^ (0.5)
  Weibull.mean <- round(b * gamma(1 + (1 / a)), 4)
  Weibull.variance <- round((b ^ 2) * (gamma(1 + 2 / a) -
                                         (gamma(1 + 1 / a) ^ 2)), 4)
  S0 <- function(t) exp( - (t / b) ^ a)     # Baseline survival
  continuous <- stats::rnorm(n, 0, 1)              # Normal variate
  binary <- stats::rbinom(n, 1, prob = 0.5) - 0.5  # Bernoulli variate
  X <- cbind(Intercept = 1, continuous, binary)
  Z <- cbind(continuous, binary)
  gammas <- c(0.4, -0.4)
  if(cure.setting == 1) betas <- c(0.75, 0.80, -0.50)  # cure 20%
  if(cure.setting == 2) betas <- c(0.30, 1, -0.75)     # cure 30%
  theta.poiss <- as.numeric(exp(X %*% betas))          # Mean of Poisson
  Cox <- as.numeric(exp(Z %*% gammas))                 # Cox part
  Nlatent <- c()                                       # Latent event times
  tobs <- c()                                          # Observed event times
  delta <- rep(1, n)                                   # Event indicator

  # Generation of latent event times
  for (i in 1:n) {
    Nlatent[i] <- stats::rpois(1, theta.poiss[i])
    if (Nlatent[i] == 0) {
      tobs[i] <- 99 #cured subject
    } else{
      #non-cured subject
      u <- stats::runif(Nlatent[i], 0, 1)
      t <- b * (-(log(1 - u)) / (Cox[i])) ^ (1 / a)
      tobs[i] <- min(t)
      while (tobs[i] > rcens) {
        u <- stats::runif(Nlatent[i], 0, 1)
        t <- b * (-(log(1 - u)) / (Cox[i])) ^ (1 / a)
        tobs[i] <- min(t)
      }
    }
  }

  censortype <- match.arg(censor)
  if (censortype == "Uniform")
    censor <- "Uniform"
  else if (censortype == "Weibull")
    censor <- "Weibull"
  else
    stop("Censoring distribution must be either Uniform or Weibull")

  if(censor == "Uniform") {
    C <- stats::runif(n, 20, 25)
  } else if(censor == "Weibull"){
    a.censor <- 3   # shape
    b.censor <- 25  # scale
    C <- stats::rweibull(n, a.censor, b.censor)
    Cbig <- which(C > xtrunk)
    while (length(Cbig) > 0) {
      C[Cbig] <- stats::rweibull(length(Cbig), a.censor, b.censor)
      Cbig <- which(C > xtrunk)
    }
  }
  censored <- which(tobs > C)
  tobs[censored] <- C[censored]
  delta[censored] <- 0

  cure.perc <- round(sum(Nlatent == 0) / n * 100, 2)
  censor.perc <- round((sum(1 - delta)) / n * 100, 2)

  if(info == TRUE){
    cat("Sample size:             ", n, "\n")
    cat("Censor setting:          ", censor, "\n")
    cat("Cure percentage is:      ", cure.perc, "%", "\n")
    cat("Censoring percentage is: ", censor.perc, "%")
  }

  survdata <- data.frame(cbind(time = tobs, status = delta, continuous, binary))

  if(KapMeier == TRUE){
    KaplanMeier <- survival::survfit(survival::Surv(time,status) ~ 1,
                                     data = survdata)
    graphics::plot(KaplanMeier, mark.time = TRUE, mark = 4, xlab = "Time")
  }


  listout <- list(n = n,
                  survdata = survdata,
                  beta.coeff = betas,
                  gamma.coeff = gammas,
                  cure.perc = cure.perc,
                  censor.perc = censor.perc,
                  censor = censor,
                  S0 = S0)

  return(invisible(listout))
}
