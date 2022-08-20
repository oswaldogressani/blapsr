#' Simulation of data for (Generalized) additive models.
#'
#' Simulation of a data set that can be used to illustrate the
#' \code{\link{amlps}} or \code{\link{gamlps}} routines to fit
#' (generalized) additive models with the Laplace-P-spline methodology.
#'
#' @param setting The simulation setting. The default is \code{setting = 1}
#'  for a setting with three smooth terms, while \code{setting = 2} is another
#'  setting with only two smooth terms. The coefficients of the linear part
#'  of the predictor are also different in the two settings.
#' @param n The sample size to simulate.
#' @param dist A character string to specify the response distribution. The
#'  default is \code{"gaussian"}. Other distributions can be \code{"poisson"},
#'  \code{"bernoulli"} and \code{"binomial"}.
#' @param scale Used to tune the noise level for Gaussian and Poisson
#'  distributions.
#' @param info Should information regarding the simulation be printed?
#' Default is true.
#'
#' @details The simulation settings contain two covariates in the linear part
#'  of the predictor, namely \emph{z1 ~ Bern(0.5)} and \emph{z2 ~ N(0,1)}. The
#'  smooth additive terms are inspired from Antoniadis et al. (2012).
#'  For Binomial data, the number of trials is fixed to 15.
#'
#' @return An object of class \emph{simgam}. Plot of a \emph{simgam} object
#'  yields a scatter plot of the generated response values.
#'
#' \item{data}{A data frame.}
#'
#' \item{f}{The true smooth functions.}
#'
#' \item{betas}{The regression coefficients of the linear part. The first
#'  term is the intercept.}
#'
#' \item{dist}{The distribution of the response.}
#'
#' @examples
#' set.seed(10)
#' sim <- simgamdata(n = 150, dist = "poisson", scale = 0.3)
#' plot(sim)
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @references Antoniadis, A., Gijbels, I., and Verhasselt, A. (2012). Variable
#'  selection in additive models using P-splines. \emph{Technometrics}
#'  \strong{54}(4): 425-438.
#'
#' @export

simgamdata <- function(setting = 1, n, dist = "gaussian", scale = 0.5,
                       info = TRUE){

  if(setting == 1) {
    f1 <- function(x) 0.5 * (2 * x - 1) ^ 2
    f2 <- function(x) cos (2 * pi * x)
    f3 <- function(x) 4 * (sin(2 * pi * x)) / (2 - sin (2 * pi * x))
    f <- list(f1, f2, f3)
    x1 <- stats::runif(n, -1, 1)
    x1[sample(1:n, size = 2)] <- c(-1, 1)
    x2 <- stats::runif(n, -1, 1)
    x2[sample(1:n, size = 2)] <- c(-1, 1)
    x3 <- stats::runif(n, -1, 1)
    x3[sample(1:n, size = 2)] <- c(-1, 1)
    z1 <- stats::rbinom(n, 1, prob = 0.5)
    z2 <- stats::rnorm(n, 0, 1)
    beta0 <- -1.45
    beta1 <-  0.25
    beta2 <- -0.90
    predictor <- beta0 + beta1 * z1 + beta2 * z2 + f1(x1) + f2(x2) + f3(x3)

    if (dist == "poisson") {
      mu <- exp(predictor * scale)
      y <- stats::rpois(n, mu)
      beta0 <- scale * beta0
      beta1 <- scale * beta1
      beta2 <- scale * beta2
    } else if (dist == "gaussian") {
      e <- stats::rnorm(n, 0, scale)
      mu <- predictor
      y <- mu + e
    } else if (dist == "bernoulli") {
      mu <- exp(predictor) / (1 + exp(predictor))
      y <- stats::rbinom(n, 1, prob = mu)
    } else if (dist == "binomial") {
      mu <- exp(predictor) / (1 + exp(predictor))
      y <- stats::rbinom(n, 15, prob = mu)
    } else stop("unknown distribution")

    if (info == TRUE) {
      cat("Setting      :", setting, "\n")
      cat("Sample size n:",  n, "\n")
      cat("Distribution :", dist, "\n")
      cat("------------------------- \n")
      cat("Covariates generated: \n")
      cat("z1 ~ Bern(0.5) \n")
      cat("z2 ~ N(0,1) \n")
      cat("xj ~ U(-1,1), j = 1,2,3 \n")
      cat("True linear coefficients:"  , c(beta0, beta1, beta2), "\n")
    }

    data <- data.frame(y = y, z1 = z1, z2 = z2, x1 = x1, x2 = x2, x3 = x3)
    listout <- list(data = data, f = f, betas = c(beta0, beta1, beta2),
                    dist = dist)
    attr(listout, "class") <- "simgam"
    return(invisible(listout))
  } else if (setting == 2) {
    f1 <- function(x) cos(2 * pi * x)
    f2 <- function(x) 3 * (x ^ 5) + 2 * sin(4 * x) + 1.5 * (x ^ 2) - 0.5
    f <- list(f1, f2)
    x1 <- stats::runif(n, -1, 1)
    x1[sample(1:n, size = 2)] <- c(-1, 1)
    x2 <- stats::runif(n, -1, 1)
    x2[sample(1:n, size = 2)] <- c(-1, 1)
    z1 <- stats::rbinom(n, 1, prob = 0.5)
    z2 <- stats::rnorm(n, 0, 1)
    beta0 <- 0.20
    beta1 <- 0.85
    beta2 <- -0.45
    predictor <- beta0 + beta1 * z1 + beta2 * z2 + f1(x1) + f2(x2)

    if (dist == "poisson") {
      mu <- exp(predictor * scale)
      y <- stats::rpois(n, mu)
      beta0 <- scale * beta0
      beta1 <- scale * beta1
      beta2 <- scale * beta2
    } else if (dist == "gaussian") {
      e <- stats::rnorm(n, 0, scale)
      mu <- predictor
      y <- mu + e
    } else if (dist == "bernoulli") {
      mu <- exp(predictor) / (1 + exp(predictor))
      y <- stats::rbinom(n, 1, prob = mu)
    } else if (dist == "binomial") {
      mu <- exp(predictor) / (1 + exp(predictor))
      y <- stats::rbinom(n, 15, prob = mu)
    } else stop("unknown distribution")

    if (info == TRUE) {
      cat("Setting      :", setting, "\n")
      cat("Sample size n:",  n, "\n")
      cat("Distribution :", dist, "\n")
      cat("------------------------- \n")
      cat("Covariates generated: \n")
      cat("z1 ~ Bern(0.5) \n")
      cat("z2 ~ N(0,1) \n")
      cat("xj ~ U(-1,1), j = 1,2 \n")
      cat("True linear coefficients:"  , c(beta0, beta1, beta2), "\n")
    }

    data <- data.frame(y = y, z1 = z1, z2 = z2, x1 = x1, x2 = x2)
    listout <- list(data = data, f = f, betas = c(beta0, beta1, beta2),
                    dist = dist)
    attr(listout, "class") <- "simgam"
    return(invisible(listout))
  } else stop("Wrong setting")
}
