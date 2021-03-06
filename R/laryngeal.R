#' Survival data of male laryngeal cancer patients.
#'
#' @docType data
#'
#' @description Survival data of male patients with larynx cancer from
#'  Section 1.8 of Klein and Moeschberger (2003).
#'
#' @usage data(laryngeal)
#'
#' @format A data frame with 90 rows and 5 columns.
#' \describe{
#'  \item{\code{stage}}{Stage of disease.}
#'  \item{\code{time}}{Time to death in months.}
#'  \item{\code{age}}{Age at diagnosis of larynx cancer.}
#'  \item{\code{diagyr}}{Year of diagnosis of larynx cancer.}
#'  \item{\code{delta}}{Event indicator, \code{1}=Dead, \code{0}=Alive.}
#' }
#'
#' @source \url{https://cran.r-project.org/package=KMsurv}
#'
#' @references Klein, J.P. and Moeschberger, M. L. (2003). Survival analysis:
#'  Techniques for Censored and Truncated Data (Second edition), Springer.
#'  ISBN 978-1-4419-2985-3
"laryngeal"
