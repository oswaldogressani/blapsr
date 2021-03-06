#' Melanoma survival data.
#'
#' @docType data
#'
#' @description Melanoma survival dataset with 205 patients suffering from
#'  skin cancer and operated for malignant melanoma at Odense University
#'  Hospital in Denmark.
#'
#' @usage data(melanoma)
#'
#' @format A data frame with 205 rows and 7 columns.
#' \describe{
#'  \item{\code{time}}{Survival time in years.}
#'  \item{\code{status}}{\code{1} Died from melanoma, \code{0} still alive or died from
#'        another event.}
#'  \item{\code{sex}}{\code{1}=Male, \code{0}=Female.}
#'  \item{\code{age}}{Age in years.}
#'  \item{\code{year}}{Year of operation.}
#'  \item{\code{thickness}}{Tumour thickness measured in mm.}
#'  \item{\code{ulcer}}{\code{1}=Presence of ulceration, \code{0}=Absence of
#'              ulceration.}
#'
#' }
#'
#' @source \url{http://www.stats.ox.ac.uk/pub/MASS4/}
#'
#' @references  Venables W.N., and Ripley, B.D. (2002). Modern Applied
#'  Statistics with S, Fourth edition. Springer, New York. ISBN 0-387-95457-0.
#' @references Andersen, P.K., Borgan, O., Gill, R.D., and Keiding, N. (1993)
#'  Statistical Models based on Counting Processes. Springer.
"melanoma"
