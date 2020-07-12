#' Specification of covariates entering the long-term part
#' in a promotion time cure model.
#' @keywords internal
#' @export

# Covariates entering the long-term survival (proba to be cured)
lt <- function(...) {
  lt.covariates <- as.list(substitute(list(...)))[-1]
  if (length(lt.covariates) == 1) {
    gsub(" + ", "+", lt.covariates, fixed = TRUE)
  } else{
    paste(unlist(lt.covariates),  collapse = "+")
  }
}

