#' Specification of covariates entering the short-term part
#' in a promotion time cure model.
#' @keywords internal
#' @export


# Covariates entering the short-term survival (Cox part)
st <- function(...) {
  st.covariates <- as.list(substitute(list(...)))[-1]
  if (length(st.covariates) == 1) {
    gsub(" + ", "+", st.covariates, fixed = TRUE)
  } else{
    paste(unlist(st.covariates),  collapse = "+")
  }
}
