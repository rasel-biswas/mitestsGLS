#' Gives the correlation value for compound symmetry correlation structure
#'
#'  (This is an internal function).
#'
#' @return A correlation value
#'
#' @import mice
#' @export

get_rho <- function(object) {
  phi <- coef(object$modelStruct$corStruct, unconstrained = FALSE)
  as.numeric(phi)
}
