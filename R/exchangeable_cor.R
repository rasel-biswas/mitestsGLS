#' Gives the compound symmetric correlation matrix
#'
#'  (This is an internal function).
#'
#' @return A compound symmetric correlation matrix of order n.
#'
#' @import mice
#' @export

exchangeable_cor <- function(n, rho) {
  mat <- matrix(rho, n, n)
  diag(mat) <- 1
  mat
}
