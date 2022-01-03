#' Gives the AR1 correlation matrix
#'
#'  (This is an internal function).
#'
#' @return A nxn AR1 correlation matrix
#'
#' @import
#' @export

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(0:(n-1), nrow = n, ncol = n, byrow = TRUE) -
                    0:(n-1))
  rho^exponent
}
