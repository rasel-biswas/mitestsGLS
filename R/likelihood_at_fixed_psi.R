#' Computes likelihood at fixed value of parameters.
#'
#'  (This is an internal function).
#'
#' @return The likelihood value at a given psi.
#'
#' @import mvtnorm
#' @export

likelihood_at_fixed_psi <- function(object, psi) {
  b <- matrix(data= psi$beta, ncol=1)
  sigma <- as.numeric(psi$sigma)
  rho <- psi$rho
  n <- dim(getVarCov(object))[1]
  N <- object$dims$N / n
  if (isTRUE(object$call$correlation[1]=="corAR1()"))
    corr_mat <- ar1_cor(n = n, rho = rho)
  if (isTRUE(object$call$correlation[1]=="corCompSymm()"))
    corr_mat <- exchangeable_cor(n = n, rho = rho)
  Sigma_hat_bar <- (sigma^2) * corr_mat
  Z <- object$dat
  y_term <- paste(object$call$model[[2]])
  id_term <- paste(object$call$correlation$form[[2]][[3]])
  formula <- as.formula(object$call[["model"]])
  ll = 0
  unique_id <- unique(as.numeric(Z[, id_term]))
  for (i in unique_id) {
    Xi <- as.matrix(model.matrix(formula,
                                 data = Z[Z[, id_term]==i, ])[1:n, 1:(length(b))])
    yi <- Z[, y_term][Z[, id_term] == i]
    mu <- (Xi %*% b)
    l <- mvtnorm::dmvnorm(x = yi, mean = mu, sigma=Sigma_hat_bar, log=TRUE)
    ll <- ll + l
  }
  return(as.numeric(ll))
}
