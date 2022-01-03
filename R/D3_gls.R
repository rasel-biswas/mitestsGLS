#' Compare two nested models using MI LRT
#'
#' The MI LRT (D3) test by Meng and Rubin (1992).
#'
#' @references
#' Meng, X. L., and D. B. Rubin. 1992.
#' Performing Likelihood Ratio Tests with Multiply-Imputed Data Sets.
#' \emph{Biometrika}, 79 (1): 103-11.
#'
#' @import mitml nlme
#' @export

D3_gls <- function(model, null_model, datalist) {
  m <- length(model)
  # evaluate log likelihood
  l1 <- map_dbl(model, function(x) logLik(x))
  l0 <- map_dbl(null_model, function(x) logLik(x))
  # pooled LR
  d_bar <- -2 * mean(l0 - l1) # mean of the m LR test statistics
  # pool the parameter estimates, psi = (beta, sigma, rho)
  psi_bar1 <- get_pooled_estimates(model)
  psi_bar0 <- get_pooled_estimates(null_model)
  # re-evaluate log likelihood at pooled parameter estimates
  for (i in 1:m) {
    model[[i]]$dat <- datalist[[i]]
    null_model[[i]]$dat <- datalist[[i]]
  }
  l1_pooled <- map_dbl(model, function(x)
    likelihood_at_fixed_psi(x, psi = psi_bar1))
  l0_pooled <- map_dbl(null_model, function(x)
    likelihood_at_fixed_psi(x, psi = psi_bar0))
  # pooled LR, based on pooled estimates
  d_tilde <- -2 * mean(l0_pooled - l1_pooled)
  # ARIV
  k <- length(psi_bar1$beta) - length(psi_bar0$beta)
  r <- ((m + 1) / (k * (m - 1))) * (d_bar - d_tilde)
  # D3
  D <- d_tilde / (k * (1 + r))
  v <- k * (m - 1)
  if (v > 4) {
    w <- 4 + (v - 4) * ((1 + (1 - 2 / v) * (1 / r)) ^ 2)
  } else {
    w <- v * (1 + 1 / k) * ((1 + 1 / r) ^ 2) / 2
  }
  pvalue = ifelse(D>=0 & r>=0, 1 - pf(D, k, w), NA)
  out <-
    data.frame(
      F.statistic = D,
      df1 = k,
      df2 = w,
      p.value = pvalue,
      RIV = r
    )
  out
}
