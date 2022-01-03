#' Computes pooled estimates from the list containing `m` fitted models.
#'
#'  (This is an internal function).
#'
#' @return A list containing the completed-data estimates and
#' their variance-covariance matrices.
#'
#' @import mice
#' @export

get_pooled_estimates <- function(model){
  #1. beta
  beta_bar <- map_dfr(model, function(x) coefficients(x)) %>%
    summarise_all(.funs = mean) %>%
    as.numeric()
  #2. sigma
  sigma_bar <- map_dbl(model, function(x) sigma(x)) %>% mean()
  #3. rho
  rho_bar <- map_dbl(model, function(x) get_rho(x)) %>% mean()
  # psi = (beta, sigma, rho)
  psi <- list(beta = beta_bar, sigma = sigma_bar, rho = rho_bar)
  return(psi)
}
