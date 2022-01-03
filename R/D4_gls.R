#' Compare two nested models using D4 test
#'
#' This is the improved MI LRT by Chan and Meng (2019).
#'
#'
#' @import mice
#' @export

D4_gls <- function(model, null_model, datalist) {
  m <- length(model)
  # evaluate log likelihood
  l1 <- map_dbl(model, function(x) logLik(x))
  l0 <- map_dbl(null_model, function(x) logLik(x))
  # pooled LR
  d_bar <- -2 * mean(l0 - l1) # mean of the m LR test statistics
  # prepare the stacked data
  id_term <- paste(model[[1]]$call$correlation$form[[2]][[3]])
  for(ii in 1:m){
    datalist[[ii]][,id_term] <-
      paste0("imp", ii, "_", datalist[[ii]][,id_term])
  }
  stacked_data <- do.call(rbind, datalist)
  stacked_data[,id_term] <- as.integer(as.factor(stacked_data[,id_term]))
  # re-evaluate log-likelihood in stacked data
  model_full <- update(model[[1]], data = stacked_data)
  model_null <- update(null_model[[1]], data = stacked_data)
  l1_s <- logLik(model_full)
  l0_s <- logLik(model_null)
  # pooled LR in stacked data
  dhat_s <- -2 * as.numeric(l0_s - l1_s) / m
  # ARIV
  k <- length(model_full$coefficients) - length(model_null$coefficients)
  r <- ((m + 1) / (k * (m - 1))) * (d_bar - dhat_s)
  r <- max(0, r)
  # D4
  D <- dhat_s / (k * (1 + r))
  v <- (k*(m-1)) * (1 + r^(-1))^2
  p_value <- pf(D, k, v, lower.tail = FALSE)
  out <-
    data.frame(
      Dm.statistic = D,
      df1 = k,
      df2 = v,
      p.value = p_value,
      RIV = r
    )
  out
}
