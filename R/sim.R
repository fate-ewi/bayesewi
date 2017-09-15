#' @export
sim_sv <- function(time, sigma2_mu, sigma_obs, gp_scale, gp_sigma_sq, state0 = 1, sigma0 = 0) {

  maxt <- length(time)
  states <- vector(mode = "numeric", length = length(time))
  sigma <- vector(mode = "numeric", length = length(time))
  Sigma <- matrix(nrow = length(time)-1, ncol = length(time)-1)
  states[1] <- state0
  sigma[1] <- 0
  pro_dev <- rnorm(length(time), 0, 1)


  for(i in 1:(maxt-1)) {
    for(j in 1:(maxt-1)) {
      Sigma[i,j] <- gp_sigma_sq * exp(-(time[i] - time[j])^2 / gp_scale);
    }
  }
  log_sigma2 <- mvtnorm::rmvnorm(1, mean = sigma2_mu, sigma = Sigma)

  for(i in 2:maxt) {
    sigma[i] <- sqrt(exp(log_sigma2[i-1]))
    states[i] <- states[i-1] * exp((0 - sigma[i]*sigma[i]/2) * time[i] + sigma[i]* pro_dev[i-1]);
  }
  y <- rnorm(length(states), states, sigma_obs);
  y
}
