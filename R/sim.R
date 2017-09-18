#' @export
# @examples
# x <- sim_sv(1:50, sigma2_mu = rep(0.5, 49), sigma_obs = 0.001,
#   gp_scale = 30, gp_sigma_sq = 0.2, state0 = 20)
# plot(x$y[-c(1:2)], type = "p")
# plot(x$states[-c(1:2)], type = "l")
# plot(x$sigma[-c(1:2)], type = "l")

sim_sv <- function(time, sigma2_mu, sigma_obs, gp_scale = 30, gp_sigma_sq = 0.2,
  state0 = 0.1, sigma0 = 0) {

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

  list(sigma = sigma, states = states, y = y, pro_dev = pro_dev)
}
