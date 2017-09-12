#' Plot results from a Bayesian EWI model
#'
#' @param fitted Object returned from fit_ewi
#' @param alpha_ci alpha level for credible intervals, defaults to 0.05 = 95% CIs
#' @param link_space Whether to plot the parameters in link or normal space. Defaults to TRUE (link space)
#'
#' @import ggplot2
#' @importFrom rstan extract
#' @importFrom broom tidy
#' @importFrom gridExtra grid.arrange
#'
#' @export
#' @examples
#' \dontrun{
#' library(bayesewi)
#' model_1 = fit_ewi(data_example, ewi_model="sv", iter = 1000, chains=3)

#' print(plot_estimates(model_1, alpha_ci = 0.1)) # make plot in log space
#' print(plot_estimates(model_1, alpha_ci = 0.1, link_space=FALSE)) # make plot in normal space
#' }
plot_estimates = function(fitted,
  alpha_ci = 0.05,
  link_space = TRUE) {
  pars = tidy(fitted$model) # note that the tidy estimates for SEs aren't same as CIs
  ewi_model = fitted$ewi_model
  n_t = fitted$data$maxt

  if (ewi_model != "sv") {
    pars = pars[startsWith(pars$term, "ar["),]
    if (link_space == FALSE) {
      pars_mcmc = extract(fitted$model, pars = "ar")[[1]]
    } else {
      pars_mcmc = extract(fitted$model, pars = "ar")[[1]]
      pars_mcmc = log(pars_mcmc) / (1 - log(pars_mcmc)) # inv-logit
    }
  } else {
    pars = pars[startsWith(pars$term, "sigma["),]
    if (link_space == FALSE) {
      pars_mcmc = extract(fitted$model, pars = "sigma")[[1]]
    } else {
      pars_mcmc = extract(fitted$model, pars = "sigma")[[1]]
      pars_mcmc = log(pars_mcmc)
    }
  }

  pars$estimate = apply(pars_mcmc, 2, median, na.rm = T)
  pars$lo = apply(pars_mcmc, 2, quantile, alpha_ci / 2, na.rm = T)
  pars$hi = apply(pars_mcmc, 2, quantile, 1 - alpha_ci / 2, na.rm = T)

  # not actually estimated in STAN
  pars$estimate[1] = NA
  pars$lo[1] = NA
  pars$hi[1] = NA

  pars$t = seq(1, nrow(pars))

  df = data.frame(
    "Time" = fitted$data$uniquet[fitted$data$x],
    "Obs" = fitted$data$y,
    "Pred" = apply(extract(fitted$model, pars = "states")[[1]], 2, median)[fitted$data$x]
  )
  g1 = ggplot(df, aes(Time, Pred)) +
    geom_line() + geom_point(aes(Time, Obs), col = "red", alpha = 0.5) +
    xlab("Time") + ylab("Observed and predicted")

  df$estimate = pars$estimate[fitted$data$x]
  df$lo = pars$lo[fitted$data$x]
  df$hi = pars$hi[fitted$data$x]
  g2 = ggplot(df, aes(Time, estimate)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "blue", alpha = 0.4) +
    geom_line(color = "blue") + xlab("Time") + ylab("Estimated parameter")

  return(grid.arrange(g1, g2, ncol = 1))
}
