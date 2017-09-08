#' Plot results from a Bayesian EWI model
#'
#' @param fitted Object returned from fit_ewi
#' @param alpha_ci alpha level for credible intervals, defaults to 0.05 = 95% CIs
#'
#' @import ggplot2
#' @importFrom broom tidy
#' @importFrom rstan extract
#' @export
#' @examples
#' \dontrun{
#' library(bayesewi)
#' model_1 = fit_ewi(data_example, ewi_model="ar", time="continuous")

#' print(plot_estimates(model_1, alpha_ci = 0.1))
#' }
plot_estimates = function(fitted, alpha_ci=0.05) {
  pars = broom::tidy(fitted$model) # note that the tidy estimates for SEs aren't same as CIs
  ewi_model = fitted$ewi_model
  n_t = fitted$data$maxt

  pars = pars[startsWith(pars$term, paste0(ewi_model,"[")), ]

  pars_mcmc = extract(fitted$model, pars = ewi_model)[[1]]
  pars$lo = apply(pars_mcmc, 2, quantile, alpha_ci/2)
  pars$hi = apply(pars_mcmc, 2, quantile, 1-alpha_ci/2)
  pars$t = seq(1, nrow(pars))
  g1 = ggplot(pars, aes(t, estimate)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill="blue", alpha=0.4) +
    geom_line(color="blue") + xlab("Time") + ylab("Estimated parameter")
  return(g1)
}
