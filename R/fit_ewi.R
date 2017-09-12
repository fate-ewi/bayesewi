#' Fit Bayesian EWI model
#'
#' @param data_frame Data frame containing the time covariarte ('x') and response ('y') with no NAs
#' @param ewi_model Which model to fit. Either the dynamic AR(1) model ('ar') or stochastic volatility model ('sv'). Defaults to 'ar'. For both cases, a gaussian process model is fit to the time-varying parameters.
#' @param iter Number of iterations in Stan sampling.
#' @param chains Number of chains in Stan sampling.
#' @param control A list of options to pass to Stan sampling.
#' @export
#'
#' @importFrom rstan sampling
#' @import Rcpp
#' @examples
#' \dontrun{
#' library(bayesewi)
#' model_1 = fit_ewi(data_example, ewi_model="ar")
#' model_2 = fit_ewi(data_example, ewi_model="ar", iter = 1000, chains=2)
#'
#' print(plot_estimates(model_2))
#' options(mc.cores = parallel::detectCores())
#' model_2 = fit_ewi(data_example, ewi_model="ar", iter = 100, chains=3)
#' }
fit_ewi <- function(data_frame,
  ewi_model = "ar",
  iter = 2000,
  chains = 4,
  control = list(adapt_delta = 0.9, max_treedepth = 20)) {
  # data list input to stan
  datalist = list(
    N = nrow(data_frame),
    y = data_frame$y,
    x = as.integer(as.numeric(as.factor(data_frame$x))),
    maxt = length(c(0, diff(
      unique(data_frame$x)
    ))),
    uniquet = unique(data_frame$x) - data_frame$x[1],
    deltat = c(0, diff(unique(data_frame$x))),
    obs_model = 1
  )
  # parameters
  pars <- c("states",
    "gp_scale",
    "gp_sigma_sq",
    "sigma_obs")
  if (ewi_model == "sv") {
    pars = c(pars, "sigma2_mu", "sigma")
    model_file = "exec/sv_model_gp.stan"
  } else {
    pars = c(pars, "ar_mu", "ar")
    model_file = "exec/ar_model_gp.stan"
  }

  mod = stan(
    file = model_file,
    data = datalist,
    pars = pars,
    iter = iter,
    chains = chains,
    control = control
  )

  return(list(
    "model" = mod,
    "data" = datalist,
    "time" = time,
    "ewi_model" = ewi_model
  ))
}
