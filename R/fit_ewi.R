#' Fit a Bayesian EWI model
#'
#' @param data_frame Data frame containing the time covariarte ('x') and response ('y') with no NAs
#' @param model Which model to fit. Either the dynamic AR(1) model ('ar') or stochastic volatility model with time varying volatility ('sv'). Defaults to 'ar'.
#' @param time Whether to fit model with discrete time ('discrete') or continuous time ('continuous'). Defaults to 'discrete'.
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
#' model_1 = fit_ewi(data_example, ewi_model="ar", time="continuous")
#' model_2 = fit_ewi(data_example, ewi_model="ar", time="discrete", iter = 1000, chains=2)
#'
#' print(plot_estimates(model_2))
#' model_3 = fit_ewi(data_example, ewi_model="sv", time="discrete", iter = 2000, chains=3)
#' }
fit_ewi <- function(data_frame, ewi_model = "ar", time = "discrete", iter = 2000,
  chains = 4, control = list(adapt_delta = 0.9, max_treedepth = 20), ar_scale = 0.1) {
  #if(model %in% c("ar","sd") == FALSE)

  if(ewi_model == "ar" & time == "discrete") {
    # fit ar model
    datalist = list(x = data_frame$x,
      y = data_frame$y,
      N = nrow(data_frame),
      obs_model = 1,
      maxt = max(data_frame$x))
    pars <- c("ar", "ar_sd", "drift", "CV", "ar_logit")
    mod = stan(file = "exec/ar_model_discrete.stan",
      data = datalist,
      pars = pars,
      iter=iter,
      chains = chains,
      control = control)
  }
  if(ewi_model == "ar" & time == "continuous") {
    # difference with continuous model is that x is differenced.
    # time also passed in
    data_frame$time = match(data_frame$x, unique(data_frame$x))
    datalist = list(N = nrow(data_frame),
      y = data_frame$y,
      maxt = length(c(0, diff(unique(data_frame$x)))),
      x = c(0, diff(unique(data_frame$x))),
      time = data_frame$time,
      obs_model = 1,
      ar_scale = ar_scale)
    pars <- c("ar", "ar_sd", "drift", "CV", "ar_logit")
    mod = stan(file = "exec/ar_model_continuous.stan",
      data = datalist,
      pars = pars,
      iter=iter,
      chains = chains,
      control = control)
  }
  if(ewi_model == "sv" & time == "discrete") {
    # fit ar model
    datalist = list(x = data_frame$x,
      y = data_frame$y,
      N = nrow(data_frame),
      obs_model = 1,
      maxt = max(data_frame$x))
    pars <- c("log_CV", "states", "sigma_obs", "CV_sd", "drift")
    mod = stan(file = "exec/sv_model_discrete.stan",
      data = datalist,
      pars = pars,
      iter=iter,
      chains = chains,
      control = control)
  }

  return(list("model"=mod, "data"=datalist, "time"=time, "ewi_model" = ewi_model))
}
