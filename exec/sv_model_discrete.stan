data {
  int<lower=1> N; // number of data points
  vector[N] y; // observed data
  int<lower=1> x[N]; // time coordinates of observed data
  int<lower=1> maxt; // time steps
  int<lower=1> obs_model; // 1 = normal, 2 gamma, etc
}
parameters {
  vector<lower=-1,upper=1>[maxt-1] pro_dev; // deviations / process model of ar parameter
  vector<lower=-1,upper=1>[maxt-1] CV_dev;
  real<lower=0> sigma_obs; // observation variance sigma
  real<lower=-1,upper=1> ar0; // fixed autocorrelation parameter
  real x0; // initial state
  real<lower=-5,upper=5> CV0; // initial sd in log space
  real drift;
  real<lower=0> CV_sd; // initial sd
}
transformed parameters {
  vector[maxt] states; // latent states
  vector[maxt] log_CV; // latent states
  vector[N] pred; // latent states
  // initial states

  states[1] = x0;
  log_CV[1] = CV0;
  // random walk in ar parameter
  for(i in 2:maxt) {
    // need to constrain the random walk, otherwise inf values are predicted
    log_CV[i] = log_CV[i-1] + CV_dev[i-1]; // sigma[i] = CV * mean
    states[i] = drift + ar0 * (states[i-1] - drift) + (states[i-1] * exp(log_CV[i])) * pro_dev[i-1];
  }
  for(i in 1:N) {
    pred[i] = states[x[i]];
  }
}
model {
  ar0 ~ student_t(3, 0, 2);
  CV0 ~ student_t(3, 0, 2);
  x0 ~ student_t(3, 0, 2);
  CV_sd ~ exponential(1);
  CV_dev ~ normal(0, CV_sd); // process deviations
  pro_dev ~ normal(0, 1); // process deviations
  drift ~ normal(0,1);

  if(obs_model == 1) {
  sigma_obs ~ student_t(3, 0, 2);
  y ~ normal(pred, sigma_obs);
  }
}
