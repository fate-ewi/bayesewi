data {
  int<lower=1> N; // number of data points
  vector[N] y; // observed data
  int<lower=1> x[N]; // time coordinates of observed data
  int<lower=1> maxt; // time steps
  int<lower=1> obs_model; // 1 = normal, 2 gamma, etc
}
parameters {
  vector[maxt-1] states_dev; // deviations / process model of ar parameter
  real<lower=0> states_sd;
  vector<lower=-1,upper=1>[maxt-1] ar_dev; // deviations / process model of ar parameter
  real<lower=0> ar_sd;
  real<lower=0> sigma_obs;
  real<lower=-1,upper=1> ar0;
  real x0;
}
transformed parameters {
  vector[maxt] ar; // parameter of ar model
  vector[maxt] states; // latent states
  vector[N] pred; // latent states
  // initial states

  states[1] = x0;
  ar[1] = ar0;
  // random walk in ar parameter
  for(i in 2:maxt) {
    // need to constrain the random walk, otherwise inf values are predicted
    ar[i] = ar[i-1] + ar_dev[i-1];
    states[i] = states[i-1] + exp(ar[i]) * states_dev[i-1];
  }
  for(i in 1:N) {
    pred[i] = states[x[i]];
  }
}
model {
  ar0 ~ student_t(3, 0, 2);
  x0 ~ student_t(3, 0, 2);
  ar_sd ~ student_t(3, 0, 2);
  ar_dev ~ normal(0, ar_sd); // process deviations
  states_sd ~ student_t(3, 0, 2);
  states_dev ~ normal(0, states_sd); // process deviations

  if(obs_model == 1) {
  sigma_obs ~ student_t(3, 0, 2);
  y ~ normal(pred, sigma_obs);
  }
}
