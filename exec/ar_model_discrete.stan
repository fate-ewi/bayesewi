data {
  int<lower=1> N; // number of data points
  vector[N] y; // observed data
  int<lower=1> x[N]; // time coordinates of observed data
  int<lower=1> maxt; // time steps
  int<lower=1> obs_model; // 1 = normal, 2 gamma, etc
  real ar_scale; // scale parameter for prior of random walk -- needs to keep logit < 20
}
parameters {
  vector[maxt-1] ar_dev; // deviations / process model of ar parameter
  real<lower=0> ar_sd;
  real<lower=0> CV;
  real ar0;
  real x0;
  real drift;
}
transformed parameters {
  vector[maxt] ar_logit; // parameter of ar model
  vector[maxt] ar; // parameter of ar model
  vector[maxt] states; // latent states
  vector[N] pred; // latent states
  // initial states
  ar_logit[1] = ar0;
  ar[1] = (-1 + 2*(exp(ar_logit[1])/(1+exp(ar_logit[1]))));
  states[1] = x0;
  // random walk in ar parameter
  for(i in 2:maxt) {
    // need to constrain the random walk, otherwise inf values are predicted
    //ar[i] = fmax(fmin(ar[i-1] + ar_dev[i-1], 1), -1);
    ar_logit[i] = ar_logit[i-1] + ar_dev[i-1];
    ar[i] = -1 + 2*exp(ar_logit[i])/(1+exp(ar_logit[i]));
    states[i] = drift + ar[i] * (states[i-1] - drift);
  }
  for(i in 1:N) {
    pred[i] = states[x[i]];
  }
}
model {
  ar0 ~ student_t(3, 0, 2);
  x0 ~ student_t(3, 0, 2);
  ar_sd ~ normal(0, ar_scale);
  ar_dev ~ normal(0, ar_sd); // process deviations
  drift ~ normal(0, 1);

  if(obs_model == 1) {
  CV ~ exponential(1);
  y ~ normal(pred, CV*pred);
  }
}
