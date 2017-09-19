data {
  int<lower=1> N; // number of data points
  vector[N] y; // observed data
  int<lower=1> x[N]; // time coordinates of observed data
  int<lower=1> maxt; // time steps
  matrix[maxt, maxt] distmat2; // matrix of squared distances of *unique* time steps
  vector[maxt] deltat; // difference in time for uniquet
  int<lower=1> obs_model; // 1 = normal, 2 gamma, etc
}
parameters {
  real<lower=0> sigma_obs; // observation variance sigma
  vector[maxt-1] log_sigma2; // log_sigma2
  vector[maxt-1] pro_dev;
  real x0; // initial state
  //real drift;
  real sigma2_mu; // mean sigma2
  real<lower=0> gp_sigma_sq;
  real<lower=0> gp_scale;
}
transformed parameters {
  vector[maxt-1] zeros; // mean for cov function
  vector[N] pred; // predictions
  cov_matrix[maxt-1] Sigma; // cov matrix for gp of log(sigma2)
  vector[maxt] states;
  vector[maxt] sigma;

  // use gp model to estimate variances (log sigma2)
  for(i in 1:(maxt-1)) {
    zeros[i] = 0;
    for(j in 1:(maxt-1)) {
      Sigma[i,j] = gp_sigma_sq * exp(-gp_scale * distmat2[i,j]);
    }
  }

  // evolution is brownian motion, ie https://en.wikipedia.org/wiki/Geometric_Brownian_motion
  states[1] = x0;
  sigma[1] = 0;
  for(i in 2:maxt) {
    sigma[i] = sqrt(exp(log_sigma2[i-1]));
    // 0 is placeholder for drift term
    states[i] = states[i-1] * exp( (0 - sigma[i]*sigma[i]/2) * deltat[i] + sigma[i]* pro_dev[i-1]);
  }

  for(i in 1:N) {
    pred[i] = states[x[i]];
  }
}
model {
  pro_dev ~ normal(0, 1); // process deviations
  //drift ~ normal(0,1);
  x0 ~ student_t(3, 0, 2);
  gp_scale ~ student_t(3, 0, 2);
  gp_sigma_sq ~ student_t(3, 0, 2);
  sigma2_mu ~ normal(0, 3);
  log_sigma2 ~ multi_normal(zeros+sigma2_mu, Sigma);

  if(obs_model == 1) {
  sigma_obs ~ student_t(3, 0, 2);
  y ~ normal(pred, sigma_obs);
  }
}
