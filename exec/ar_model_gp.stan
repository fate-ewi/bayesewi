data {
  int<lower=1> N; // number of data points
  vector[N] y; // observed data
  int<lower=1> x[N]; // time coordinates of observed data
  int<lower=1> maxt; // time steps
  matrix[maxt, maxt] distmat2; // matrix of squared distances of *unique* time steps
  vector[maxt] deltat; // difference in time for uniquet
  int<lower=1> obs_model; // 1 = normal
}
parameters {
  real<lower=0> sigma_obs; // observation variance sigma
  vector[maxt-1] logit_ar; // logit_ar
  vector[maxt-1] pro_dev; // deviations
  real x0; // initial state
  // real ar0; // initial state
  // real drift;
  real ar_mu; // mean ar
  real<lower=0> gp_sigma_sq;
  real<lower=0> gp_scale;
  real<lower=0> sigma_pro;
}
transformed parameters {
  vector[maxt-1] zeros; // mean for cov function
  vector[N] pred; // predictions
  vector[maxt] scalet;
  cov_matrix[maxt-1] Sigma; // cov matrix for gp of log(sigma2)
  vector[maxt] states;
  vector[maxt] ar;

  // use gp model to estimate variances (log sigma2)
  for(i in 1:(maxt-1)) {
    zeros[i] = 0;
    for(j in 1:(maxt-1)) {
      Sigma[i,j] = gp_sigma_sq * exp(-gp_scale * distmat2[i,j]);
    }
  }

  // evolution is Fokker–Planck Ornstein–Uhlenbeck process, ie https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process
  states[1] = x0;
  ar[1] = 0;
  scalet[1] = 1;
  for(i in 2:maxt) {
    ar[i] = exp(logit_ar[i-1])/(1 + exp(logit_ar[i-1]));
    // 0 is placeholder for drift

    scalet[i] = sqrt( (1 - exp(-2*ar[i]*deltat[i]))/ (2*ar[i]));

    states[i] = 0 + exp(-ar[i] * deltat[i]) * (states[i-1] - 0) + pro_dev[i-1] * (scalet[i] * sigma_pro);

  }

  for(i in 1:N) {
    pred[i] = states[x[i]];
  }
}
model {
  // drift ~ normal(0,1);
  x0 ~ student_t(3, 0, 2);
  // ar0 ~ student_t(3, 0, 2);
  gp_scale ~ student_t(3, 0, 2);
  gp_sigma_sq ~ student_t(3, 0, 2);
  ar_mu ~ normal(0, 3);
  logit_ar ~ multi_normal(zeros+ar_mu, Sigma);
  sigma_pro ~ student_t(3, 0, 2);
  pro_dev ~ normal(0, 1);

  if(obs_model == 1) {
  sigma_obs ~ student_t(3, 0, 2);
  y ~ normal(pred, sigma_obs);
  }
}
