functions {
  real[] LV(real t,
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {
    real r = x_r[1];
    real K = x_r[2];
    real m = x_r[3];
    real a = theta[1];
    real c = theta[2];
    real prey = y[1];
    real pred = y[2];
    real dydt[2];
    dydt[1] = (r - r*prey / K - a * pred) * prey;
    dydt[2] = (c * a * prey - m ) * pred;
     return dydt;
  }
}

data {
  int<lower=1> T;
  //real<lower=0> prey0;
  //real<lower=0> pred0;
  real <lower=0> y0[2];
  real y[T, 2];
  real t0;
  real ts[T];
  real r;
  real K;
  real m; 
  real a_mean;
  real a_sd;
  real c_mean;
  real c_sd;
  real sigma1_mean;
  real sigma1_sd;
  real sigma2_mean;
  real sigma2_sd;

}

transformed data {
  real x_r[3];  // r and K are fixed parameters
  x_r[1]=r; // intrinsic growrth estimated from monoculture
  x_r[2]=K; // carrying capacity estimated by monoculture
  x_r[3]=m; // conversion rate from the literature
  int x_i[0];   
}

parameters {
  real<lower=0> theta[2];
  real<lower=0> sigma[2];
}

model {
  theta[1] ~ normal(a_mean, a_sd);
  theta[2] ~ noraml(c_mean, c_sd)
  sigma[1] ~ normal(sigma1_mean, sigma1_sd);
  sigma[2] ~ normal(sigma2_mean, sigma2_sd);
  // ode
  real y_hat[T, 2];
  y_hat = integrate_ode_rk45(LV, y0, t0, ts, theta, x_r, x_i);
  // likelihood
  for (t in 1:T) {
      y[t,1] ~ normal(y_hat[t, 1], sigma[1]);
      y[t,2] ~normal(y_hat[t, 2], sigma[2]);
  }
}

generated quantities{
  real y_pred[T,2];
  real z_pred[T,2];
  y_pred = integrate_ode_rk45(LV, y0, t0, ts, theta, x_r, x_i );
  for (t in 1:T) {
    z_pred[t,1] = y_pred[t,1] + normal_rng(0,sigma[1]);
    z_pred[t,2] = y_pred[t,2] + normal_rng(0,sigma[2]);
  }
}

