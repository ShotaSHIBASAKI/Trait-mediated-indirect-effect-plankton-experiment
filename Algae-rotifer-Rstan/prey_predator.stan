
functions {
  real[] LV(real t,
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {
    real r = x_r[1];
    real K = x_r[2];
    real c = x_r[3];
    real a = theta[1];
    real m = theta[2];
    real prey = y[1];
    real pred = y[2];
    real dydt[2];
    dydt[1] = (r - r*prey / K - a * pred) * prey;
    dydt[2] = (c * a *  prey - m ) * pred;
     return dydt;
  }
}

// functions {
//   // non-dim ver
//   real[] LV(real t,
//              real[] y,
//              real[] theta,
//              real[] x_r,
//              int[] x_i) {
//     real r = x_r[1];
//     real K = x_r[2];
//     real a = theta[1];
//     real c = theta[2];
//     real m = theta[3];
//     real alpha = m/(r*c*K);
//     real beta = m/r;
//     real gamma = K*a*c/m;
//     real prey = y[1]/K; // non-dim
//     real pred = y[2]*(c*a*K*r)/m; // non dim
//     real dydt[2];
//     dydt[1] = (1 - prey -  alpha*pred) * prey;
//     dydt[2] = beta*(gamma * prey - 1) * pred;
//      return dydt;
//   }
// }

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
  real c;
}

transformed data {
  real x_r[3];  // r and K are fixed parameters
  x_r[1]=r; // intrinsic growrth estimated from monoculture
  x_r[2]=K;// carrying capacity estimated by monoculture
  x_r[3]=c;// set due to the multicoliniarity
  int x_i[0];   
}

parameters {
  real<lower=0> theta[2];
  real<lower=0> sigma[2];
}

model {
  // priors
  theta[1] ~ exponential(0.5); // attack rate
  theta[2] ~ exponential(0.5); // mortality
  sigma[1] ~ normal(0,0.05);
  sigma[2] ~ normal(0,0.05);
  //ode
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


