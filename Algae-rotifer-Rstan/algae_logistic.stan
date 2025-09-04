//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  real[] logisticgrowth(real t,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i
                  ) {
    real dydt[x_i[1]];
    for (i in 1:x_i[1]){
      dydt[i] = theta[1] * y[i] * (1-y[i]/theta[2]);
    }
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=1> n_samples;
  real y0[n_samples];
  real z[T,n_samples];
  real t0;
  real ts[T];
}
transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1] = n_samples;
}
parameters {
  real<lower=0> theta[2];
  real<lower=0> sigma;
}
model {
  real y_hat[T,n_samples];
  theta[1] ~ cauchy(0,2.5);
  theta[2] ~ normal(1e4,1e3);
  sigma ~ normal(0,0.05);
  y_hat = integrate_ode_rk45(logisticgrowth, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T) {
    for (i in 1:n_samples) {
      z[t,i] ~ normal(y_hat[t,i], sigma);
    }
  }
}
generated quantities{
  real y_pred[T,n_samples];
  real z_pred[T,n_samples];
  y_pred = integrate_ode_rk45(logisticgrowth, y0, t0, ts, theta, x_r, x_i );
  for (t in 1:T) {
    for(i in 1:n_samples){
      z_pred[t,i] = y_pred[t,i] + normal_rng(0,sigma);
    }
  }
}

