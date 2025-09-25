// Logistic function with shift

functions {
  real[] logisticgrowth(real t,
                        real[] y,
                        real[] theta,
                        real[] x_r,
                        int[] x_i) {
    int n = x_i[1];
    real dydt[n];
    for (i in 1:n) {
      dydt[i] = theta[2*i - 1] * y[i] * (1 - y[i] / theta[2*i]);
    }
    return dydt;
  }
}

data {
  int<lower=1> T;
  int<lower=1> n_samples;
  real<lower=0> y0[n_samples];
  real<lower=0> z[T, n_samples];
  real<lower=0> t0;
  real<lower=0> ts[T];
  int<lower=0, upper=1> NaCl[n_samples]; // 0: 無塩, 1: 有塩
}

transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1] = n_samples;
}

parameters {
  real<lower=0> mu_r;  
  real<lower=0> mu_K;  
  real delta_r; // shift in r due to salt  
  real delta_K; // shift in K due to salt
  real<lower=0> sigma;
}

transformed parameters {
  real theta[2 * n_samples];
  real r_cond[2];  
  real K_cond[2];  
  r_cond[1] = mu_r;           // NaCl = 0  
  r_cond[2] = mu_r + delta_r; // NaCl = 1  
  K_cond[1] = mu_K;  
  K_cond[2] = mu_K + delta_K;
  for (i in 1:n_samples) {
    real r_i = r_cond[NaCl[i] + 1]; // NaCl[i] = 0 or 1 → index 1 or 2
    real K_i = K_cond[NaCl[i] + 1];
    theta[2*i - 1] = r_i;
    theta[2*i]     = K_i;
  }
}

model {
  real y_hat[T, n_samples];

  r_cond ~ cauchy(0, 2.5);
  K_cond ~ normal(1.0, 0.2);
  delta_r ~ normal(0, 0.25);  // 
  delta_K ~ normal(0, 0.25); 
  sigma ~ normal(0, 0.5); // half normal since sigma is positive 

  y_hat = integrate_ode_rk45(logisticgrowth, y0, t0, ts, theta, x_r, x_i);

  for (t in 1:T) {
    for (i in 1:n_samples) {
      z[t, i] ~ lognormal(log(fmax(y_hat[t, i], 1e-3)), sigma);
    }
  }
}

generated quantities {
  real y_pred[T, n_samples];
  real z_pred[T, n_samples];

  y_pred = integrate_ode_rk45(logisticgrowth, y0, t0, ts, theta, x_r, x_i);

  for (t in 1:T) {
    for (i in 1:n_samples) {
      z_pred[t, i] = lognormal_rng(log(fmax(y_pred[t, i], 1e-3)), sigma);
    }
  }
}

