functions {
  real[] LV(real t,
            real[] y,
            real[] theta,
            real[] x_r,
            int[] x_i) {
    real r = x_r[1];
    real K = x_r[2];
    real a = theta[1];
    real m = theta[2];
    real c = theta[3];
    real prey = y[1];
    real pred = y[2];
    real dydt[2];
    dydt[1] = r * prey * (1 -  prey / K) - a * prey * pred;
    dydt[2] = (c * a * prey - m) * pred;
    return dydt;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=2> salt[N];
  int<lower=1, upper=2> prey_type[N];
  real y0[N, 2];
  real y[N, T, 2];
  real ts[T];
  real r[2];
  real K[2];
  real m_mu[2];
  real m_sigma[2];
}

transformed data {
  int x_i[0];
  real x_r[N, 2];
  for (n in 1:N) {
    x_r[n, 1] = r[salt[n]];
    x_r[n, 2] = K[salt[n]];
  }
}

parameters {
  real <lower=0> a_row[2];
  real<lower=0> c_raw;
  real log_m[2];  
  real<lower=0> sigma[2];
}


transformed parameters {
  real a[N];
  real m[N];
  real c[N];
  for (n in 1:N) {
     m[n] = exp(log_m[salt[n]]);
    a[n] = ((salt[n] == 2 ? a_row[2] : a_row[1]));// corresponding to the diff in scaling parms
    c[n] = c_raw;
  }
}

model {
  a_row ~ normal(0.8, 0.25); // attack rate would be  smaller than Cv. The mean is based on initial guess 
  c_raw ~ normal(1.7, 0.2); //  rough approximation from Becks et al (2010)
  log_m[1] ~ normal(m_mu[1], m_sigma[1]);
  log_m[2] ~ normal(m_mu[2], m_sigma[2]);
  sigma ~ normal(0, 0.5); // half normal since sigma is positive

  for (n in 1:N) {
    real theta[3] = {a[n], m[n],  c[n]};
    real y_hat[T, 2] = integrate_ode_rk45(LV, y0[n], 0.0, ts, theta, x_r[n], x_i);
    for (t in 1:T) {
      y[n, t, 1] ~ lognormal(log(fmax(y_hat[t, 1],1e-3)), sigma[1]); // avoid negative scale value
      y[n, t, 2] ~ lognormal(log(fmax(y_hat[t, 2],1e-3)), sigma[2]);
    }
  }
}

generated quantities {
  real y_pred[N, T, 2];
  real z_pred[N, T, 2];
  
  for (n in 1:N) {
     real theta[3] = {a[n],  m[n], c[n]};
    y_pred[n] = integrate_ode_rk45(LV, y0[n], 0.0, ts, theta, x_r[n], x_i);
    for (t in 1:T) {
      z_pred[n, t, 1] = lognormal_rng(log(fmax(y_pred[n, t, 1], 1e-3)), sigma[1]);
      z_pred[n, t, 2] = lognormal_rng(log(fmax(y_pred[n, t, 2], 1e-3)), sigma[2]);
    }
  }
}
