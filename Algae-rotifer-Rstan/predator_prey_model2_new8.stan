// Fit to Cv and Rotifer data
//Model 2  conversion rate may change under salinity
functions {
  real[] LV(real t,
            real[] y,
            real[] theta,
            real[] x_r,
            int[] x_i) {
    real r = x_r[1];
    real K = x_r[2];
   // real c=x_r[3];
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
  real log_m;            // baseline log m
  real log_a;
  real c_raw;

  real<lower=0> sigma[2];
}


transformed parameters {
  real a[N];
  real m[N];
  real c[N];
  for (n in 1:N) {
    m[n] = exp(log_m);
    a[n] = exp(log_a);
    c[n] = c_raw;
    } 
  }

model {
  log_m ~ normal(log(0.1), 1); // mortality rate wihtout salt would be low
  log_a ~ normal(log(1), 0.5); // rough guess
  c_raw ~ normal(1.5, 0.2); //  rough approximation from Becks et al (2010)
  sigma ~ normal(0, 0.5); // half normal since sigma is positive

  for (n in 1:N) {
    real theta[3] = {a[n], m[n], c[n]}; #
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
  vector[N*T*2] log_lik;//use for WAIC

  {
    int k = 1;

    for (n in 1:N) {

      real theta[3] = {a[n], m[n], c[n]};

      y_pred[n] =
        integrate_ode_rk45(
          LV,
          y0[n],
          0.0,
          ts,
          theta,
          x_r[n],
          x_i
        );

      for (t in 1:T) {

        // prey
        log_lik[k] =
          lognormal_lpdf(
            y[n,t,1] |
            log(fmax(y_pred[n,t,1], 1e-3)),
            sigma[1]
          );
        k += 1;

        // predator
        log_lik[k] =
          lognormal_lpdf(
            y[n,t,2] |
            log(fmax(y_pred[n,t,2], 1e-3)),
            sigma[2]
          );
        k += 1;

        // posterior predictive samples
        z_pred[n,t,1] =
          lognormal_rng(
            log(fmax(y_pred[n,t,1], 1e-3)),
            sigma[1]
          );

        z_pred[n,t,2] =
          lognormal_rng(
            log(fmax(y_pred[n,t,2], 1e-3)),
            sigma[2]
          );
      }
    }
  }
}
