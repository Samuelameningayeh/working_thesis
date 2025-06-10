functions {
  vector seir(real t,
              vector y,
              real beta,
              real sigma,
              real gamma,
              real N) {
    vector[5] dydt;
    real S = y[1];
    real E = y[2];
    real I = y[3];
    real R = y[4];

    dydt[1] = -beta * I * S / N;              // dS/dt
    dydt[2] = beta * I * S / N - sigma * E;   // dE/dt
    dydt[3] = sigma * E - gamma * I;          // dI/dt
    dydt[4] = gamma * I;                      // dR/dt
    dydt[5] = sigma * E;                      // Cumulative incidence

    return dydt;
  }
}
data {
  int<lower=1> N_countries;                // Number of countries (3)
  int<lower=1> n_days;                     // Number of time points
  array[N_countries] vector[5] y0;         // Initial conditions for each country
  real t0;                                 // Initial time
  array[n_days] real t;                    // Time points
  array[N_countries] int N;                // Population sizes
  array[N_countries, n_days] int<lower=0> cases; // Observed cases for each country
  array[N_countries] real beta_value;
  real gamma_value;
  real sigma_value;
  real<lower=0,upper=1> reporting_rate;
}
parameters {
  vector<lower=0>[N_countries] beta;       // Country-specific transmission rates
  real<lower=0> sigma;                     // Shared recovery rate
  real<lower=0> gamma;                     // Shared progression rate (E to I)
  real<lower=0> phi_inv;                   // Negative binomial overdispersion
}
transformed parameters {    
  array[N_countries, n_days] vector[5] y;   // SEIR states for each country
  array[N_countries, n_days] real incidence; // Incidence for each country
  real<lower=0> phi = 1.0 / phi_inv;       // Negative binomial dispersion
  array[N_countries, n_days] real adj_incidence;
  
  //reporting_rate = 0.8;

  for (c in 1:N_countries) {
    // Solve ODE for each country
    y[c] = ode_rk45(seir, y0[c], t0, t, beta[c], sigma, gamma, N[c] + 0.0);
    
    // Compute incidence
    incidence[c, 1] = y[c, 1, 5];         // Initial incidence
    adj_incidence[c, 1] = reporting_rate*incidence[c, 1]+0.00001;
    for (i in 2:n_days) {
      incidence[c, i] = y[c, i, 5] - y[c, i-1, 5];
      adj_incidence[c, i] = reporting_rate*incidence[c, i]+0.00001;
    }
  }
}
model {
  // Priors
  for (i in 1:N_countries){
    beta[i] ~ lognormal(log(beta_value[i]), 0.5);        // Infection rate prior
  }
  sigma ~ lognormal(log(sigma_value), 0.5);      // Prior mean: 4-day recovery period
  gamma ~ lognormal(log(gamma_value), 0.5);      // Prior mean: 6-day incubation period
  phi_inv ~ exponential(5);
  //reporting_rate ~ beta(8, 2);             // Centered around 0.8

  // Likelihood
  for (c in 1:N_countries) {
    for (i in 1:N_countries){
      cases[c, i] ~ neg_binomial_2(adj_incidence[c,i], phi);
  }
  }
}
generated quantities {
  vector[N_countries] R0 = beta / gamma;   // Country-specific R0
  real recovery_time = 1.0 / gamma;        // Shared recovery time
  real incubation_period = 1.0 / sigma;    // Shared incubation period
  array[N_countries, n_days] real pred_incidence; // Predicted cases

  for (c in 1:N_countries) {
    for (i in 1:N_countries){
      pred_incidence[c, i] = neg_binomial_2_rng(adj_incidence[c, i], phi);
  }
  }
}