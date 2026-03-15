// inst/stan/restricted_glam.stan
#include include/glam_likelihood.stan

data {
  int<lower=1> N;                // Total trials
  int<lower=1> J;                // Number of subjects
  array[N] int<lower=1,upper=J> subject; 
  array[N] int<lower=0,upper=1> choice;  
  vector[N] rt;
  vector[N] v_left;
  vector[N] v_right;
  vector[N] g_left;
  vector[N] g_right;
  real s_fixed;                  
}

parameters {
  // Gamma is excluded as it is restricted to 1.0
  real mu_v_raw;
  real mu_sigma_raw;
  real mu_s_raw;

  vector<lower=0>[3] tau;

  vector[J] v_offset;
  vector[J] sigma_offset;
  vector[J] s_offset;
}

transformed parameters {
  vector<lower=0>[J] v;
  vector<lower=0>[J] sigma;
  vector<lower=0>[J] s;
  real gamma_restricted = 1.0; // The restriction

  for (j in 1:J) {
    v[j] = exp(mu_v_raw + tau[1] * v_offset[j]);
    sigma[j] = exp(mu_sigma_raw + tau[2] * sigma_offset[j]);
    
    if (s_fixed > 0) {
      s[j] = s_fixed;
    } else {
      s[j] = exp(mu_s_raw + tau[3] * s_offset[j]);
    }
  }
}

model {
  // Priors
  mu_v_raw ~ normal(0, 1);
  mu_sigma_raw ~ normal(0, 1);
  mu_s_raw ~ normal(-8, 1);
  tau ~ cauchy(0, 2.5);
  
  v_offset ~ std_normal();
  sigma_offset ~ std_normal();
  s_offset ~ std_normal();

  // Likelihood
  for (n in 1:N) {
    int sj = subject[n];
    // Call the shared function with gamma = 1.0
    real drift = calculate_glam_drift(v_right[n], v_left[n], g_right[n], g_left[n], 
                                     gamma_restricted, v[sj], s[sj]);
    
    rt[n] ~ normal(1.0 / drift, sigma[sj]); 
    choice[n] ~ bernoulli_logit(drift); 
  }
}