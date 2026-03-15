functions {
  real glam_drift(real v_this, real v_other, real g_this, real g_other, real gamma, real nu, real s) {
    real signal = (v_this * (g_this + gamma * g_other)) - (v_other * (g_other + gamma * g_this));
    return nu * s * (1 / (1 + exp(-signal))); 
  }
}

data {
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1,upper=J> subject; 
  array[N] int<lower=0,upper=1> choice; 
  vector[N] rt;
  vector[N] v_left;
  vector[N] v_right;
  vector[N] g_left;
  vector[N] g_right;
  real s_fixed; // If > 0, use this; if -1, estimate.
}

parameters {
  real mu_v_raw;
  real mu_gamma_raw;
  real mu_sigma_raw;
  real mu_s_raw;
  vector<lower=0>[4] tau;
  vector[J] v_offset;
  vector[J] gamma_offset;
  vector[J] sigma_offset;
  vector[J] s_offset;
}

transformed parameters {
  vector<lower=0>[J] v;
  vector<lower=0,upper=1>[J] gamma;
  vector<lower=0>[J] sigma;
  vector<lower=0>[J] s;

  for (j in 1:J) {
    v[j] = exp(mu_v_raw + tau[1] * v_offset[j]);
    gamma[j] = inv_logit(mu_gamma_raw + tau[2] * gamma_offset[j]);
    sigma[j] = exp(mu_sigma_raw + tau[3] * sigma_offset[j]);
    s[j] = (s_fixed > 0) ? s_fixed : exp(mu_s_raw + tau[4] * s_offset[j]);
  }
}

model {
  // Priors: Adjusted for better recovery
  mu_v_raw ~ normal(1, 0.5);       // Pull mu_v slightly higher
  mu_gamma_raw ~ normal(-0.8, 0.5); // Prior centered around 0.3 in logit space
  mu_sigma_raw ~ normal(0, 0.5);
  mu_s_raw ~ normal(-8, 0.5);
  
  tau ~ normal(0, 1);              // Tighter SD priors for group-level offsets
  
  v_offset ~ std_normal();
  gamma_offset ~ std_normal();
  sigma_offset ~ std_normal();
  s_offset ~ std_normal();
  
  for (n in 1:N) {
    int sj = subject[n];
    
    // Calculate the raw signal for Choice
    real signal = (v_right[n] * (g_right[n] + gamma[sj] * g_left[n])) - 
    (v_left[n] * (g_left[n] + gamma[sj] * g_right[n]));
    
    // // Calculate magnitude-based drift for RT
    // real drift = v[sj] * s[sj] * (abs(signal) + 0.01);
    // 
    // // Ensure drift is never exactly zero
    // if (drift <= 0) drift = 1e-9;
    // 1. Calculate raw drift
    real drift_raw = v[sj] * s[sj] * (abs(signal) + 0.01);
    
    // 2. FORCE a minimum drift for numerical stability
    // 1e-7 ensures log(1/drift) never exceeds ~16
    real drift = fmax(drift_raw, 1e-7);

    // Likelihoods
    rt[n] ~ lognormal(log(1.0 / drift), sigma[sj]); 
    choice[n] ~ bernoulli_logit(signal * 0.5); // Choice follows the sign of the signal
  }
}
