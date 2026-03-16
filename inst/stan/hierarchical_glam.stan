// inst/stan/hierarchical_glam.stan

#include glam_likelihood.stan

data {
  int<lower=1> N;                
  int<lower=1> J;                
  int<lower=2> K;                
  array[N] int<lower=1,upper=J> subject; 
  array[N] int<lower=1,upper=K> choice;  
  vector[N] rt;                  
  matrix[N, K] v;                
  matrix[N, K] g;                
  real s_fixed;                  
  real<lower=0> min_rt;          
}

parameters {
  real mu_v_raw;                 
  real mu_gamma_raw;             
  real mu_sigma_raw;             
  real mu_s_raw;                 
  real mu_tau_raw;               

  // Constraints to prevent "Funnel" geometry issues
  vector<lower=0, upper=1>[5] tau_hier; 
  vector<lower=-3, upper=3>[J] v_offset;
  vector<lower=-3, upper=3>[J] gamma_offset;
  vector<lower=-3, upper=3>[J] sigma_offset;
  vector<lower=-3, upper=3>[J] s_offset;
  vector<lower=-3, upper=3>[J] tau_offset;
}

transformed parameters {
  vector<lower=0>[J] v_sub;      
  vector<lower=0,upper=1>[J] gamma_sub; 
  vector<lower=0>[J] sigma_sub;  
  vector<lower=0, upper=10>[J] s_sub;  // HARD CAP on sensitivity for stability
  vector<lower=0>[J] ndt_sub;    

  for (j in 1:J) {
    v_sub[j] = exp(mu_v_raw + tau_hier[1] * v_offset[j]);
    gamma_sub[j] = inv_logit(mu_gamma_raw + tau_hier[2] * gamma_offset[j]);
    sigma_sub[j] = exp(mu_sigma_raw + tau_hier[3] * sigma_offset[j]);
    
    // Sensitivity: cap at 10 to prevent logistic "cliffs"
    real s_val = exp(mu_s_raw + tau_hier[4] * s_offset[j]);
    s_sub[j] = (s_fixed > 0) ? s_fixed : fmin(s_val, 10.0);
    
    // NDT: Centered more naturally
    real raw_ndt = exp(mu_tau_raw + tau_hier[5] * tau_offset[j]);
    ndt_sub[j] = fmin(raw_ndt, min_rt * 0.8); // Slightly more room (80% of min_rt)
  }
}

model {
  // --- 1. Regularizing Priors ---
  // These guide the model toward plausible ranges
  mu_v_raw ~ normal(0, 0.5);
  mu_gamma_raw ~ normal(-0.5, 0.5); 
  
  // Sigma: Prior for roughly 0.1 - 0.3 range
  mu_sigma_raw ~ normal(-1.5, 0.5); 
  
  // Sensitivity: Prior for roughly 1.0 - 5.0 range
  mu_s_raw ~ normal(1.0, 0.5); 
  
  // NDT: Prior for roughly 200ms - 400ms range
  mu_tau_raw ~ normal(-1.2, 0.3); 
  
  tau_hier ~ normal(0, 0.5);
  v_offset ~ std_normal();
  gamma_offset ~ std_normal();
  sigma_offset ~ std_normal();
  s_offset ~ std_normal();
  tau_offset ~ std_normal();

  // --- 2. Likelihood ---
  for (n in 1:N) {
    int sj = subject[n];
    vector[K] drifts;
    for (k in 1:K) {
      real evidence_this = v[n, k] * (g[n, k] + gamma_sub[sj] * (1.0 - g[n, k]));
      real max_others = -1e9;
      for (i in 1:K) {
        if (i != k) {
          real evidence_other = v[n, i] * (g[n, i] + gamma_sub[sj] * (1.0 - g[n, i]));
          if (evidence_other > max_others) max_others = evidence_other;
        }
      }
      drifts[k] = calculate_glam_drift(evidence_this - max_others, v_sub[sj], s_sub[sj]);
    }
    target += glam_race_lpdf(rt[n] | choice[n], drifts, 1.0, ndt_sub[sj]);
  }
}
