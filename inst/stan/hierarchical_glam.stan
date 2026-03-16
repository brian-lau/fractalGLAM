// inst/stan/hierarchical_glam.stan

// Include the shared likelihood function
#include glam_likelihood.stan

data {
  int<lower=1> N;                // Total number of trials across all subjects
  int<lower=1> J;                // Number of subjects
  array[N] int<lower=1,upper=J> subject; // Subject index for each trial
  array[N] int<lower=0,upper=1> choice;  // Binary choice (1 = Right, 0 = Left)
  vector[N] rt;                  // Response time in milliseconds (or seconds)
  vector[N] v_left;              // Value of the left item (Z-scored)
  vector[N] v_right;             // Value of the right item (Z-scored)
  vector[N] g_left;              // Relative gaze to the left item (0-1)
  vector[N] g_right;             // Relative gaze to the right item (0-1)
  real s_fixed;                  // Fixed scaling constant (if > 0); else estimated
}

parameters {
  // Hyper-parameters (Group-level means on the unconstrained scale)
  real mu_v_raw;                 // Mean of log(velocity)
  real mu_gamma_raw;             // Mean of logit(gamma)
  real mu_sigma_raw;             // Mean of log(noise)
  real mu_s_raw;                 // Mean of log(scaling) if estimated

  // Hierarchical scales (Width of individual differences)
  vector<lower=0>[4] tau;

  // Individual-level offsets for Non-Centered Parameterization (NCP)
  // This helps the sampler navigate the "funnel" of hierarchical models.
  vector[J] v_offset;
  vector[J] gamma_offset;
  vector[J] sigma_offset;
  vector[J] s_offset;
}

transformed parameters {
  vector<lower=0>[J] v;          // Subject-level velocity
  vector<lower=0,upper=1>[J] gamma; // Subject-level attentional discount
  vector<lower=0>[J] sigma;      // Subject-level log-normal noise
  vector<lower=0>[J] s;          // Subject-level scaling constant

  for (j in 1:J) {
    // Non-Centered Parameterization: Param = exp(Group_Mean + Scale * Offset)
    v[j] = exp(mu_v_raw + tau[1] * v_offset[j]);
    gamma[j] = inv_logit(mu_gamma_raw + tau[2] * gamma_offset[j]);
    sigma[j] = exp(mu_sigma_raw + tau[3] * sigma_offset[j]);
    
    // Use fixed s if provided in data; otherwise estimate hierarchically
    s[j] = (s_fixed > 0) ? s_fixed : exp(mu_s_raw + tau[4] * s_offset[j]);
  }
}

model {
  // --- 1. Priors ---
  // Informative priors based on typical GLAM parameter ranges
  mu_v_raw ~ normal(1, 0.5);     // Velocity usually centers around exp(1) ≈ 2.7
  mu_gamma_raw ~ normal(-0.8, 0.5); // Gamma centers around logit(-0.8) ≈ 0.3
  mu_sigma_raw ~ normal(0, 0.5);
  mu_s_raw ~ normal(-8, 0.5);    // Scaling constant is often very small (e.g., 3e-4)
  
  tau ~ normal(0, 1);            // Standard half-normal for hierarchical scales
  
  // NCP offsets must follow a standard normal distribution
  v_offset ~ std_normal();
  gamma_offset ~ std_normal();
  sigma_offset ~ std_normal();
  s_offset ~ std_normal();

  // --- 2. Likelihood ---
  for (n in 1:N) {
    int sj = subject[n];
    
    // Call the SHARED function for RT drift
    // Note: We use v_right as 'this' and v_left as 'other' to define 
    // the signal relative to the Right choice.
    real drift_raw = calculate_glam_drift(v_right[n], v_left[n], g_right[n], g_left[n], 
                                          gamma[sj], v[sj], s[sj]);
    
    real drift = fmax(drift_raw, 1e-7); 

    // We still need the raw signal for the Choice bernoulli_logit
    real signal = (v_right[n] * (g_right[n] + gamma[sj] * g_left[n])) - 
                  (v_left[n] * (g_left[n] + gamma[sj] * g_right[n]));

    // Likelihoods
    rt[n] ~ lognormal(log(1.0 / drift), sigma[sj]); 
    choice[n] ~ bernoulli_logit(signal * 0.5);
  }
}
