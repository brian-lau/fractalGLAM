// inst/stan/include/glam_likelihood.stan

functions {
  /**
  * Log-Survival Function for Single-Boundary Wiener Process
  * Calculates log(P(T > t)), the probability an accumulator has NOT hit 
  * the boundary by time t.
  */
  real wiener_race_survival_lccdf(real t, real a, real delta, real sigma) {
    // Safety check: if time is effectively zero, survival is 100% (log(1)=0)
    if (t < 1e-7) return 0.0; 
    
    // Using fmax(t, 1e-7) inside the calculation for numerical stability
    return std_normal_lcdf((a - delta * t) / (sigma * sqrt(fmax(t, 1e-7))));
  }
  
  /**
  * Inverse Gaussian Log-Density
  * Correct likelihood for a Wiener process with a single absorbing boundary.
  */
  real inverse_gaussian_lpdf(real x, real mu, real lambda) {
    if (x <= 0) return -1e10;
    return 0.5 * log(lambda / (2 * pi() * pow(x, 3))) - 
    (lambda * square(x - mu)) / (2 * square(mu) * x);
  }
  
  /**
  * Linear Drift (Stable for 2-item sets)
  */
  real calculate_glam_drift_simple(real abs_ev, real nu) {
    return nu * abs_ev;
  }
  
  /**
  * Scaled Logistic Drift (Robust for N-item sets)
  */
  real calculate_glam_drift_multinomial(real rel_signal, real abs_ev, real nu, real s) {
    return nu * abs_ev / (1.0 + exp(-s * rel_signal));
  }
  
  /**
  * GLAM Race Likelihood (N-Item Inverse Gaussian)
  * This matches the single-boundary simulation logic in R.
  */
  real glam_race_lpdf(real rt, int winner_idx, vector drifts, real boundary, real tau, real sigma) {
    int K_items = num_elements(drifts);
    real lp = 0;
    real t_adj = rt - tau;
    
    // Safety guard for Non-Decision Time
    if (t_adj <= 0) return -1e10;
    
    for (k in 1:K_items) {
      real d = fmax(drifts[k], 1e-4); // Ensure positive drift
      
      if (k == winner_idx) {
        // 1. Density of the Winner (Inverse Gaussian)
        // mu = a/d, lambda = a^2/sigma^2
        real mu = boundary / d;
        real lambda = square(boundary) / square(sigma);
        lp += inverse_gaussian_lpdf(t_adj | mu, lambda);
      } else {
        // 2. Log-Survival of the Losers
        // Vertical bar (|) is required here for the _lccdf suffix
        lp += wiener_race_survival_lccdf(t_adj | boundary, d, sigma);
      }
    }
    return lp;
  }
}
