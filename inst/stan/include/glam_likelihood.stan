// inst/stan/include/glam_likelihood.stan

functions {
  /**
   * Wiener CDF Approximation
   * Calculates the probability of an accumulator being below boundary at time t.
   */
  real wiener_cdf_approx(real t, real alpha, real tau, real beta, real delta) {
    if (t <= tau) return 0.0;
    real t_adj = t - tau;
    // Large-time approximation for the Wiener CDF with single boundary
    return 1.0 - exp(-2.0 * delta * alpha * (1.0 - beta) - 0.5 * delta^2 * t_adj);
  }

  /**
   * Logistic Drift Rate Mapping (glambox spec)
   */
  real calculate_glam_drift(real signal, real nu, real s) {
    return nu / (1.0 + exp(-s * signal));
  }

  /**
   * GLAM Race Likelihood (N-Item)
   */
  real glam_race_lpdf(real rt, int winner_idx, vector drifts, real boundary, real tau) {
    int K_items = num_elements(drifts);
    real lp = 0;
    
    // Clamp drifts to a numerically stable range (prevents gradients from exploding)
    vector[K_items] safe_drifts;
    for(k in 1:K_items) safe_drifts[k] = fmin(fmax(drifts[k], 1e-4), 50.0);
    
    // 1. Likelihood of the Winner
    lp += wiener_lpdf(rt | boundary, tau, 0.5, safe_drifts[winner_idx]);
    
    // 2. Likelihood of the Losers (The Competition)
    for (j in 1:K_items) {
      if (j != winner_idx) {
        real cdf = wiener_cdf_approx(rt, boundary, tau, 0.5, safe_drifts[j]);
        lp += log(fmax(1.0 - cdf, 1e-10)); 
      }
    }
    return lp;
  }
}