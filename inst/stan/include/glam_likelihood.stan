// inst/stan/include/glam_likelihood.stan

functions {
  /**
  * Core GLAM Drift Calculation
  * @param v_this Value of the item being considered
  * @param v_other Value of the alternative item
  * @param g_this Gaze proportion for the item being considered
  * @param g_other Gaze proportion for the alternative item
  * @param gamma Attentional discount parameter (0-1)
  * @param nu Velocity/Gain parameter
  * @param s Scaling constant
  * @return The magnitude-based drift rate
  */
  real calculate_glam_drift(real v_this, real v_other, real g_this, real g_other, 
                            real gamma, real nu, real s) {
      
      // 1. Calculate the relative decision signal
      real signal = (v_this * (g_this + gamma * g_other)) - 
                    (v_other * (g_other + gamma * g_this));
      
      // 2. Return magnitude-based drift with stability floor
      return nu * s * (abs(signal) + 0.01); 
  }
}
