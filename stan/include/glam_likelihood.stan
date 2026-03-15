/**
 * Core GLAM Likelihood Function
 * Calculates the decision signal and resulting drift rate.
 */
real calculate_glam_drift(real v_this, real v_other, real g_this, real g_other, real gamma, real nu, real s) {
    // 1. Calculate the relative signal (same as before)
    real signal = (v_this * (g_this + gamma * g_other)) - (v_other * (g_other + gamma * g_this));
    
    // 2. Magnitude-based drift (The Fix for the Hump)
    // We use the absolute signal to ensure RT is slow at 0 and fast at +/- extremes.
    // The + 0.01 floor prevents division by zero.
    return nu * s * (fabs(signal) + 0.01); 
}
