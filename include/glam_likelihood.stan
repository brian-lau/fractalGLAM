/**
 * Core GLAM Likelihood Function
 * Calculates the decision signal and resulting drift rate.
 */
/*real calculate_glam_drift(real v_this, real v_other, real g_this, real g_other, real gamma, real nu, real s) {
    // The Gaze-weighted decision signal (r)
    // r = nu * [ (v_this * (g_this + gamma * g_other)) - (v_other * (g_other + gamma * g_this)) ]
    real signal = (v_this * (g_this + gamma * g_other)) - (v_other * (g_other + gamma * g_this));
    
    // Scale by velocity (nu) and scaling constant (s)
    // We use a sigmoid transform to ensure the drift rate is mapped appropriately
    // Use 1e-7 to prevent the drift from ever being exactly zero
    return nu * s * (1.0 / (1.0 + exp(-signal))) + 1e-7;
}*/

real calculate_glam_drift(real v_this, real v_other, real g_this, real g_other, real gamma, real nu, real s) {
    // 1. Calculate the relative signal (same as before)
    real signal = (v_this * (g_this + gamma * g_other)) - (v_other * (g_other + gamma * g_this));
    
    // 2. Magnitude-based drift (The Fix for the Hump)
    // We use the absolute signal to ensure RT is slow at 0 and fast at +/- extremes.
    // The + 0.01 floor prevents division by zero.
    return nu * s * (fabs(signal) + 0.01); 
}
