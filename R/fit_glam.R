#' Fit GLAM via Hierarchical Bayesian Inference
#' @param data Data frame from prep_glam_data.
#' @param method "meanfield" (VB) or "sample" (MCMC).
#' @param s_val Numeric value for fixed s, or NULL to estimate.
#' @param ... Additional arguments passed to the cmdstanr fit methods (e.g., iter, chains).
#' @import cmdstanr
#' @export
fit_glam <- function(data, method = "meanfield", s_val = 3e-4, ...) {
  
  # 1. Robust Path Check: Only attempt to set if path isn't already valid
  # This avoids the dir.exists() 'invalid filename' error
  has_path <- tryCatch({
    # If this works, the path is already set correctly
    dir.exists(cmdstanr::cmdstan_path())
  }, error = function(e) FALSE)
  
  if (!has_path) {
    # Only try to set if it's not already there. 
    # If this fails, we want it to fail gracefully with a clear message.
    tryCatch({
      cmdstanr::set_cmdstan_path()
    }, error = function(e) {
      stop("CmdStan path is not set and the default location could not be verified. 
            Please run cmdstanr::set_cmdstan_path('/your/path') in your console.")
    })
  }
  
  # 2. Proceed with model loading
  model_path <- system.file("stan", "hierarchical_glam.stan", package = "fractalGLAM")
  include_dir <- system.file("stan", "include", package = "fractalGLAM")
  
  # mod <- cmdstanr::cmdstan_model(model_path, include_paths = include_dir)
  
  # Robustly load the model
  mod <- cmdstanr::cmdstan_model(
    model_path,
    include_paths = include_dir, # Ensure this is the folder path
    force_recompile = TRUE # Force recompilation to ensure the model is up-to-date
    )
  
  stan_list <- list(
    N = nrow(data),
    J = length(unique(data$subject_id)),
    subject = as.numeric(as.factor(data$subject_id)),
    choice = data$choice,
    rt = data$rt,
    v_left = data$value_left,
    v_right = data$value_right,
    g_left = data$gaze_left,
    g_right = data$gaze_right,
    s_fixed = ifelse(is.null(s_val), -1.0, s_val)
  )
  
  if (method == "meanfield") {
    return(mod$variational(data = stan_list, ...))
  } else {
    return(mod$sample(data = stan_list, ...))
  }
}
