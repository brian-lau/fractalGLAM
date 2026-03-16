#' Fit GLAM via Hierarchical Bayesian Inference
#' 
#' This function interfaces with CmdStan to fit the Gaze-weighted Linear 
#' Accumulator Model. It supports both Variational Bayes (fast approximation) 
#' and MCMC (full posterior sampling).
#' 
#' @param data Data frame prepared via \code{prep_glam_data}.
#' @param method Character, either "meanfield" (Variational Bayes) or 
#'        "sample" (Hamiltonian Monte Carlo). Defaults to "meanfield".
#' @param s_val Numeric. A fixed scaling constant (e.g., 3e-4). If NULL, 
#'        the model will attempt to estimate 's' as a free parameter.
#' @param force_recompile Logical. If TRUE, forces Stan to recompile the model 
#'        executable. Useful if the underlying Stan code has changed. 
#'        Defaults to FALSE for speed.
#' @param ... Additional arguments passed to the cmdstanr fit methods 
#'        (e.g., iter, chains, parallel_chains, adapt_delta).
#' @importFrom cmdstanr cmdstan_model cmdstan_path set_cmdstan_path
#' @export
fit_glam <- function(data, method = "meanfield", s_val = 3e-4, 
                     force_recompile = FALSE, ...) {
  
  # 1. Verify CmdStan Environment
  # Checks if the CmdStan path is set; if not, attempts to find it.
  has_path <- tryCatch({
    dir.exists(cmdstanr::cmdstan_path())
  }, error = function(e) FALSE)
  
  if (!has_path) {
    tryCatch({
      cmdstanr::set_cmdstan_path()
    }, error = function(e) {
      stop("CmdStan path is not set. Please run cmdstanr::set_cmdstan_path().")
    })
  }
  
  # 2. Locate Stan Model Files
  # We look in the package installation folder first, then fallback to 
  # 'inst/stan' for active development (e.g., using devtools::load_all()).
  model_path <- system.file("stan", "hierarchical_glam.stan", package = "fractalGLAM")
  include_dir <- system.file("stan", "include", package = "fractalGLAM")
  
  if (model_path == "") {
    # Fallback for development environments
    model_path <- "inst/stan/hierarchical_glam.stan"
    include_dir <- "inst/stan/include"
    
    if (!file.exists(model_path)) {
      stop("Could not locate hierarchical_glam.stan. Ensure the package is installed or the path is correct.")
    }
  }
  
  # 3. Load/Compile Model
  # Compilation is skipped if the binary already exists and force_recompile is FALSE.
  mod <- cmdstanr::cmdstan_model(
    model_path,
    include_paths = include_dir,
    force_recompile = force_recompile
  )
  
  # 4. Prepare Data List for Stan
  # Translates the R data frame into the specific list format required by Stan.
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
    s_fixed = ifelse(is.null(s_val), -1.0, s_val) # -1 triggers estimation in Stan
  )
  
  # 5. Execute Inference
  if (method == "meanfield") {
    return(mod$variational(data = stan_list, ...))
  } else {
    return(mod$sample(data = stan_list, ...))
  }
}