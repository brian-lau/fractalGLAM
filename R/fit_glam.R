#' Fit N-Item GLAM via Hierarchical Bayesian Inference
#' 
#' @param data Data frame prepared via \code{prep_glam_data}.
#' @param method Character, "meanfield" or "sample".
#' @param s_val Numeric. Fixed sensitivity. If NULL, estimated as free parameter.
#' @param force_recompile Logical. Defaults to FALSE.
#' @param ... Additional arguments for cmdstanr.
#' @export
fit_glam <- function(data, method = "meanfield", s_val = NULL, 
                     force_recompile = FALSE, ...) {
  
  # 1. Detect Item Columns
  v_cols <- sort(colnames(data)[grepl("^v_\\d+$", colnames(data))])
  g_cols <- sort(colnames(data)[grepl("^g_\\d+$", colnames(data))])
  K_items <- length(v_cols)
  
  # 2. Locate Stan Model
  model_path <- system.file("stan", "hierarchical_glam.stan", package = "fractalGLAM")
  include_dir <- system.file("stan", "include", package = "fractalGLAM")
  
  if (model_path == "") {
    model_path <- "inst/stan/hierarchical_glam.stan"
    include_dir <- "inst/stan/include"
  }
  
  mod <- cmdstanr::cmdstan_model(
    model_path,
    include_paths = include_dir,
    force_recompile = force_recompile
  )
  
  # 3. Prepare Data List (Multinomial Format)
  stan_list <- list(
    N = nrow(data),
    J = length(unique(data$subject_id)),
    K = K_items,
    subject = as.numeric(as.factor(data$subject_id)),
    choice = as.integer(data$choice),
    rt = as.numeric(data$rt),
    v = as.matrix(data[, v_cols]),
    g = as.matrix(data[, g_cols]),
    s_fixed = ifelse(is.null(s_val), -1.0, s_val),
    min_rt = min(data$rt) # ADD THIS LINE
  )
  
  # 4. Execute Inference
  if (method == "meanfield") {
    return(mod$variational(data = stan_list, ...))
  } else {
    return(mod$sample(data = stan_list, ...))
  }
}