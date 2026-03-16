#' Generate Predictions from a Fitted GLAM
#' 
#' This function generates simulated choices and response times using the 
#' parameters estimated by a GLAM fit. It can use either the posterior means 
#' for a point-estimate prediction or sample from the full posterior 
#' distribution to capture model uncertainty.
#' 
#' @param fit A cmdstanr fit object.
#' @param newdata A data frame prepared via \code{prep_glam_data}.
#' @param n_sims Number of simulations per trial.
#' @param type Character, either "mean" to use posterior means or 
#'        "posterior" to sample from the full posterior distribution.
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate summarise across starts_with left_join slice filter select n
#' @importFrom posterior as_draws_df
#' @export
predict_glam <- function(fit, newdata, n_sims = 10, type = "mean") {
  # Define variables for global check
  param_raw <- subject_idx <- param <- value <- value_right <- gaze_right <- 
    gamma <- gaze_left <- value_left <- v <- s <- signal <- drift <- NULL
  
  # 1. Extract Subject-Level Parameters
  # We extract the subject-indexed parameters (e.g., v[1], gamma[1])
  draws <- fit$draws(c("v", "gamma", "sigma", "s")) %>% 
    posterior::as_draws_df()
  
  subs <- unique(newdata$subject_id)
  
  if (type == "mean") {
    # Use the average of the posterior for each parameter (Point Estimate)
    param_summary <- draws %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = "param_raw", values_to = "value") %>%
      dplyr::filter(grepl("\\[", .data$param_raw)) %>% 
      dplyr::mutate(
        param = gsub("\\[.*", "", .data$param_raw),
        subject_idx = as.numeric(gsub(".*\\[|\\].*", "", .data$param_raw)),
        subject_id = subs[.data$subject_idx]
      ) %>%
      dplyr::select(-.data$param_raw, -.data$subject_idx) %>%
      tidyr::pivot_wider(names_from = .data$param, values_from = .data$value)
    
    pred_data <- newdata %>%
      dplyr::left_join(param_summary, by = "subject_id") %>%
      dplyr::slice(rep(1:dplyr::n(), each = n_sims))
    
  } else {
    # Sample n_sims random draws from the posterior for each trial
    # This captures the 'Total Uncertainty' of the model.
    draw_indices <- sample(1:nrow(draws), n_sims, replace = TRUE)
    
    pred_data <- lapply(draw_indices, function(i) {
      this_draw <- draws[i, ] %>%
        tidyr::pivot_longer(dplyr::everything(), names_to = "param_raw", values_to = "value") %>%
        dplyr::filter(grepl("\\[", .data$param_raw)) %>% 
        dplyr::mutate(
          param = gsub("\\[.*", "", .data$param_raw),
          subject_idx = as.numeric(gsub(".*\\[|\\].*", "", .data$param_raw)),
          subject_id = subs[.data$subject_idx]
        ) %>%
        dplyr::select(-.data$param_raw, -.data$subject_idx) %>%
        tidyr::pivot_wider(names_from = .data$param, values_from = .data$value)
      
      newdata %>% dplyr::left_join(this_draw, by = "subject_id")
    }) %>% dplyr::bind_rows()
  }
  
  # 2. Vectorized Simulation Engine
  # This logic MUST match the calculate_glam_drift function in Stan exactly.
  preds <- pred_data %>%
    dplyr::mutate(
      # A. Recalculate the relative decision signal
      signal = (.data$value_right * (.data$gaze_right + .data$gamma * .data$gaze_left)) - 
        (.data$value_left * (.data$gaze_left + .data$gamma * .data$gaze_right)),
      
      # B. Calculate Magnitude-based Drift (The 'Hump' Logic)
      # s is the scaling constant (e.g., 3e-4), v is velocity.
      # We add 0.01 to prevent infinite RTs at the point of indifference.
      drift_raw = .data$v * .data$s * (abs(.data$signal) + 0.01),
      drift = pmax(drift_raw, 1e-7), 
      
      # C. Simulate outcomes
      # Choice follows the sign of the signal
      pred_choice = stats::rbinom(dplyr::n(), 1, stats::plogis(.data$signal * 0.5)), 
      # RT follows the log-normal distribution defined in Stan
      pred_rt = stats::rlnorm(dplyr::n(), log(1 / .data$drift), .data$sigma)
    )
  
  return(preds)
}