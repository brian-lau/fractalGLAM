#' Generate Predictions from a Fitted GLAM
#' @param fit A cmdstanr fit object.
#' @param newdata A data frame prepared via prep_glam_data.
#' @param n_sims Number of simulations per trial.
#' @param type Character, "mean" to use posterior means or "posterior" to sample.
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate summarise across starts_with left_join slice filter select n
#' @export
predict_glam <- function(fit, newdata, n_sims = 10, type = "mean") {
  # Define variables for check
  param_raw <- subject_idx <- param <- value <- value_right <- gaze_right <- 
    gamma <- gaze_left <- value_left <- v <- s <- signal <- drift <- NULL
  
  # 1. Extract subject-level parameters only
  # We look specifically for parameters with brackets [ ]
  draws <- fit$draws(c("v", "gamma", "sigma", "s")) %>% 
    posterior::as_draws_df()
  
  subs <- unique(newdata$subject_id)
  
  param_summary <- draws %>%
    dplyr::summarise(across(everything(), mean)) %>%
    tidyr::pivot_longer(everything(), names_to = "param_raw", values_to = "value") %>%
    # FIX: Filter for indexed parameters only to avoid 'NAs introduced by coercion'
    dplyr::filter(grepl("\\[", .data$param_raw)) %>% 
    mutate(
      param = gsub("\\[.*", "", .data$param_raw),
      subject_idx = as.numeric(gsub(".*\\[|\\].*", "", .data$param_raw)),
      subject_id = subs[.data$subject_idx]
    ) %>%
    dplyr::select(-.data$param_raw, -.data$subject_idx) %>%
    tidyr::pivot_wider(names_from = .data$param, values_from = .data$value)
  
  # 2. Vectorized simulation
  preds <- newdata %>%
    dplyr::left_join(param_summary, by = "subject_id") %>%
    dplyr::slice(rep(1:dplyr::n(), each = n_sims)) %>%
    mutate(
      # Recalculate signal using estimated gamma
      signal = (.data$value_right * (.data$gaze_right + .data$gamma * .data$gaze_left)) - 
        (.data$value_left * (.data$gaze_left + .data$gamma * .data$gaze_right)),
      # Calculate drift based on estimated v and s
      drift = .data$v * .data$s * (1 / (1 + exp(-.data$signal))),
      # Simulate outcomes
      pred_choice = stats::rbinom(dplyr::n(), 1, stats::plogis(.data$drift * 500)), # Matching Stan scale
      pred_rt = stats::rlnorm(dplyr::n(), log(1 / .data$drift), .data$sigma)
    )
  
  return(preds)
}