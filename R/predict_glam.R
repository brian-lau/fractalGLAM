#' Generate Predictions from an N-Item Fitted GLAM
#' 
#' Simulates a stochastic race between K accumulators for each trial.
#' Matches the Stan model's "Start at Zero" logic and Dual-Path drift architecture.
#' 
#' @param fit A cmdstanr fit object.
#' @param newdata A data frame prepared via \code{prep_glam_data}.
#' @param n_sims Number of simulations per trial.
#' @param type Character, "mean" to use posterior means or "posterior" to sample.
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate summarise across starts_with left_join slice filter select n bind_rows pick
#' @importFrom posterior as_draws_df
#' @export
predict_glam <- function(fit, newdata, n_sims = 10, type = "mean") {
  # Define variables for global check
  param_raw <- subject_idx <- param <- value <- subject_id <- NULL
  
  # 1. Identify item columns
  v_cols <- sort(colnames(newdata)[grepl("^v_\\d+$", colnames(newdata))])
  g_cols <- sort(colnames(newdata)[grepl("^g_\\d+$", colnames(newdata))])
  K <- length(v_cols)
  
  # 2. Extract Subject-Level Parameters
  draws <- fit$draws(c("v_sub", "gamma_sub", "sigma_sub", "s_sub", "ndt_sub")) %>% 
    posterior::as_draws_df()
  
  subs <- unique(newdata$subject_id)
  
  if (type == "mean") {
    params <- draws %>%
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
    
    sim_data <- newdata %>%
      dplyr::left_join(params, by = "subject_id") %>%
      dplyr::slice(rep(1:dplyr::n(), each = n_sims))
  } else {
    draw_indices <- sample(1:nrow(draws), n_sims, replace = TRUE)
    sim_data <- lapply(draw_indices, function(i) {
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
  
  # 3. Race Simulation Engine
  boundary <- 1.0
  dt <- 0.002 
  
  preds <- sim_data %>%
    dplyr::mutate(sim_idx = 1:dplyr::n()) %>%
    dplyr::group_by(.data$sim_idx) %>%
    dplyr::mutate(
      drifts = list({
        v_vals <- unlist(dplyr::pick(dplyr::all_of(v_cols)))
        g_vals <- unlist(dplyr::pick(dplyr::all_of(g_cols)))
        
        # Absolute Evidence (A)
        A <- v_vals * (g_vals + .data$gamma_sub * (1 - g_vals))
        
        if (K == 2) {
          # Path A: Simple linear drift for 2-item stability
          .data$v_sub * A
        } else {
          # Path B: Scaled Logistic Race for N-item multinomial
          r <- sapply(1:K, function(k) A[k] - max(A[-k]))
          s_capped <- pmin(.data$s_sub, 10.0) #
          .data$v_sub * A / (1 + exp(-s_capped * r))
        }
      }),
      
      race_result = list({
        d_vec <- unlist(.data$drifts)
        noise_sd <- .data$sigma_sub * sqrt(dt)
        
        # START AT ZERO logic
        # In this single-boundary race, hitting zero does NOT end the trial.
        evidence <- rep(0.001, K) 
        time_steps <- 0
        max_steps <- 5.0 / dt 
        
        winner <- NA
        while(is.na(winner) && time_steps < max_steps) {
          time_steps <- time_steps + 1
          # Wiener step: Evidence_t+1 = Evidence_t + drift*dt + noise
          evidence <- evidence + d_vec * dt + rnorm(K, 0, noise_sd)
          
          # Check for boundary crossing (Top boundary only)
          if(any(evidence >= boundary)) {
            winner <- which.max(evidence >= boundary)
          }
        }
        
        list(
          pred_choice = ifelse(is.na(winner), sample(1:K, 1), winner),
          pred_rt = (time_steps * dt) + .data$ndt_sub
        )
      })
    ) %>%
    dplyr::mutate(
      pred_choice = .data$race_result[[1]]$pred_choice,
      pred_rt = .data$race_result[[1]]$pred_rt
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$drifts, -.data$race_result, -.data$sim_idx)
  
  return(preds)
}