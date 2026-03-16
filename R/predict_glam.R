#' Generate Predictions from an N-Item Fitted GLAM
#' 
#' Simulates a stochastic race between K accumulators for each trial using
#' the logistic drift mapping and Wiener process logic.
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
  
  # 2. Extract Subject-Level Parameters from the Fit
  # Note: These names match the 'transformed parameters' block in Stan
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
    # Full posterior sampling mode
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
  
  # 3. Race Simulation Engine (Wiener Process)
  # Boundary is fixed at 1.0 per glambox spec
  boundary <- 1.0
  dt <- 0.002 # Simulation time step (2ms for speed)
  
  preds <- sim_data %>%
    dplyr::mutate(sim_idx = 1:dplyr::n()) %>%
    dplyr::group_by(.data$sim_idx) %>%
    dplyr::mutate(
      # A. Calculate Drifts for all K items
      # Logic: Absolute Evidence A_k = v_k * (g_k + gamma * (1-g_k))
      drifts = list({
        v_vals <- unlist(dplyr::pick(dplyr::all_of(v_cols)))
        g_vals <- unlist(dplyr::pick(dplyr::all_of(g_cols)))
        
        A <- v_vals * (g_vals + .data$gamma_sub * (1 - g_vals))
        
        # Relative signal r_k = A_k - max(A_others)
        r <- sapply(1:K, function(k) A[k] - max(A[-k]))
        
        # Logistic drift mapping: d = v / (1 + exp(-s * r))
        .data$v_sub / (1 + exp(-.data$s_sub * r))
      }),
      
      # B. Perform Stochastic Race
      race_result = list({
        d_vec <- unlist(.data$drifts)
        noise_sd <- .data$sigma_sub * sqrt(dt)
        
        evidence <- rep(boundary/2, K) # Start at half-way (beta=0.5)
        time_steps <- 0
        max_steps <- 5.0 / dt # Cap at 5 seconds
        
        winner <- NA
        while(is.na(winner) && time_steps < max_steps) {
          time_steps <- time_steps + 1
          # Wiener step: Evidence_t+1 = Evidence_t + drift*dt + noise
          evidence <- evidence + d_vec * dt + rnorm(K, 0, noise_sd)
          
          # Check for boundary crossing
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