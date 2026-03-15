#' @import ggplot2
#' @import dplyr
#' @importFrom stats glm coef
#' @importFrom rlang .data
NULL

#' Plot Choice Psychometric Curves
#' 
#' @param data Prepared data frame (observed or predicted).
#' @param bin_gaze Logical. If TRUE, splits curves by gaze advantage (Left vs Right).
#' @importFrom ggplot2 ggplot aes geom_smooth facet_wrap labs theme_minimal
#' @export
plot_choice_curves <- function(data, bin_gaze = FALSE) {
  delta_v <- choice <- gaze_diff <- gaze_bin <- subject_id <- NULL
  
  # Calculate predictors
  plot_df <- data %>%
    dplyr::mutate(
      delta_v = .data$value_right - .data$value_left,
      gaze_diff = .data$gaze_right - .data$gaze_left
    )
  
  # Optional gaze binning
  if (bin_gaze) {
    plot_df <- plot_df %>%
      dplyr::mutate(gaze_bin = ifelse(gaze_diff > 0, "Gaze Right", "Gaze Left"))
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_v, y = choice, color = gaze_bin))
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_v, y = choice))
  }
  
  # Build the plot
  p <- p +
    ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                         se = TRUE, alpha = 0.1) +
    ggplot2::labs(
      x = "Value Difference (Right - Left)",
      y = "Probability of Choosing Right",
      title = "Choice Psychometric Function",
      color = "Gaze Advantage"
    ) +
    ggplot2::theme_minimal()
  
  # Facet if multiple subjects exist
  if (length(unique(data$subject_id)) > 1) {
    p <- p + ggplot2::facet_wrap(~subject_id)
  }
  
  return(p)
}

#' Plot RT Chronometric Curves with Optional Raw Data
#' 
#' @param data Prepared data frame (usually predictions).
#' @param observed_data Optional. The raw data frame used for fitting.
#' @param bin_gaze Logical. If TRUE, splits curves by gaze advantage.
#' @importFrom ggplot2 ggplot aes geom_smooth geom_point facet_wrap labs theme_minimal
#' @export
plot_rt_curves <- function(data, observed_data = NULL, bin_gaze = FALSE) {
  delta_v <- rt <- gaze_diff <- gaze_bin <- subject_id <- type <- NULL
  
  # Identify RT column in prediction data
  rt_col <- if ("pred_rt" %in% colnames(data)) "pred_rt" else "rt"
  
  # Process prediction/main data
  plot_df <- data %>%
    dplyr::mutate(
      delta_v = .data$value_right - .data$value_left,
      gaze_diff = .data$gaze_right - .data$gaze_left,
      type = "Model Prediction"
    )
  
  # Base plot with the prediction curve
  if (bin_gaze) {
    plot_df <- plot_df %>%
      dplyr::mutate(gaze_bin = ifelse(gaze_diff > 0, "Gaze Right", "Gaze Left"))
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_v, y = .data[[rt_col]], color = gaze_bin))
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_v, y = .data[[rt_col]]))
  }
  
  # Add observed data as points if provided
  if (!is.null(observed_data)) {
    obs_df <- observed_data %>%
      dplyr::mutate(delta_v = .data$value_right - .data$value_left)
    
    # We add points first so they sit behind the trend line
    p <- p + ggplot2::geom_point(data = obs_df, ggplot2::aes(x = delta_v, y = rt), 
                                 alpha = 0.2, size = 1, inherit.aes = FALSE, color = "grey50")
  }
  
  # Add the smooth trend line
  p <- p +
    ggplot2::geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    ggplot2::labs(
      x = "Value Difference (Right - Left)",
      y = "Response Time",
      title = "RT Chronometric Function",
      subtitle = "Points = Observed Data; Lines = Model Fit",
      color = "Gaze Advantage"
    ) +
    ggplot2::theme_minimal()
  
  # Facet if multiple subjects exist
  if (length(unique(data$subject_id)) > 1) {
    p <- p + ggplot2::facet_wrap(~subject_id, scales = "free_y")
  }
  
  return(p)
}

#' Plot Individual Gamma Estimates vs. Choice Slope
#' 
#' @param fit A cmdstanr fit object.
#' @param data The prepared data frame used for fitting.
#' @importFrom stats glm coef
#' @importFrom dplyr mutate group_by summarise left_join across starts_with everything
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_smooth geom_point geom_hline labs theme_minimal
#' @export
plot_gamma_performance <- function(fit, data) {
  # Define variables for check
  gamma_est <- slope <- subject_id <- subject_idx <- value_left <- value_right <- delta_v <- choice <- NULL
  
  # 1. Extract Gamma Estimates
  draws <- fit$draws("gamma") %>% 
    posterior::as_draws_df() %>%
    dplyr::summarise(dplyr::across(dplyr::starts_with("gamma["), mean)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "raw", values_to = "gamma_est") %>%
    dplyr::mutate(subject_idx = as.numeric(gsub(".*\\[|\\].*", "", .data$raw)))
  
  # 2. Calculate Behavioral Slopes
  # FIX: use base::scale or just scale(), not stats::scale
  # slopes <- data %>%
  #   dplyr::mutate(delta_v = .data$value_right - .data$value_left) %>%
  #   dplyr::group_by(.data$subject_id) %>%
  #   dplyr::summarise(
  #     # We scale delta_v to make the logit slope interpretable on the plot
  #     slope = stats::coef(stats::glm(choice ~ base::scale(delta_v), family = "binomial"))[2]
  #   ) %>%
  #   dplyr::mutate(subject_idx = as.numeric(as.factor(.data$subject_id)))

  slopes <- data %>%
    mutate(
      delta_v = value_right - value_left,
      gaze_diff = gaze_right - gaze_left
    ) %>%
    group_by(subject_id) %>%
    summarise(
      # We now look at how Gaze Difference drives choice, 
      # which is the behavioral signature gamma is capturing.
      slope = stats::coef(stats::glm(choice ~ base::scale(gaze_diff) + base::scale(delta_v), 
                                     family = "binomial"))[2]
    ) %>%
    dplyr::mutate(subject_idx = as.numeric(as.factor(.data$subject_id)))
  
  
  # 3. Combine and Plot
  plot_df <- dplyr::left_join(slopes, draws, by = "subject_idx")
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$gamma_est, y = .data$slope)) +
    ggplot2::geom_smooth(method = "lm", color = "black", alpha = 0.2, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(color = factor(.data$subject_id)), size = 4) +
    ggplot2::labs(
      x = "Model Estimate (Gamma)",
      y = "Observed Gaze Bias Effect (Logistic Slope)",
      title = "Validation: Model Estimates vs. Actual Behavior",
      subtitle = "Higher Gamma should correlate with higher behavioral gaze sensitivity",
      color = "Subject"
    ) +
    ggplot2::theme_minimal()
}
#' #' Plot Individual Gamma Estimates vs. Choice Slope
#' #' 
#' #' @param fit A cmdstanr fit object.
#' #' @param data The prepared data frame used for fitting.
#' #' @export
#' plot_gamma_performance <- function(fit, data) {
#'   # Define variables for check
#'   gamma_est <- slope <- subject_id <- subject_idx <- value_left <- value_right <- delta_v <- choice <- NULL
#'   
#'   draws <- fit$draws("gamma") %>% 
#'     posterior::as_draws_df() %>%
#'     dplyr::summarise(dplyr::across(dplyr::starts_with("gamma["), mean)) %>%
#'     tidyr::pivot_longer(dplyr::everything(), names_to = "raw", values_to = "gamma_est") %>%
#'     dplyr::mutate(subject_idx = as.numeric(gsub(".*\\[|\\].*", "", .data$raw)))
#'   
#'   slopes <- data %>%
#'     dplyr::mutate(delta_v = .data$value_right - .data$value_left) %>%
#'     dplyr::group_by(.data$subject_id) %>%
#'     dplyr::summarise(
#'       slope = stats::coef(stats::glm(choice ~ delta_v, family = "binomial"))["delta_v"]
#'     ) %>%
#'     dplyr::mutate(subject_idx = as.numeric(as.factor(.data$subject_id)))
#'   
#'   plot_df <- dplyr::left_join(slopes, draws, by = "subject_idx")
#'   
#'   ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$gamma_est, y = .data$slope)) +
#'     ggplot2::geom_smooth(method = "lm", color = "black", alpha = 0.2) +
#'     ggplot2::geom_point(ggplot2::aes(color = factor(.data$subject_id)), size = 3) +
#'     ggplot2::theme_minimal()
#' }

#' Plot Predicted vs Observed Psychometric Functions
#' 
#' @param observed Data frame of observed data.
#' @param predicted Data frame of predicted data.
#' @export
plot_glam_fit <- function(observed, predicted) {
  delta_val <- rt <- pred_rt <- type <- value_right <- value_left <- NULL
  
  rt_obs <- observed %>%
    dplyr::mutate(delta_val = .data$value_right - .data$value_left) %>%
    dplyr::group_by(.data$delta_val) %>%
    dplyr::summarise(rt = mean(.data$rt), type = "Observed")
  
  rt_pred <- predicted %>%
    dplyr::mutate(delta_val = .data$value_right - .data$value_left) %>%
    dplyr::group_by(.data$delta_val) %>%
    dplyr::summarise(rt = mean(.data$pred_rt), type = "Predicted")
  
  dplyr::bind_rows(rt_obs, rt_pred) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$delta_val, y = .data$rt, color = .data$type)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal()
}

#' Visualize Gaze Bias Effect
#'
#' @param data Prepared data frame.
#' @export
plot_gaze_bias_effect <- function(data) {
  gaze_diff <- gaze_right <- gaze_left <- val_diff <- value_right <- value_left <- val_bin <- choice <- NULL
  
  data %>%
    dplyr::mutate(
      gaze_diff = .data$gaze_right - .data$gaze_left,
      val_diff = .data$value_right - .data$value_left,
      val_bin = ifelse(.data$val_diff > 0, "Right Higher Value", "Left Higher Value")
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$gaze_diff, y = .data$choice, color = .data$val_bin)) +
    ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    ggplot2::theme_minimal()
}

#' Plot Hierarchical Parameter Posteriors
#' 
#' @param fit A cmdstanr fit object.
#' @export
plot_hierarchical_params <- function(fit) {
  params_to_plot <- c("mu_v_raw", "mu_gamma_raw", "mu_sigma_raw")
  bayesplot::mcmc_areas(fit$draws(params_to_plot), prob = 0.95)
}