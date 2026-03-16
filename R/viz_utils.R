#' @import ggplot2
#' @import dplyr
#' @importFrom stats glm coef
#' @importFrom rlang .data
NULL

#' Professional Theme for fractalGLAM
#' @keywords internal
theme_glam <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      legend.position = "bottom"
    )
}

#' Plot Choice Psychometric Curves
#' 
#' Visualizes the probability of choosing the right item as a function of 
#' the value difference, optionally split by gaze advantage.
#' 
#' @param data Prepared data frame (observed or predicted).
#' @param bin_gaze Logical. If TRUE, splits curves by gaze advantage.
#' @export
plot_choice_curves <- function(data, bin_gaze = FALSE) {
  delta_v <- choice <- gaze_diff <- gaze_bin <- subject_id <- NULL
  
  plot_df <- data %>%
    dplyr::mutate(
      delta_v = .data$value_right - .data$value_left,
      gaze_diff = .data$gaze_right - .data$gaze_left
    )
  
  if (bin_gaze) {
    plot_df <- plot_df %>%
      dplyr::mutate(gaze_bin = ifelse(gaze_diff > 0, "Gaze Right", "Gaze Left"))
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_v, y = choice, color = gaze_bin))
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_v, y = choice))
  }
  
  p + ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                           se = TRUE, alpha = 0.2) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(x = "Value Difference (R - L)", y = "Pr(Choose Right)",
                  title = "Choice Psychometric Function") +
    theme_glam() +
    ggplot2::facet_wrap(~subject_id)
}

#' Plot RT Chronometric Curves
#' 
#' @param data Prepared data frame (usually predictions).
#' @param observed_data Optional. The raw data frame used for fitting.
#' @param bin_gaze Logical. If TRUE, splits curves by gaze advantage.
#' @export
plot_rt_curves <- function(data, observed_data = NULL, bin_gaze = FALSE) {
  delta_v <- rt <- gaze_diff <- gaze_bin <- subject_id <- type <- NULL
  rt_col <- if ("pred_rt" %in% colnames(data)) "pred_rt" else "rt"
  
  plot_df <- data %>%
    dplyr::mutate(
      delta_v = .data$value_right - .data$value_left,
      gaze_diff = .data$gaze_right - .data$gaze_left
    )
  
  if (bin_gaze) {
    plot_df <- plot_df %>%
      dplyr::mutate(gaze_bin = ifelse(gaze_diff > 0, "Gaze Right", "Gaze Left"))
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_v, y = .data[[rt_col]], color = gaze_bin))
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = delta_v, y = .data[[rt_col]]))
  }
  
  if (!is.null(observed_data)) {
    obs_df <- observed_data %>% dplyr::mutate(delta_v = .data$value_right - .data$value_left)
    p <- p + ggplot2::geom_point(data = obs_df, ggplot2::aes(x = delta_v, y = rt), 
                                 alpha = 0.2, size = 1, inherit.aes = FALSE, color = "grey50")
  }
  
  p + ggplot2::geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(x = "Value Difference (R - L)", y = "Response Time (ms)",
                  title = "RT Chronometric Function") +
    theme_glam() +
    ggplot2::facet_wrap(~subject_id, scales = "free_y")
}

#' Plot Binned RT Fit Comparison (Posterior Predictive Check)
#' 
#' A high-level diagnostic plot that bins trials by value difference to compare 
#' observed RT 'humps' against model predictions.
#' 
#' @param observed Prepared data frame of actual trials.
#' @param predicted Data frame of simulations from \code{predict_glam}.
#' @param n_bins Number of bins for the delta_v axis.
#' @export
plot_rt_fit <- function(observed, predicted, n_bins = 10) {
  delta_val <- rt <- type <- m_delta <- m_rt <- se_rt <- bin <- NULL
  
  obs_clean <- observed %>% 
    dplyr::mutate(delta_val = .data$value_right - .data$value_left, type = "Observed")
  
  rt_col <- if ("pred_rt" %in% colnames(predicted)) "pred_rt" else "rt"
  pred_clean <- predicted %>% 
    dplyr::mutate(delta_val = .data$value_right - .data$value_left, 
                  rt = .data[[rt_col]], type = "Predicted")
  
  combined <- dplyr::bind_rows(
    dplyr::select(obs_clean, delta_val, rt, type),
    dplyr::select(pred_clean, delta_val, rt, type)
  ) %>%
    dplyr::mutate(bin = ggplot2::cut_number(delta_val, n_bins)) %>%
    dplyr::group_by(bin, type) %>%
    dplyr::summarise(m_delta = mean(delta_val), m_rt = mean(rt),
                     se_rt = sd(rt) / sqrt(dplyr::n()), .groups = "drop")
  
  ggplot2::ggplot(combined, ggplot2::aes(x = m_delta, y = m_rt, color = type, group = type)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = m_rt - se_rt, ymax = m_rt + se_rt), width = 0.1) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = c("Observed" = "#E41A1C", "Predicted" = "#377EB8")) +
    ggplot2::labs(x = "Value Difference (R - L)", y = "Mean RT (ms)",
                  title = "Binned RT Posterior Predictive Check") +
    theme_glam()
}

#' Plot Individual Gamma Estimates vs. Choice Slope
#' 
#' @param fit A cmdstanr fit object.
#' @param data The prepared data frame used for fitting.
#' @export
plot_gamma_performance <- function(fit, data) {
  gamma_est <- slope <- subject_id <- subject_idx <- delta_v <- gaze_diff <- NULL
  
  draws <- fit$draws("gamma") %>% 
    posterior::as_draws_df() %>%
    dplyr::summarise(dplyr::across(dplyr::starts_with("gamma["), mean)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "raw", values_to = "gamma_est") %>%
    dplyr::mutate(subject_idx = as.numeric(gsub(".*\\[|\\].*", "", .data$raw)))
  
  slopes <- data %>%
    dplyr::mutate(delta_v = .data$value_right - .data$value_left,
                  gaze_diff = .data$gaze_right - .data$gaze_left) %>%
    dplyr::group_by(.data$subject_id) %>%
    dplyr::summarise(
      slope = stats::coef(stats::glm(choice ~ base::scale(gaze_diff) + base::scale(delta_v), 
                                     family = "binomial"))[2]
    ) %>%
    dplyr::mutate(subject_idx = as.numeric(as.factor(.data$subject_id)))
  
  plot_df <- dplyr::left_join(slopes, draws, by = "subject_idx")
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$gamma_est, y = .data$slope)) +
    ggplot2::geom_smooth(method = "lm", color = "black", alpha = 0.2, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(color = factor(.data$subject_id)), size = 4) +
    ggplot2::labs(x = "Model Estimate (Gamma)", y = "Observed Gaze Bias Sensitivity",
                  title = "Model Validation: Estimates vs. Behavior") +
    theme_glam()
}

#' Visualize Gaze Bias Effect
#'
#' @param data Prepared data frame.
#' @export
plot_gaze_bias_effect <- function(data) {
  gaze_diff <- val_diff <- val_bin <- choice <- NULL
  data %>%
    dplyr::mutate(
      gaze_diff = .data$gaze_right - .data$gaze_left,
      val_diff = .data$value_right - .data$value_left,
      val_bin = ifelse(.data$val_diff > 0, "Right Higher Value", "Left Higher Value")
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$gaze_diff, y = .data$choice, color = .data$val_bin)) +
    ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    ggplot2::labs(x = "Gaze Difference (R - L)", y = "Pr(Choose Right)",
                  title = "Behavioral Gaze Bias Effect") +
    theme_glam()
}

#' Visualize Group-Level Parameter Posteriors
#' 
#' @param fit A cmdstanr fit object.
#' @export
plot_hierarchical_params <- function(fit) {
  params <- c("mu_v_raw", "mu_gamma_raw", "mu_sigma_raw")
  bayesplot::mcmc_areas(fit$draws(params), prob = 0.95) +
    ggplot2::labs(title = "Group-Level Hyper-parameter Posteriors",
                  subtitle = "95% Credible Intervals") +
    theme_glam()
}