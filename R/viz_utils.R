#' @import ggplot2
#' @import dplyr
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

#' Plot N-Item Choice Psychometric Curves
#' 
#' Visualizes the probability of choosing each item as a function of its 
#' relative value (Value_i - max(Value_others)).
#' 
#' @param data Prepared data frame (observed or predicted).
#' @export
plot_choice_curves <- function(data) {
  # Identify item columns
  v_cols <- sort(colnames(data)[grepl("^v_\\d+$", colnames(data))])
  K <- length(v_cols)
  
  # Calculate relative values for each item
  plot_list <- lapply(1:K, function(k) {
    this_v <- data[[v_cols[k]]]
    other_v <- apply(data[, v_cols[-k]], 1, max)
    
    data.frame(
      subject_id = data$subject_id,
      rel_value = this_v - other_v,
      chosen = as.numeric(data$choice == k),
      item = paste0("Item ", k-1)
    )
  })
  
  plot_df <- dplyr::bind_rows(plot_list)
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$rel_value, y = .data$chosen, color = .data$item)) +
    ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                         se = TRUE, alpha = 0.1) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(x = "Relative Value (Item_i - max_others)", y = "Pr(Choice)",
                  title = "Multinomial Choice Functions",
                  color = "Option") +
    theme_glam() +
    ggplot2::facet_wrap(~subject_id)
}

#' Plot Binned RT Fit Comparison (N-Item Race)
#' 
#' @param observed Prepared data frame of actual trials.
#' @param predicted Data frame of simulations from \code{predict_glam}.
#' @param n_bins Number of bins for the RT distribution.
#' @export
plot_rt_fit <- function(observed, predicted, n_bins = 15) {
  type <- rt <- m_rt <- se_rt <- NULL
  
  obs_clean <- observed %>% dplyr::mutate(type = "Observed")
  
  rt_col <- if ("pred_rt" %in% colnames(predicted)) "pred_rt" else "rt"
  pred_clean <- predicted %>% 
    dplyr::mutate(rt = .data[[rt_col]], type = "Predicted")
  
  combined <- dplyr::bind_rows(
    dplyr::select(obs_clean, rt, type),
    dplyr::select(pred_clean, rt, type)
  )
  
  ggplot2::ggplot(combined, ggplot2::aes(x = rt, fill = type)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), 
                            bins = n_bins, alpha = 0.5, position = "identity") +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::scale_fill_manual(values = c("Observed" = "#E41A1C", "Predicted" = "#377EB8")) +
    ggplot2::labs(x = "Response Time (s)", y = "Density",
                  title = "RT Posterior Predictive Check",
                  subtitle = "Race Model Predicted vs. Observed Densities") +
    theme_glam()
}

#' Visualize Gaze Bias Effect (N-Item)
#'
#' @param data Prepared data frame.
#' @export
plot_gaze_bias_effect <- function(data) {
  # Using Gaze Advantage: Gaze_i - average(Gaze_others)
  g_cols <- sort(colnames(data)[grepl("^g_\\d+$", colnames(data))])
  K <- length(g_cols)
  
  plot_list <- lapply(1:K, function(k) {
    this_g <- data[[g_cols[k]]]
    other_g <- rowMeans(data[, g_cols[-k], drop = FALSE])
    
    data.frame(
      gaze_adv = this_g - other_g,
      chosen = as.numeric(data$choice == k),
      item = paste0("Item ", k-1)
    )
  })
  
  plot_df <- dplyr::bind_rows(plot_list)
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$gaze_adv, y = .data$chosen)) +
    ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "black") +
    ggplot2::labs(x = "Gaze Advantage (Item_i - mean_others)", y = "Pr(Choice)",
                  title = "Attentional Bias Effect") +
    theme_glam()
}