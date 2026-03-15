#' Prepare Trial Data for GLAM
#' @param df A data frame containing trial data.
#' @param scale_range Numeric vector of length 2 for rescaling Z-scores.
#' @importFrom dplyr mutate across starts_with
#' @importFrom scales rescale
#' @importFrom rlang .data
#' @export
#'
prep_glam_data <- function(df, scale_range = c(1, 10)) {
  # Define variables for check
  total_gaze <- gaze_left <- gaze_right <- NULL
  # 1. Validate Z-scores
  mean_val <- mean(df$value_left, na.rm = TRUE)
  if (abs(mean_val) > 0.1) {
    stop(paste("Input values must be Z-score standardized per subject. Found mean:", round(mean_val, 3)))
  }
  
  # 2. Rescale to positive range for the accumulator
  df <- df %>%
    mutate(across(starts_with("value_"), 
                  ~ scales::rescale(.x, to = scale_range)))
  
  # 3. Gaze Normalization (Relative Gaze)
  df <- df %>%
    mutate(total_gaze = gaze_left + gaze_right,
           gaze_left = gaze_left / total_gaze,
           gaze_right = gaze_right / total_gaze) %>%
    dplyr::select(-total_gaze)
  
  return(df)
}