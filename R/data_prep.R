#' Prepare Trial Data for GLAM
#' 
#' This function performs essential cleaning and normalization required by the 
#' Gaze-weighted Linear Accumulator Model. It validates that values are Z-scored,
#' rescales them to a positive range for the drift process, and ensures gaze
#' proportions sum to 1.
#' 
#' @param df A data frame containing trial data. Must include: subject_id, rt, 
#'        choice, value_left, value_right, gaze_left, gaze_right.
#' @param scale_range Numeric vector of length 2. The GLAM typically requires 
#'        positive values (e.g., 1 to 10) to define the drift rate.
#' @importFrom dplyr mutate across starts_with select
#' @importFrom scales rescale
#' @importFrom rlang .data
#' @export
prep_glam_data <- function(df, scale_range = c(1, 10)) {
  # Define variables for global check
  total_gaze <- gaze_left <- gaze_right <- NULL
  
  # 1. Check for required columns
  required_cols <- c("subject_id", "rt", "choice", "value_left", 
                     "value_right", "gaze_left", "gaze_right")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(paste("The input data frame is missing required columns:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # 2. Validate Z-scores
  # GLAM requires values to be relative to the subject's mean to properly 
  # estimate the 'v' (velocity) parameter across a population.
  mean_val <- mean(df$value_left, na.rm = TRUE)
  if (abs(mean_val) > 0.1) {
    stop(paste("Input values must be Z-score standardized per subject (mean approx 0).", 
               "Found mean:", round(mean_val, 3)))
  }
  
  # 3. Rescale to positive range
  # The accumulator process requires a positive evidence 'floor'.
  df <- df %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("value_"), 
                                ~ scales::rescale(.x, to = scale_range)))
  
  # 4. Gaze Normalization (Relative Gaze)
  # Ensures that the sum of gaze proportions for any given trial equals 1.0.
  df <- df %>%
    dplyr::mutate(total_gaze = .data$gaze_left + .data$gaze_right,
                  gaze_left = .data$gaze_left / .data$total_gaze,
                  gaze_right = .data$gaze_right / .data$total_gaze) %>%
    dplyr::select(-.data$total_gaze)
  
  return(df)
}