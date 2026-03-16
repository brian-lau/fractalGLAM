#' Prepare N-Item Trial Data for GLAM
#' 
#' Generalizes preprocessing for an arbitrary number of items. Detects item
#' columns based on naming patterns (v_0, v_1... and g_0, g_1...).
#' 
#' @param df Data frame containing subject_id, rt, choice, and item columns.
#' @param scale_range Numeric vector of length 2 for value rescaling.
#' @importFrom dplyr mutate across select starts_with
#' @importFrom scales rescale
#' @export
prep_glam_data <- function(df, scale_range = c(1, 10)) {
  
  # 1. Detect number of items
  # We look for columns starting with 'v_' followed by a digit
  v_cols <- sort(colnames(df)[grepl("^v_\\d+$", colnames(df))])
  g_cols <- sort(colnames(df)[grepl("^g_\\d+$", colnames(df))])
  n_items <- length(v_cols)
  
  if (n_items < 2) {
    stop("Could not detect at least 2 item columns (v_0, v_1...). Check naming convention.")
  }
  
  if (length(v_cols) != length(g_cols)) {
    stop("Mismatch between number of value columns and gaze columns.")
  }
  
  # 2. Validate Z-scores (across all items)
  # Reference glambox assumes values are relative to the individual's mean
  all_values <- as.vector(as.matrix(df[, v_cols]))
  if (abs(mean(all_values, na.rm = TRUE)) > 0.1) {
    warning("Input values do not appear to be Z-scored. Standardizing now...")
    # Optional: implement auto-standardization here if preferred
  }
  
  # 3. Rescale and Normalize
  # Positive rescaling is required for the drift accumulation floor.
  df <- df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(v_cols), 
                                ~ scales::rescale(.x, to = scale_range)))
  
  # Gaze must sum to 1 across all N items for every trial.
  total_gaze <- rowSums(df[, g_cols])
  df[, g_cols] <- df[, g_cols] / total_gaze
  
  # 4. Final Structure Check
  required_base <- c("subject_id", "rt", "choice")
  missing <- setdiff(required_base, colnames(df))
  if (length(missing) > 0) stop(paste("Missing core columns:", paste(missing, collapse = ", ")))
  
  return(df)
}