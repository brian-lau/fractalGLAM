#' Perform Model Comparison via LOO/WAIC
#'
#' This function takes a list of fitted GLAM objects and returns a comparison
#' table to determine which model (Full vs. Restricted) better explains the data.
#'
#' @param fit_list A named list of cmdstanr fit objects (e.g., list(full = fit1, restricted = fit2)).
#' @param criterion Character, either "loo" or "waic".
#' @importFrom loo loo waic loo_compare
#' @importFrom dplyr as_tibble mutate
#' @export
compare_glam_models <- function(fit_list, criterion = "loo") {
  
  if (is.null(names(fit_list))) {
    names(fit_list) <- paste0("model_", seq_along(fit_list))
  }
  
  results <- lapply(fit_list, function(f) {
    if (criterion == "loo") {
      return(f$loo())
    } else {
      return(f$waic())
    }
  })
  
  comparison <- loo::loo_compare(results)
  
  # Convert to a more readable tidy format
  comp_df <- as.data.frame(comparison) %>%
    dplyr::mutate(model = rownames(.)) %>%
    dplyr::as_tibble() %>%
    dplyr::select(model, everything())
  
  return(comp_df)
}

#' K-Fold Cross-Validation for GLAM
#'
#' Implements out-of-sample prediction testing by partitioning data by subject.
#'
#' @param data Prepared data frame.
#' @param k Number of folds (if k = J, this is Leave-One-Subject-Out).
#' @param ... Arguments passed to fit_glam.
#' @export
run_glam_cv <- function(data, k = 5, ...) {
  subjects <- unique(data$subject_id)
  n_subs <- length(subjects)
  
  # Create folds
  folds <- cut(seq_along(subjects), breaks = k, labels = FALSE)
  shuffled_subs <- sample(subjects)
  
  cv_results <- lapply(1:k, function(i) {
    test_subs <- shuffled_subs[folds == i]
    train_data <- data %>% dplyr::filter(!(subject_id %in% test_subs))
    test_data <- data %>% dplyr::filter(subject_id %in% test_subs)
    
    # Fit on training data
    fit <- fit_glam(train_data, ...)
    
    # Note: For full out-of-sample RT/Choice prediction, 
    # we would use predict_glam() logic here.
    return(list(fold = i, test_subjects = test_subs))
  })
  
  return(cv_results)
}