test_that("Data preprocessing handles Z-scores and normalization", {
  test_data <- data.frame(
    subject_id = rep(1, 5),
    value_left = c(-1, 0, 1, 0.5, -0.5), # Z-scored
    gaze_left = c(0.4, 0.5, 0.6, 0.1, 0.9),
    gaze_right = c(0.4, 0.5, 0.4, 0.8, 0.1)
  )
  
  # Check for error on non-Z-scored data
  bad_data <- test_data
  bad_data$value_left <- bad_data$value_left + 10
  expect_error(prep_glam_data(bad_data))
  
  # Check normalization
  clean_data <- prep_glam_data(test_data)
  expect_equal(all(clean_data$gaze_left + clean_data$gaze_right == 1), TRUE)
})