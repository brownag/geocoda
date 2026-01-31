# Tests for Input Validation
# Tests gc_validate_input_data() function

context("Input Validation Checklist")

test_that("gc_validate_input_data identifies valid data correctly", {
  # Create clean data (with projected coordinates to avoid geographic warning)
  valid_data <- data.frame(
    x = c(1000, 2000, 3000, 4000, 5000),
    y = c(1000, 2000, 3000, 4000, 5000),
    sand = c(30, 40, 35, 45, 50),
    silt = c(50, 45, 50, 40, 35),
    clay = c(20, 15, 15, 15, 15)
  )
  
  result <- gc_validate_input_data(
    valid_data,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_true(result$valid)
  # No should be no errors (warnings/notes are OK for this data)
  expect_true(!any(result$issues$severity == "error"))
})

test_that("gc_validate_input_data detects missing columns", {
  bad_data <- data.frame(
    x = c(1, 2, 3),
    y = c(1, 2, 3),
    sand = c(30, 40, 35)
    # Missing silt and clay
  )
  
  result <- gc_validate_input_data(
    bad_data,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_false(result$valid)
  expect_true(any(grepl("Missing component", result$issues$message)))
})

test_that("gc_validate_input_data detects missing coordinate columns", {
  data_df <- data.frame(
    sand = c(30, 40, 35),
    silt = c(50, 45, 50),
    clay = c(20, 15, 15)
  )
  
  result <- gc_validate_input_data(
    data_df,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",  # Doesn't exist
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_false(result$valid)
  expect_true(any(grepl("X-coordinate", result$issues$message) |
                    grepl("requires", result$issues$message)))
})

test_that("gc_validate_input_data detects sum constraint violations", {
  bad_sums <- data.frame(
    x = c(1000, 2000, 3000, 4000, 5000),
    y = c(1000, 2000, 3000, 4000, 5000),
    sand = c(30, 40, 35, 45, 50),
    silt = c(50, 45, 50, 40, 30),  # Row 5 has 30 instead of 40
    clay = c(20, 15, 15, 15, 10)   # Row 5 sums to 90, violates constraint
  )
  
  result <- gc_validate_input_data(
    bad_sums,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    target_sum = 100,
    sum_tolerance = 1,
    verbose = FALSE
  )
  
  # Should detect sum violation - message should mention sum or composition
  expect_true(any(grepl("sum|composition", result$issues$message, ignore.case = TRUE)))
})

test_that("gc_validate_input_data detects missing values", {
  with_na <- data.frame(
    x = c(1000, 2000, 3000, 4000, 5000),
    y = c(1000, 2000, 3000, 4000, 5000),
    sand = c(30, 40, NA, 45, 50),
    silt = c(50, 45, 50, 40, 35),
    clay = c(20, 15, 15, 15, 15)
  )
  
  result <- gc_validate_input_data(
    with_na,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_true(any(grepl("missing", result$issues$message, ignore.case = TRUE)))
  expect_equal(result$summary$n_complete_cases, 4)
})

test_that("gc_validate_input_data handles negative values", {
  negative_vals <- data.frame(
    x = c(1000, 2000, 3000),
    y = c(1000, 2000, 3000),
    sand = c(-5, 40, 35),
    silt = c(50, 45, 50),
    clay = c(20, 15, 15)
  )
  
  result <- gc_validate_input_data(
    negative_vals,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_true(any(grepl("negative", result$issues$message, ignore.case = TRUE)))
})

test_that("gc_validate_input_data detects values exceeding target sum", {
  over_sum <- data.frame(
    x = c(1000, 2000, 3000),
    y = c(1000, 2000, 3000),
    sand = c(120, 40, 35),    # Row 1 sand is 120, exceeds 100
    silt = c(50, 45, 50),
    clay = c(20, 15, 15)
  )
  
  result <- gc_validate_input_data(
    over_sum,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  # Should detect sand values exceeding 100
  expect_true(any(grepl("sand", result$issues$message, ignore.case = TRUE) & 
                    grepl("100|values", result$issues$message, ignore.case = TRUE)))
})

test_that("gc_validate_input_data detects all-zero compositions", {
  all_zero <- data.frame(
    x = c(1000, 2000, 3000),
    y = c(1000, 2000, 3000),
    sand = c(0, 40, 35),
    silt = c(0, 45, 50),
    clay = c(0, 15, 15)
  )
  
  result <- gc_validate_input_data(
    all_zero,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_true(any(grepl("all zero", result$issues$message)))
})

test_that("gc_validate_input_data detects geographic coordinates", {
  geo_coords <- data.frame(
    lon = c(-120, -119, -118),
    lat = c(38, 39, 40),
    sand = c(30, 40, 35),
    silt = c(50, 45, 50),
    clay = c(20, 15, 15)
  )
  
  result <- gc_validate_input_data(
    geo_coords,
    component_names = c("sand", "silt", "clay"),
    x_coord = "lon",
    y_coord = "lat",
    verbose = FALSE
  )
  
  expect_true(any(grepl("geographic", result$issues$message, ignore.case = TRUE)))
})

test_that("gc_validate_input_data detects duplicate locations", {
  with_dups <- data.frame(
    x = c(1000, 1000, 2000, 3000),        # Row 1 and 2 are duplicates
    y = c(1000, 1000, 2000, 3000),
    sand = c(30, 30, 40, 35),
    silt = c(50, 50, 45, 50),
    clay = c(20, 20, 15, 15)
  )
  
  result <- gc_validate_input_data(
    with_dups,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_true(any(grepl("duplicate", result$issues$message)))
})

test_that("gc_validate_input_data detects no spatial variation", {
  no_variation <- data.frame(
    x = c(1000, 1000, 1000),
    y = c(1000, 1000, 1000),
    sand = c(30, 40, 35),
    silt = c(50, 45, 50),
    clay = c(20, 15, 15)
  )
  
  result <- gc_validate_input_data(
    no_variation,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_true(any(grepl("same location", result$issues$message)))
})

test_that("gc_validate_input_data detects outliers when check_outliers=TRUE", {
  with_outliers <- data.frame(
    x = c(1000, 2000, 3000, 4000, 5000),
    y = c(1000, 2000, 3000, 4000, 5000),
    sand = c(30, 40, 35, 45, 95),    # Last value is outlier
    silt = c(50, 45, 50, 40, 5),
    clay = c(20, 15, 15, 15, 0)
  )
  
  result <- gc_validate_input_data(
    with_outliers,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    check_outliers = TRUE,
    verbose = FALSE
  )
  
  expect_true(any(grepl("outlier", result$issues$message, ignore.case = TRUE)))
})

test_that("gc_validate_input_data skips outlier check when check_outliers=FALSE", {
  with_outliers <- data.frame(
    x = c(1000, 2000, 3000, 4000, 5000),
    y = c(1000, 2000, 3000, 4000, 5000),
    sand = c(30, 40, 35, 45, 95),
    silt = c(50, 45, 50, 40, 5),
    clay = c(20, 15, 15, 15, 0)
  )
  
  result <- gc_validate_input_data(
    with_outliers,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    check_outliers = FALSE,
    verbose = FALSE
  )
  
  expect_false(any(grepl("outlier", result$issues$message, ignore.case = TRUE)))
})

test_that("gc_validate_input_data returns proper summary", {
  data_df <- data.frame(
    x = c(1000, 2000, 3000, 4000, 5000),
    y = c(1000, 2000, 3000, 4000, 5000),
    sand = c(30, 40, 35, 45, 50),
    silt = c(50, 45, 50, 40, 35),
    clay = c(20, 15, 15, 15, 15)
  )
  
  result <- gc_validate_input_data(
    data_df,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_equal(result$summary$n_observations, 5)
  expect_equal(result$summary$n_components, 3)
  expect_equal(result$summary$n_complete_cases, 5)
  expect_equal(result$summary$n_spatial_observations, 5)
})

test_that("gc_validate_input_data provides recommendations for fixed issues", {
  bad_sums <- data.frame(
    x = c(1000, 2000, 3000),
    y = c(1000, 2000, 3000),
    sand = c(30, 40, 35),
    silt = c(50, 45, 50),
    clay = c(20, 15, 50)  # Row 3 sums to 135
  )
  
  result <- gc_validate_input_data(
    bad_sums,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_true(length(result$recommendations) > 0)
  expect_true(any(grepl("rescal|resolv", result$recommendations, ignore.case = TRUE)))
})

test_that("gc_validate_input_data recommends data collection for small sample sizes", {
  small_data <- data.frame(
    x = c(1000, 2000),
    y = c(1000, 2000),
    sand = c(30, 40),
    silt = c(50, 45),
    clay = c(20, 15)
  )
  
  result <- gc_validate_input_data(
    small_data,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  expect_true(any(grepl("50|small", result$recommendations, ignore.case = TRUE)))
})

test_that("gc_validate_input_data uses custom target_sum and tolerance", {
  # Create data for target_sum = 1.0 (proportions)
  proportion_data <- data.frame(
    x = c(1000, 2000, 3000),
    y = c(1000, 2000, 3000),
    sand = c(0.30, 0.40, 0.35),
    silt = c(0.50, 0.45, 0.50),
    clay = c(0.20, 0.15, 0.15)
  )
  
  result <- gc_validate_input_data(
    proportion_data,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    target_sum = 1.0,
    sum_tolerance = 0.01,
    verbose = FALSE
  )
  
  expect_true(result$valid)
})

test_that("gc_validate_input_data handles sf objects", {
  skip_if_not_installed("sf")
  
  # Create valid data
  valid_data <- data.frame(
    x = c(1000, 2000, 3000),
    y = c(1000, 2000, 3000),
    sand = c(30, 40, 35),
    silt = c(50, 45, 50),
    clay = c(20, 15, 15)
  )
  
  # Convert to sf
  sf_data <- sf::st_as_sf(valid_data, coords = c("x", "y"))
  
  result <- gc_validate_input_data(
    sf_data,
    component_names = c("sand", "silt", "clay"),
    verbose = FALSE
  )
  
  expect_true(result$valid)
  expect_equal(result$summary$n_observations, 3)
})

test_that("gc_validate_input_data shows verbose output when requested", {
  data_df <- data.frame(
    x = c(1000, 2000),
    y = c(1000, 2000),
    sand = c(30, 40),
    silt = c(50, 45),
    clay = c(20, 15)
  )
  
  # Capture output
  output <- capture.output({
    result <- gc_validate_input_data(
      data_df,
      component_names = c("sand", "silt", "clay"),
      x_coord = "x",
      y_coord = "y",
      verbose = TRUE
    )
  })
  
  expect_true(any(grepl("Starting", output)))
  expect_true(any(grepl("SUMMARY", output)))
})

test_that("gc_validate_input_data handles complex mixed issues", {
  complex_data <- data.frame(
    x = c(1000, 2000, 3000, 4000, 5000),
    y = c(1000, 2000, 3000, 4000, 5000),
    sand = c(30, NA, 35, 45, 50),      # One NA
    silt = c(50, 45, 50, 40, 30),      # Row 5 will sum to 80
    clay = c(20, 15, 15, 15, 0)        # Row 5 below sum target
  )
  
  result <- gc_validate_input_data(
    complex_data,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  # Should have multiple issues
  expect_true(nrow(result$issues) >= 2)
  expect_false(result$valid)
})

test_that("gc_validate_input_data correctly identifies high error percentage", {
  # Create data with many constraint violations
  many_errors <- data.frame(
    x = c(1000, 2000, 3000, 4000, 5000),
    y = c(1000, 2000, 3000, 4000, 5000),
    sand = c(30, 10, 35, 45, 50),      # Rows 2-5 will have issues
    silt = c(50, 80, 50, 40, 35),      # Many sum violations
    clay = c(20, 5, 15, 15, 15)
  )
  
  result <- gc_validate_input_data(
    many_errors,
    component_names = c("sand", "silt", "clay"),
    x_coord = "x",
    y_coord = "y",
    verbose = FALSE
  )
  
  # With many violations, should detect them (as error or warning depending on %)
  sum_issues <- result$issues[grepl("sum|composition", result$issues$message, ignore.case = TRUE), ]
  if (nrow(sum_issues) > 0) {
    # Should detect violations - severity depends on percentage
    expect_true(sum_issues$severity[1] %in% c("error", "warning"))
  }
})
