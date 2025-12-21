test_that("gc_assess_stationarity returns list with expected structure", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = rnorm(n, mean = 0.5, sd = 0.8),
    ilr2 = rnorm(n, mean = -0.2, sd = 0.6)
  )

  result <- gc_assess_stationarity(data, method = "biplot", plot = FALSE)

  expect_type(result, "list")
  expect_true(all(c("stationary", "method_used", "pca_loadings", "pca_scores", 
                     "spatial_variance", "recommendation", "summary") %in% names(result)))
  expect_type(result$stationary, "logical")
  expect_equal(result$method_used, "biplot")
  expect_type(result$recommendation, "character")
})

test_that("gc_assess_stationarity local method works", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = rnorm(n, mean = 0.5, sd = 0.8),
    ilr2 = rnorm(n, mean = -0.2, sd = 0.6)
  )

  result <- gc_assess_stationarity(data, method = "local", plot = FALSE)

  expect_type(result, "list")
  expect_equal(result$method_used, "local")
  expect_type(result$stationary, "logical")
})

test_that("gc_assess_stationarity works with sf objects", {
  skip_if_not_installed("sf")
  set.seed(42)
  n <- 50
  data <- sf::st_as_sf(
    data.frame(
      x = runif(n, 0, 100),
      y = runif(n, 0, 100),
      ilr1 = rnorm(n, mean = 0.5, sd = 0.8),
      ilr2 = rnorm(n, mean = -0.2, sd = 0.6),
      id = 1:n
    ),
    coords = c("x", "y"),
    crs = "local"
  )
  
  # Remove geometry to convert back to data frame, keeping x/y as regular columns
  data_df <- as.data.frame(data)
  data_df$x <- sf::st_coordinates(data)[, 1]
  data_df$y <- sf::st_coordinates(data)[, 2]

  result <- gc_assess_stationarity(data_df, method = "biplot", plot = FALSE)

  expect_type(result$stationary, "logical")
  expect_equal(nrow(result$pca_scores), n)
})

test_that("gc_assess_stationarity requires at least 2 ILR dimensions", {
  data <- data.frame(
    x = c(1, 2, 3),
    y = c(4, 5, 6),
    ilr1 = c(0.1, 0.2, 0.3)
  )

  expect_error(gc_assess_stationarity(data), "at least 2 ILR columns")
})

test_that("gc_assess_gaussianity anderson method works", {
  samples <- data.frame(
    sand = c(20, 25, 30, 22, 18, 35, 28, 32, 26, 24),
    silt = c(60, 55, 50, 58, 62, 45, 52, 48, 54, 56),
    clay = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20)
  )

  ilr_vals <- compositions::ilr(compositions::acomp(samples))

  result <- gc_assess_gaussianity(ilr_vals, method = "anderson", plot = FALSE)

  expect_type(result, "list")
  expect_true(all(c("gaussian", "method_used", "p_values", "recommendation", "summary") %in% names(result)))
  expect_type(result$gaussian, "logical")
  expect_equal(result$method_used, "anderson")
  expect_equal(length(result$p_values), 2)
})

test_that("gc_assess_gaussianity shapiro method works", {
  samples <- data.frame(
    sand = c(20, 25, 30, 22, 18, 35, 28, 32, 26, 24),
    silt = c(60, 55, 50, 58, 62, 45, 52, 48, 54, 56),
    clay = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20)
  )

  ilr_vals <- compositions::ilr(compositions::acomp(samples))

  result <- gc_assess_gaussianity(ilr_vals, method = "shapiro", plot = FALSE)

  expect_type(result$gaussian, "logical")
  expect_equal(result$method_used, "shapiro")
})

test_that("gc_assess_gaussianity mardia method works", {
  samples <- data.frame(
    sand = c(20, 25, 30, 22, 18, 35, 28, 32, 26, 24),
    silt = c(60, 55, 50, 58, 62, 45, 52, 48, 54, 56),
    clay = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20)
  )

  ilr_vals <- compositions::ilr(compositions::acomp(samples))

  result <- gc_assess_gaussianity(ilr_vals, method = "mardia", plot = FALSE)

  expect_type(result$gaussian, "logical")
  expect_equal(result$method_used, "mardia")
  expect_true(!is.na(result$skewness))
  expect_true(!is.na(result$kurtosis))
})

test_that("gc_assess_gaussianity with data frame input", {
  samples <- data.frame(
    sand = c(20, 25, 30, 22, 18, 35, 28, 32, 26, 24),
    silt = c(60, 55, 50, 58, 62, 45, 52, 48, 54, 56),
    clay = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20)
  )

  ilr_vals <- compositions::ilr(compositions::acomp(samples))

  result <- gc_assess_gaussianity(as.data.frame(ilr_vals), method = "anderson", plot = FALSE)

  expect_type(result$gaussian, "logical")
})

test_that("gc_assess_gaussianity produces recommendations", {
  samples <- data.frame(
    sand = c(20, 25, 30, 22, 18, 35, 28, 32, 26, 24),
    silt = c(60, 55, 50, 58, 62, 45, 52, 48, 54, 56),
    clay = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20)
  )

  ilr_vals <- compositions::ilr(compositions::acomp(samples))

  result <- gc_assess_gaussianity(ilr_vals, method = "anderson", plot = FALSE)

  expect_type(result$recommendation, "character")
  expect_gt(nchar(result$recommendation), 0)
})

test_that("gc_apply_anamorphosis transforms simulated values", {
  set.seed(42)
  
  ref_data <- matrix(rnorm(100, 0.5, 0.8), nrow = 50, ncol = 2)
  sim_data <- matrix(rnorm(200, 0, 1), nrow = 100, ncol = 2)
  
  colnames(ref_data) <- c("ilr1", "ilr2")
  colnames(sim_data) <- c("ilr1", "ilr2")

  result <- gc_apply_anamorphosis(sim_data, ref_data, despike = TRUE)

  expect_equal(dim(result), dim(sim_data))
  expect_equal(colnames(result), colnames(sim_data))
  
  # Check that values are within reference range
  ref_range <- apply(ref_data, 2, range)
  for (i in seq_len(ncol(result))) {
    expect_true(all(result[, i] >= ref_range[1, i] & result[, i] <= ref_range[2, i]))
  }
})

test_that("gc_apply_anamorphosis preserves dimension count", {
  set.seed(42)
  
  ref_data <- matrix(rnorm(90, 0, 1), nrow = 30, ncol = 3)
  sim_data <- matrix(rnorm(150, 0, 1), nrow = 50, ncol = 3)

  result <- gc_apply_anamorphosis(sim_data, ref_data)

  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), nrow(sim_data))
})

test_that("gc_apply_anamorphosis checks dimension compatibility", {
  ref_data <- matrix(rnorm(50, 0, 1), nrow = 25, ncol = 2)
  sim_data <- matrix(rnorm(90, 0, 1), nrow = 30, ncol = 3)

  expect_error(gc_apply_anamorphosis(sim_data, ref_data), 
               "same number of dimensions")
})

test_that("gc_apply_anamorphosis despike works", {
  set.seed(42)
  
  ref_data <- matrix(rnorm(100, 0, 1), nrow = 50, ncol = 2)
  # Create mostly normal data with extreme values at the end
  extreme_vals <- matrix(rnorm(200, 0, 1), nrow = 100, ncol = 2)
  extreme_vals[99, 1] <- -10
  extreme_vals[100, 2] <- 10
  sim_data <- extreme_vals

  result <- gc_apply_anamorphosis(sim_data, ref_data, despike = TRUE)
  
  # Check output dimensions
  expect_equal(dim(result), dim(sim_data))
  
  # Check that result is a matrix
  expect_true(is.matrix(result))
})

test_that("Integration: stationarity -> gaussianity -> anamorphosis", {
  set.seed(42)
  
  # Create sample data with spatial structure
  n <- 50
  data <- data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100),
    ilr1 = rnorm(n, 0.5, 0.8),
    ilr2 = rnorm(n, -0.2, 0.6)
  )

  # Assess stationarity
  stationarity_result <- gc_assess_stationarity(data, method = "biplot", plot = FALSE)
  expect_type(stationarity_result$stationary, "logical")
  
  # Assess Gaussianity
  ilr_vals <- data[, c("ilr1", "ilr2")]
  gaussianity_result <- gc_assess_gaussianity(ilr_vals, method = "anderson", plot = FALSE)
  expect_type(gaussianity_result$gaussian, "logical")
  
  # Apply anamorphosis if needed
  sim_data <- matrix(rnorm(100, 0, 1), nrow = 50, ncol = 2)
  colnames(sim_data) <- c("ilr1", "ilr2")
  back_transformed <- gc_apply_anamorphosis(sim_data, as.matrix(ilr_vals))
  expect_equal(dim(back_transformed), dim(sim_data))
})
