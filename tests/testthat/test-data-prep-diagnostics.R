test_that("gc_handle_zeros processes data without errors", {
  # Create test data with some zeros
  comp_data <- data.frame(
    sand = c(40, 35, 0, 42, 38),
    silt = c(35, 40, 45, 36, 40),
    clay = c(25, 25, 55, 22, 22)
  )

  # Test with default method (mzero)
  result <- gc_handle_zeros(comp_data)

  expect_type(result, "list")
  expect_true(all(c("imputed_data", "n_zeros_imputed", "imputation_rate", 
                     "method_used", "row_status") %in% names(result)))
  expect_equal(nrow(result$imputed_data), nrow(comp_data))
  expect_equal(ncol(result$imputed_data), ncol(comp_data))
  expect_true(result$n_zeros_imputed > 0)
  expect_true(result$imputation_rate > 0)
})

test_that("gc_handle_zeros identifies problematic rows", {
  comp_data <- data.frame(
    sand = c(40, 35, 0, 42),
    silt = c(35, 40, 45, 36),
    clay = c(25, 25, 55, 22)
  )

  result <- gc_handle_zeros(comp_data, method = "mzero")
  
  expect_equal(length(result$row_status), nrow(comp_data))
  expect_equal(result$row_status[1], factor("observed", levels = c("observed", "imputed", "failed")))
  # Row 3 with zero may be marked as failed if imputation fails (depends on zCompositions version)
  expect_true(result$row_status[3] %in% c(factor("imputed", levels = c("observed", "imputed", "failed")),
                                           factor("failed", levels = c("observed", "imputed", "failed"))))
})

test_that("gc_handle_zeros validates input", {
  # Empty data frame
  expect_error(gc_handle_zeros(data.frame()), "must have at least one row")

  # Single column
  expect_error(gc_handle_zeros(data.frame(x = c(1, 2))), "at least 2 columns")

  # Non-data frame input
  expect_error(gc_handle_zeros(matrix(1:4, nrow = 2)), "must be a data frame")
})

test_that("gc_fit_vgm with correct.diagonal parameter", {
  constraints <- list(
    SAND = list(min = 15, max = 35),
    SILT = list(min = 45, max = 65),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 50, seed = 45)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create sample locations
  set.seed(42)
  n_side <- 5
  sample_locs <- expand.grid(
    x = seq(0, 80, length.out = n_side),
    y = seq(0, 80, length.out = n_side)
  )
  sample_locs$x <- sample_locs$x + rnorm(nrow(sample_locs), 0, 2)
  sample_locs$y <- sample_locs$y + rnorm(nrow(sample_locs), 0, 2)

  n_samples <- nrow(sample_locs)
  samples_subset <- samples[sample(nrow(samples), n_samples, replace = TRUE), ]
  ilr_samples <- compositions::ilr(compositions::acomp(samples_subset))
  colnames(ilr_samples) <- paste0("ilr", seq_len(ncol(ilr_samples)))

  sample_data_df <- cbind(sample_locs, as.data.frame(ilr_samples))

  # Test with default correct.diagonal
  result <- gc_fit_vgm(ilr_params, sample_data_df, aggregate = TRUE)
  
  expect_s3_class(result, "variogramModel")
  expect_true(!is.null(attr(result, "lmc_admissibility")))
  
  # Test with custom correct.diagonal
  result2 <- gc_fit_vgm(ilr_params, sample_data_df, aggregate = TRUE, 
                        correct.diagonal = 1.05)
  
  expect_s3_class(result2, "variogramModel")
  # Check that the adjusted sill is stored in fitted_params
  fitted_params <- attr(result2, "fitted_params")
  expect_true("psill_corrected" %in% colnames(fitted_params))
})

test_that("gc_fit_vgm validates correct.diagonal parameter", {
  comp_data <- data.frame(
    sand = c(40, 35, 30, 42),
    silt = c(35, 40, 45, 36),
    clay = c(25, 25, 25, 22)
  )
  ilr_params <- gc_ilr_params(comp_data)
  
  sample_data <- data.frame(
    x = c(0, 10, 20, 30),
    y = c(0, 10, 20, 30),
    ilr1 = c(0.1, -0.1, 0.0, 0.2),
    ilr2 = c(0.05, 0.1, -0.05, 0.0)
  )

  # Invalid correct.diagonal (less than 1.0)
  expect_error(gc_fit_vgm(ilr_params, sample_data, correct.diagonal = 0.9),
               "correct.diagonal must be numeric and >= 1.0")
})

test_that("gc_validate_conditioning returns expected structure", {
  constraints <- list(
    SAND = list(min = 10, max = 40),
    SILT = list(min = 40, max = 70),
    CLAY = list(min = 10, max = 30)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 30, seed = 42)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create conditioning data
  set.seed(1)
  coords <- data.frame(
    x = runif(10, 0, 100),
    y = runif(10, 0, 100)
  )
  ilr_data <- compositions::ilr(compositions::acomp(samples[1:10, ]))
  colnames(ilr_data) <- paste0("ilr", seq_len(ncol(ilr_data)))
  
  # Create data frame (not sf) to avoid geometry issues in test
  conditioning_data_df <- data.frame(coords, as.data.frame(ilr_data))

  # Build model
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )
  
  # Convert to sf for model building
  conditioning_data_sf <- sf::st_as_sf(
    conditioning_data_df,
    coords = c("x", "y"),
    crs = "local"
  )
  
  model <- gc_ilr_model(ilr_params, vgm_model, data = conditioning_data_sf)

  # Validate conditioning using data frame input
  validation <- gc_validate_conditioning(model, conditioning_data_df)

  expect_type(validation, "list")
  expect_true(all(c("predictions_at_obs", "observed_values", "residuals",
                     "error_metrics", "overall_metrics") %in% names(validation)))
  expect_equal(nrow(validation$residuals), nrow(conditioning_data_df))
  expect_equal(nrow(validation$error_metrics), 2) # 2 ILR dimensions
})

test_that("gc_validate_conditioning computes correct error metrics", {
  # Create simple test case with known compositions
  comp_data <- data.frame(
    SAND = c(40, 35, 30),
    SILT = c(35, 40, 45),
    CLAY = c(25, 25, 25)
  )
  
  conditioning_data <- data.frame(
    x = c(0, 50, 100),
    y = c(0, 50, 100)
  )

  # Transform to ILR
  ilr_data <- compositions::ilr(compositions::acomp(comp_data))
  colnames(ilr_data) <- c("ilr1", "ilr2")
  conditioning_data <- cbind(conditioning_data, as.data.frame(ilr_data))

  ilr_params <- gc_ilr_params(comp_data)
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )
  
  # Convert to sf for model building
  conditioning_data_sf <- sf::st_as_sf(
    conditioning_data,
    coords = c("x", "y"),
    crs = "local"
  )
  
  model <- gc_ilr_model(ilr_params, vgm_model, data = conditioning_data_sf)
  validation <- gc_validate_conditioning(model, conditioning_data)

  # Check that error metrics are computed
  expect_equal(nrow(validation$error_metrics), 2)
  expect_true("RMSE" %in% colnames(validation$error_metrics))
})

test_that("gc_sim_composition accepts nmax parameter", {
  constraints <- list(
    SAND = list(min = 10, max = 40),
    SILT = list(min = 40, max = 70),
    CLAY = list(min = 10, max = 30)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 20, seed = 42)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create variogram model
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )

  # Build model without conditioning data
  model <- gc_ilr_model(ilr_params, vgm_model)

  # Create simulation grid
  grid_df <- expand.grid(x = seq(0, 100, by = 20), y = seq(0, 100, by = 20))

  # Test simulation with nmax
  result <- gc_sim_composition(model, grid_df, nsim = 1, 
                               target_names = c("sand", "silt", "clay"),
                               nmax = 5)

  expect_s4_class(result, "SpatRaster")
  expect_equal(length(names(result)), 3) # 3 components × 1 simulation
  expect_true(all(grepl("sand|silt|clay", names(result))))
})

test_that("gc_sim_composition validates nmax parameter", {
  comp_data <- data.frame(
    sand = c(40, 35, 30),
    silt = c(35, 40, 45),
    clay = c(25, 25, 25)
  )
  ilr_params <- gc_ilr_params(comp_data)
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )
  model <- gc_ilr_model(ilr_params, vgm_model)

  grid_df <- expand.grid(x = seq(0, 100, by = 50), y = seq(0, 100, by = 50))

  # Invalid nmax (not integer)
  expect_error(gc_sim_composition(model, grid_df, nmax = 2.5),
               "nmax must be a positive integer or NULL")

  # Invalid nmax (negative)
  expect_error(gc_sim_composition(model, grid_df, nmax = -1),
               "nmax must be a positive integer or NULL")

  # Valid nmax should work
  result <- gc_sim_composition(model, grid_df, nmax = 3)
  expect_s4_class(result, "SpatRaster")
})

test_that("gc_sim_composition works without nmax (backward compatibility)", {
  constraints <- list(
    SAND = list(min = 10, max = 40),
    SILT = list(min = 40, max = 70),
    CLAY = list(min = 10, max = 30)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 20, seed = 42)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )

  model <- gc_ilr_model(ilr_params, vgm_model)
  grid_df <- expand.grid(x = seq(0, 100, by = 50), y = seq(0, 100, by = 50))

  # Call without nmax parameter (should use all data)
  result <- gc_sim_composition(model, grid_df, nsim = 1,
                               target_names = c("sand", "silt", "clay"))

  expect_s4_class(result, "SpatRaster")
  expect_equal(length(names(result)), 3)
})
