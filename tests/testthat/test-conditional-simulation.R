test_that("gc_ilr_model accepts data for conditioning", {
  # Setup
  constraints <- list(
    SAND = list(min = 10, max = 40),
    SILT = list(min = 40, max = 70),
    CLAY = list(min = 10, max = 30)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 30, seed = 42)
  samples <- samples_res$samples

  # Estimate ILR parameters
  ilr_params <- gc_ilr_params(samples)

  # Create grid for model building
  model_grid <- gc_expand_bounds(constraints, step = 10, target_sum = 100)
  model_res <- gc_resample_compositions(model_grid, method = "uniform", n = 20, seed = 43)
  model_samples <- model_res$samples

  # Transform samples to ILR space for conditioning data
  ilr_data <- compositions::ilr(compositions::acomp(model_samples))
  colnames(ilr_data) <- paste0("ilr", seq_len(ncol(ilr_data)))

  # Create sf object with coordinates
  set.seed(1)
  coords <- data.frame(
    x = runif(nrow(ilr_data), 0, 100),
    y = runif(nrow(ilr_data), 0, 100)
  )
  conditioning_data <- sf::st_as_sf(
    cbind(coords, as.data.frame(ilr_data)),
    coords = c("x", "y"),
    crs = "local"
  )

  # Get variogram
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )

  # Build model WITH conditioning data
  model_conditional <- gc_ilr_model(ilr_params, vgm_model, data = conditioning_data)
  expect_s3_class(model_conditional, "gstat")

  # Build model WITHOUT conditioning data (default)
  model_unconditional <- gc_ilr_model(ilr_params, vgm_model, data = NULL)
  expect_s3_class(model_unconditional, "gstat")
})

test_that("gc_fit_vgm handles multiple ILR dimensions", {
  constraints <- list(
    SAND = list(min = 15, max = 35),
    SILT = list(min = 45, max = 65),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 50, seed = 45)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create sample locations with better spatial structure
  # Use a regular grid with some jitter for realistic spatial correlation
  set.seed(42)
  n_side <- 5 # 5x5 grid = 25 samples
  sample_locs <- expand.grid(
    x = seq(0, 80, length.out = n_side),
    y = seq(0, 80, length.out = n_side)
  )
  # Add small random jitter to locations
  sample_locs$x <- sample_locs$x + rnorm(nrow(sample_locs), 0, 2)
  sample_locs$y <- sample_locs$y + rnorm(nrow(sample_locs), 0, 2)

  # Create ILR data at locations with spatial correlation
  # Bootstrap from the grid to ensure realistic compositions
  n_samples <- nrow(sample_locs)
  samples_subset <- samples[sample(nrow(samples), n_samples, replace = TRUE), ]
  ilr_samples <- compositions::ilr(compositions::acomp(samples_subset))
  colnames(ilr_samples) <- paste0("ilr", seq_len(ncol(ilr_samples)))

  # Combine into data frame with x, y columns (not as sf)
  sample_data_df <- cbind(sample_locs, as.data.frame(ilr_samples))

  # Test: optimize with aggregate = FALSE
  result_per_dim <- gc_fit_vgm(
    ilr_params,
    data = sample_data_df,
    aggregate = FALSE
  )

  expect_type(result_per_dim, "list")
  expect_true(all(c("fitted_vgms", "empirical_vgms", "fitted_params") %in% names(result_per_dim)))
  expect_type(result_per_dim$fitted_vgms, "list")
  expect_equal(length(result_per_dim$fitted_vgms), 2) # 3 components = 2 ILR

  # Test: optimize with aggregate = TRUE
  result_agg <- gc_fit_vgm(
    ilr_params,
    data = sample_data_df,
    aggregate = TRUE
  )

  expect_s3_class(result_agg, c("variogramModel", "data.frame"))
  expect_true("psill" %in% colnames(result_agg))
  expect_true("range" %in% colnames(result_agg))
})

test_that("Empirical variogram fitting captures dimension-specific structure", {
  constraints <- list(
    SAND = list(min = 20, max = 30),
    SILT = list(min = 50, max = 60),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 25, seed = 46)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create spatial dataset with stronger correlative structure
  set.seed(123)
  n_side <- 6 # 6x6 grid = 36 samples for better variogram fitting
  coords <- expand.grid(
    x = seq(0, 80, length.out = n_side),
    y = seq(0, 80, length.out = n_side)
  )
  # Add small jitter for realistic variation
  coords$x <- coords$x + rnorm(nrow(coords), 0, 1.5)
  coords$y <- coords$y + rnorm(nrow(coords), 0, 1.5)

  # Create spatially correlated ILR values with clear structure
  ilr1_vals <- 0.5 * sin((coords$x + coords$y) / 20) + rnorm(nrow(coords), sd = 0.15)
  ilr2_vals <- 0.3 * cos(coords$x / 15) + rnorm(nrow(coords), sd = 0.2)
  ilr3_vals <- 0.1 * (coords$x + coords$y) / 100 + rnorm(nrow(coords), sd = 0.15)

  # Create as data frame with x, y columns (not sf)
  sample_data_df <- data.frame(
    x = coords$x,
    y = coords$y,
    ilr1 = ilr1_vals,
    ilr2 = ilr2_vals,
    ilr3 = ilr3_vals
  )

  # Fit variograms
  result <- gc_fit_vgm(
    ilr_params,
    data = sample_data_df,
    aggregate = FALSE
  )

  expect_equal(length(result$empirical_vgms), 2)

  # Extract fitted parameters
  fitted_params <- result$fitted_params
  expect_s3_class(fitted_params, "data.frame")
  expect_true(all(c("ilr_id", "range", "nugget", "psill") %in% colnames(fitted_params)))
  expect_equal(nrow(fitted_params), 2)

  expect_true(all(fitted_params$range > 0))
  expect_true(all(fitted_params$psill > 0))
  expect_true(all(fitted_params$nugget >= 0))
})

test_that("Function signatures are correct", {
  constraints <- list(
    SAND = list(min = 15, max = 35),
    SILT = list(min = 45, max = 65),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 20, seed = 47)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  sim_grid_df <- data.frame(
    x = seq(0, 100, by = 25),
    y = seq(0, 100, by = 25)
  )
  sim_grid <- sf::st_as_sf(sim_grid_df, coords = c("x", "y"), crs = "local")

  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )

  # Call gc_ilr_model with data parameter
  model <- gc_ilr_model(ilr_params, vgm_model, data = NULL)
  expect_s3_class(model, "gstat")

  # Call gc_sim_composition with required target_names
  result <- gc_sim_composition(
    model, sim_grid,
    nsim = 1,
    target_names = c("sand", "silt", "clay")
  )
  expect_s4_class(result, "SpatRaster")
})
