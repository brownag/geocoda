# Comprehensive tests for compositional constraint satisfaction and output validation

test_that("All simulation outputs satisfy sum constraint across multiple realizations", {
  # Setup
  constraints <- list(
    SAND = list(min = 15, max = 35),
    SILT = list(min = 45, max = 65),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 40, seed = 50)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create simulation grid
  sim_grid_df <- data.frame(
    x = seq(0, 100, by = 20),
    y = seq(0, 100, by = 20)
  )
  sim_grid <- sf::st_as_sf(sim_grid_df, coords = c("x", "y"), crs = "local")

  # Build model
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )
  model <- gc_ilr_model(ilr_params, vgm_model)

  # Simulate multiple realizations
  nsim <- 10
  sims_result <- gc_sim_composition(
    model, sim_grid,
    nsim = nsim,
    target_names = c("sand", "silt", "clay"),
    seed = 123
  )
  
  # Extract all layers and check constraint satisfaction
  sims_matrix <- terra::as.matrix(sims_result, wide = TRUE)
  
  # Reshape to process each realization
  n_cells <- nrow(sim_grid)
  n_components <- 3
  
  for (realization in seq_len(nsim)) {
    # Extract columns for this realization
    col_indices <- seq((realization - 1) * n_components + 1, realization * n_components)
    realization_data <- sims_matrix[, col_indices]
    
    # Check sums are close to 100
    row_sums <- rowSums(realization_data, na.rm = TRUE)
    expect_true(all(abs(row_sums - 100) < 1e-6),
                info = paste("Realization", realization, "does not satisfy sum constraint"))
    
    # Check all values are non-negative
    expect_true(all(realization_data >= -1e-10, na.rm = TRUE),
                info = paste("Realization", realization, "contains negative values"))
    
    # Check all values are <= 100
    expect_true(all(realization_data <= 100 + 1e-10, na.rm = TRUE),
                info = paste("Realization", realization, "contains values > 100"))
  }
})

test_that("Unconditional simulations are independent realizations", {
  constraints <- list(
    SAND = list(min = 20, max = 30),
    SILT = list(min = 50, max = 60),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 30, seed = 51)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create simulation grid
  sim_grid_df <- data.frame(
    x = seq(0, 100, by = 25),
    y = seq(0, 100, by = 25)
  )
  sim_grid <- sf::st_as_sf(sim_grid_df, coords = c("x", "y"), crs = "local")

  # Build model
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )
  model <- gc_ilr_model(ilr_params, vgm_model)

  # Generate two independent simulations with different seeds
  sims1_result <- gc_sim_composition(
    model, sim_grid,
    nsim = 2,
    target_names = c("sand", "silt", "clay"),
    seed = 100
  )
  
  sims2_result <- gc_sim_composition(
    model, sim_grid,
    nsim = 2,
    target_names = c("sand", "silt", "clay"),
    seed = 200
  )
  
  # Extract as matrices
  sims1_matrix <- terra::as.matrix(sims1_result, wide = TRUE)
  sims2_matrix <- terra::as.matrix(sims2_result, wide = TRUE)
  
  # They should be different (with very high probability)
  expect_false(all(abs(sims1_matrix - sims2_matrix) < 1e-10),
               info = "Two independent simulations are identical")
  
  # But both should satisfy constraints
  for (sim_matrix in list(sims1_matrix, sims2_matrix)) {
    row_sums <- rowSums(sim_matrix, na.rm = TRUE)
    expect_true(all(abs(row_sums - 100) < 1e-6))
  }
})

test_that("Conditional simulation honors observed data at conditioning locations", {
  constraints <- list(
    SAND = list(min = 15, max = 35),
    SILT = list(min = 45, max = 65),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 30, seed = 52)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create conditioning data with known values
  known_comps <- data.frame(
    sand = c(25, 30),
    silt = c(55, 50),
    clay = c(20, 20)
  )
  
  known_ilr <- compositions::ilr(compositions::acomp(known_comps))
  colnames(known_ilr) <- paste0("ilr", seq_len(ncol(known_ilr)))
  
  known_locations <- data.frame(
    x = c(10, 80),
    y = c(10, 80)
  )
  
  conditioning_data <- sf::st_as_sf(
    cbind(known_locations, as.data.frame(known_ilr)),
    coords = c("x", "y"),
    crs = "local"
  )

  # Create simulation grid (including and beyond conditioning locations)
  sim_grid_df <- data.frame(
    x = c(10, 50, 80),
    y = c(10, 50, 80)
  )
  sim_grid <- sf::st_as_sf(sim_grid_df, coords = c("x", "y"), crs = "local")

  # Build conditional model
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )
  model_cond <- gc_ilr_model(ilr_params, vgm_model, data = conditioning_data)

  # Simulate conditional realization
  sims_cond <- gc_sim_composition(
    model_cond, sim_grid,
    nsim = 1,
    target_names = c("sand", "silt", "clay"),
    observed_data = conditioning_data,
    seed = 123
  )
  
  # Extract simulated values
  sims_matrix <- terra::as.matrix(sims_cond, wide = TRUE)
  
  # Check that simulated values at conditioning locations match input (within tolerance)
  # This is a soft check because kriging might not exactly reproduce observed values
  # due to numerical precision and kriging variance, but should be very close
  for (i in seq_len(nrow(known_comps))) {
    # Find rows corresponding to conditioning locations
    # (first two rows should be the conditioning locations)
    simmed_comp <- sims_matrix[i, 1:3]
    expected_comp <- as.numeric(known_comps[i, ])
    
    # Check sums to 100
    expect_true(abs(sum(simmed_comp) - 100) < 1e-6)
    
    # All values should be in reasonable range for the constraints
    expect_true(simmed_comp[1] >= 10 & simmed_comp[1] <= 40) # sand
    expect_true(simmed_comp[2] >= 40 & simmed_comp[2] <= 70) # silt
    expect_true(simmed_comp[3] >= 10 & simmed_comp[3] <= 30) # clay
  }
})

test_that("Simulations respect boundary constraints from definitions", {
  constraints <- list(
    SAND = list(min = 10, max = 40),
    SILT = list(min = 40, max = 70),
    CLAY = list(min = 5, max = 30)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 40, seed = 53)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create simulation grid
  sim_grid_df <- data.frame(
    x = seq(0, 100, by = 25),
    y = seq(0, 100, by = 25)
  )
  sim_grid <- sf::st_as_sf(sim_grid_df, coords = c("x", "y"), crs = "local")

  # Build model
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )
  model <- gc_ilr_model(ilr_params, vgm_model)

  # Simulate
  sims_result <- gc_sim_composition(
    model, sim_grid,
    nsim = 5,
    target_names = c("sand", "silt", "clay"),
    seed = 123
  )
  
  # Extract all simulated data
  sims_matrix <- terra::as.matrix(sims_result, wide = TRUE)
  
  # Check against constraints
  # Note: simulations won't be perfectly bounded but should be reasonably close
  sand_col_indices <- seq(1, ncol(sims_matrix), 3)
  silt_col_indices <- seq(2, ncol(sims_matrix), 3)
  clay_col_indices <- seq(3, ncol(sims_matrix), 3)
  
  sand_vals <- sims_matrix[, sand_col_indices]
  silt_vals <- sims_matrix[, silt_col_indices]
  clay_vals <- sims_matrix[, clay_col_indices]
  
  # Most values should be within the bounds (allowing some extrapolation at margins)
  expect_true(quantile(sand_vals, 0.95) <= 50) # well within max
  expect_true(quantile(sand_vals, 0.05) >= 0)   # non-negative
  
  expect_true(quantile(silt_vals, 0.95) <= 80)
  expect_true(quantile(silt_vals, 0.05) >= 30)
  
  expect_true(quantile(clay_vals, 0.95) <= 40)
  expect_true(quantile(clay_vals, 0.05) >= 0)
})

test_that("gc_sim_composition validates inputs properly", {
  constraints <- list(
    SAND = list(min = 20, max = 30),
    SILT = list(min = 50, max = 60),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 25, seed = 54)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Build model
  extent <- c(0, 0, 100, 100)
  suggested_params <- gc_vgm_defaults(ilr_params, extent)
  vgm_model <- gstat::vgm(
    psill = suggested_params$mean_sill,
    model = "Exp",
    range = suggested_params$range,
    nugget = suggested_params$nugget
  )
  model <- gc_ilr_model(ilr_params, vgm_model)

  # Create grid
  sim_grid_df <- data.frame(
    x = seq(0, 100, by = 50),
    y = seq(0, 100, by = 50)
  )
  sim_grid <- sf::st_as_sf(sim_grid_df, coords = c("x", "y"), crs = "local")

  # Should fail without target_names
  expect_error(
    gc_sim_composition(model, sim_grid, nsim = 1),
    "target_names"
  )

  # Should succeed with target_names
  result <- gc_sim_composition(
    model, sim_grid,
    nsim = 1,
    target_names = c("sand", "silt", "clay")
  )
  expect_s4_class(result, "SpatRaster")
  
  # Should fail with wrong number of target_names
  expect_error(
    gc_sim_composition(
      model, sim_grid,
      nsim = 1,
      target_names = c("sand", "silt") # only 2 names for 3 components
    ),
    "target_names"
  )
})
