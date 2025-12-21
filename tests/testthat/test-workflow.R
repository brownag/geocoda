test_that("gc_expand_bounds generates valid grid", {
  constraints <- list(
    SAND = list(min = 0, max = 40),
    SILT = list(min = 50, max = 80),
    CLAY = list(min = 10, max = 20)
  )

  grid <- gc_expand_bounds(constraints, step = 1.0, target_sum = 100)

  # Check that grid has correct column names
  expect_equal(colnames(grid), c("SAND", "SILT", "CLAY"))

  # Check that all rows sum to 100 (within tolerance)
  row_sums <- rowSums(grid)
  expect_true(all(abs(row_sums - 100) < 1e-6))

  # Check that all values are within bounds
  expect_true(all(grid$SAND >= 0 & grid$SAND <= 40))
  expect_true(all(grid$SILT >= 50 & grid$SILT <= 80))
  expect_true(all(grid$CLAY >= 10 & grid$CLAY <= 20))

  # Check that grid is non-empty
  expect_gt(nrow(grid), 0)
})

test_that("gc_expand_bounds rejects impossible constraints", {
  # Min values sum to more than target
  bad_constraints <- list(
    SAND = list(min = 40, max = 60),
    SILT = list(min = 40, max = 60),
    CLAY = list(min = 30, max = 40)
  )

  expect_error(
    gc_expand_bounds(bad_constraints, target_sum = 100),
    "Minimum values sum"
  )
})

test_that("gc_resample_compositions generates correct number of samples", {
  constraints <- list(
    SAND = list(min = 0, max = 40),
    SILT = list(min = 50, max = 80),
    CLAY = list(min = 10, max = 20)
  )

  grid <- gc_expand_bounds(constraints, step = 1.0, target_sum = 100)

  result <- gc_resample_compositions(grid, n = 100, method = "uniform", seed = 42)

  expect_type(result, "list")
  expect_named(result, c("samples", "method", "n"))
  expect_equal(nrow(result$samples), 100)
  expect_equal(result$method, "uniform")
  expect_equal(result$n, 100)
})

test_that("gc_ilr_params produces correct output structure", {
  samples <- data.frame(
    sand = c(20, 25, 30, 22),
    silt = c(60, 55, 50, 58),
    clay = c(20, 20, 20, 20)
  )

  params <- gc_ilr_params(samples)

  expect_type(params, "list")
  expect_named(params, c("mean", "cov", "names", "base_class"))

  # ILR mean should have length D-1 = 2
  expect_length(params$mean, 2)

  # Covariance should be 2x2
  expect_equal(dim(params$cov), c(2, 2))

  # Names should match input
  expect_equal(params$names, c("sand", "silt", "clay"))

  # Base class should be "acomp"
  expect_equal(params$base_class, "acomp")
})

test_that("gc_ilr_model creates valid gstat object", {
  samples <- data.frame(
    sand = c(20, 25, 30, 22, 18),
    silt = c(60, 55, 50, 58, 62),
    clay = c(20, 20, 20, 20, 20)
  )

  params <- gc_ilr_params(samples)

  library(gstat)
  vgm_template <- vgm(psill = 1, model = "Exp", range = 30, nugget = 0.01)

  model <- gc_ilr_model(params, variogram_model = vgm_template)

  # Should be a gstat object
  expect_s3_class(model, "gstat")
})

test_that("Full workflow: expand, bootstrap, estimate, build, simulate", {
  # Define constraints
  constraints <- list(
    SAND = list(min = 0, max = 40),
    SILT = list(min = 50, max = 80),
    CLAY = list(min = 10, max = 20)
  )

  # Expand grid
  grid <- gc_expand_bounds(constraints, step = 2.0, target_sum = 100)
  expect_gt(nrow(grid), 0)

  # Bootstrap samples
  set.seed(123)
  boot_result <- gc_resample_compositions(grid, n = 500, method = "uniform")
  expect_equal(nrow(boot_result$samples), 500)

  # Estimate parameters
  params <- gc_ilr_params(boot_result$samples)
  expect_length(params$mean, 2)

  # Build model
  library(gstat)
  vgm_template <- vgm(psill = 1, model = "Exp", range = 30, nugget = 0.01)
  model <- gc_ilr_model(params, variogram_model = vgm_template)
  expect_s3_class(model, "gstat")

  # Create spatial grid as sf object
  x.range <- seq(0, 100, by = 10)
  y.range <- seq(0, 100, by = 10)
  grid_df <- expand.grid(x = x.range, y = y.range)
  grid_sf <- sf::st_as_sf(grid_df, coords = c("x", "y"), crs = "local")

  # Simulate
  sims <- gc_sim_composition(
    model,
    grid_sf,
    nsim = 1,
    target_names = c("sand", "silt", "clay")
  )

  # Check output structure
  expect_s4_class(sims, "SpatRaster")
  expect_equal(terra::nlyr(sims), 3)

  # Check that layer names are correct
  expected_names <- c("sand.sim1", "silt.sim1", "clay.sim1")
  expect_equal(names(sims), expected_names)

  # Extract values and check composition
  vals <- as.data.frame(terra::values(sims))
  row_sums <- rowSums(vals)

  # Check sums are close to 100
  expect_true(all(abs(row_sums - 100) < 0.01))

  # Check all values are non-negative
  expect_true(all(as.matrix(vals) >= -1e-6))

  # Check all values are <= 100
  expect_true(all(as.matrix(vals) <= 100 + 1e-6))
})

test_that("Multiple realizations are properly stacked", {
  # Smaller example for speed
  constraints <- list(
    SAND = list(min = 10, max = 30),
    SILT = list(min = 50, max = 70),
    CLAY = list(min = 10, max = 30)
  )

  grid <- gc_expand_bounds(constraints, step = 5.0, target_sum = 100)
  set.seed(456)
  boot_result <- gc_resample_compositions(grid, n = 200, method = "uniform")
  params <- gc_ilr_params(boot_result$samples)

  library(gstat)
  vgm_template <- vgm(psill = 1, model = "Exp", range = 30, nugget = 0.01)
  model <- gc_ilr_model(params, variogram_model = vgm_template)

  # Small spatial grid as sf object
  x.range <- seq(0, 50, by = 10)
  y.range <- seq(0, 50, by = 10)
  grid_df <- expand.grid(x = x.range, y = y.range)
  grid_sf <- sf::st_as_sf(grid_df, coords = c("x", "y"), crs = "local")

  # Simulate multiple realizations
  sims <- gc_sim_composition(
    model,
    grid_sf,
    nsim = 3,
    target_names = c("sand", "silt", "clay")
  )

  # Should have 3 components * 3 realizations = 9 layers
  expect_equal(terra::nlyr(sims), 9)

  # Check naming pattern: comp.sim<N>
  expected_names <- c(
    "sand.sim1", "silt.sim1", "clay.sim1",
    "sand.sim2", "silt.sim2", "clay.sim2",
    "sand.sim3", "silt.sim3", "clay.sim3"
  )
  expect_equal(names(sims), expected_names)

  # Check sums for each realization
  vals <- as.data.frame(terra::values(sims))

  # Sum by realization
  for (sim_idx in 1:3) {
    col_indices <- grep(paste0("sim", sim_idx), colnames(vals))
    sim_cols <- vals[, col_indices]
    sums <- rowSums(sim_cols)
    expect_true(all(abs(sums - 100) < 0.01))
  }
})

test_that("gc_vgm_defaults returns reasonable values", {
  samples <- data.frame(
    sand = c(20, 25, 30, 22),
    silt = c(60, 55, 50, 58),
    clay = c(20, 20, 20, 20)
  )

  params <- gc_ilr_params(samples)
  extent <- c(0, 0, 100, 100)

  suggestions <- gc_vgm_defaults(params, extent)

  expect_type(suggestions, "list")
  expect_named(suggestions, c("range", "nugget", "mean_sill"))

  # Range should be positive
  expect_gt(suggestions$range, 0)

  # Nugget should be positive
  expect_gt(suggestions$nugget, 0)

  # Mean sill should be positive
  expect_gt(suggestions$mean_sill, 0)

  # Nugget should be smaller than sill
  expect_lt(suggestions$nugget, suggestions$mean_sill)
})
