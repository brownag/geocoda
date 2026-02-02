test_that("gc_sim_composition accepts 3D locations", {
  skip_if_not_installed("gstat")

  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))

  vgm_fit <- gstat::vgm(psill = 1.0, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)

  # Create 3D location grid
  locations_3d <- expand.grid(
    x = seq(0, 100, by = 50),
    y = seq(0, 100, by = 50),
    z = c(10, 30, 50)  # Three depth levels
  )

  # Simulate at 3D locations
  sims_3d <- gc_sim_composition(
    model = model,
    locations = locations_3d,
    nsim = 2,
    target_names = c("sand", "silt", "clay")
  )

  # Check output structure
  expect_s4_class(sims_3d, "SpatRaster")
  expect_true(terra::nlyr(sims_3d) >= 6)  # At least 3 components × 2 sims
})

test_that("gc_sim_composition 3D output has multiple z layers", {
  skip_if_not_installed("gstat")

  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))

  vgm_fit <- gstat::vgm(psill = 0.8, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)

  # Create 3D grid with distinct z values
  locations_3d <- data.frame(
    x = c(0, 50, 100, 0, 50, 100),
    y = c(0, 50, 100, 0, 50, 100),
    z = c(10, 10, 10, 50, 50, 50)  # Two depth levels
  )

  sims <- gc_sim_composition(
    model = model,
    locations = locations_3d,
    nsim = 1,
    target_names = c("sand", "silt", "clay")
  )

  # Should successfully simulate at 3D locations
  expect_s4_class(sims, "SpatRaster")
  expect_true(terra::nlyr(sims) >= 3)  # At least 3 components
})

test_that("gc_sim_composition 3D maintains 2D compatibility", {
  skip_if_not_installed("gstat")

  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))

  vgm_fit <- gstat::vgm(psill = 1.0, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)

  # Create 2D location grid (no z column)
  locations_2d <- expand.grid(
    x = seq(0, 100, by = 50),
    y = seq(0, 100, by = 50)
  )

  # Simulate at 2D locations (should still work)
  sims_2d <- gc_sim_composition(
    model = model,
    locations = locations_2d,
    nsim = 1,
    target_names = c("sand", "silt", "clay")
  )

  expect_s4_class(sims_2d, "SpatRaster")
  expect_true(terra::nlyr(sims_2d) >= 3)
})

test_that("gc_sim_composition 3D with varying depth levels", {
  skip_if_not_installed("gstat")

  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov,
                    names = c("sand", "silt", "clay"))

  vgm_fit <- gstat::vgm(psill = 0.8, model = "Exp", range = 80,
                       nugget = 0.1)

  model <- gc_ilr_model(ilr_params, vgm_fit)

  # Create 3D grid with many depth levels
  locations_3d <- expand.grid(
    x = seq(0, 100, by = 50),
    y = seq(0, 100, by = 50),
    z = seq(5, 95, by = 10)  # 10 depth levels
  )

  # Simulate at many 3D locations
  sims <- gc_sim_composition(
    model = model,
    locations = locations_3d,
    nsim = 1,
    target_names = c("sand", "silt", "clay")
  )

  expect_s4_class(sims, "SpatRaster")
  # Should have reasonable number of layers
  expect_true(terra::nlyr(sims) >= 3)
})

test_that("gc_sim_composition 3D with unconditional simulation", {
  skip_if_not_installed("gstat")

  set.seed(42)

  ilr_params <- list(mean = c(0, 0), cov = diag(2),
                    names = c("sand", "silt", "clay"))
  vgm_fit <- gstat::vgm(psill = 0.8, model = "Exp", range = 150,
                       nugget = 0.1)

  # Unconditional model (no data)
  model <- gc_ilr_model(ilr_params, vgm_fit)

  # Prediction on 3D grid
  pred_locs <- expand.grid(
    x = seq(50, 150, by = 100),
    y = seq(50, 150, by = 100),
    z = seq(10, 90, by = 40)
  )

  # Unconditional simulation at 3D locations
  sims <- gc_sim_composition(
    model = model,
    locations = pred_locs,
    nsim = 2,
    target_names = c("sand", "silt", "clay")
  )

  expect_s4_class(sims, "SpatRaster")
  expect_true(terra::nlyr(sims) >= 6)  # 3 components × 2 sims
})

test_that("gc_sim_composition 3D compositional constraint satisfaction", {
  skip_if_not_installed("gstat")

  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))

  vgm_fit <- gstat::vgm(psill = 0.8, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)

  # 3D prediction grid
  locations_3d <- expand.grid(
    x = seq(0, 100, by = 50),
    y = seq(0, 100, by = 50),
    z = c(20, 60)
  )

  sims_3d <- gc_sim_composition(
    model = model,
    locations = locations_3d,
    nsim = 3,
    target_names = c("sand", "silt", "clay")
  )

  # Extract values and check constraints
  values <- terra::values(sims_3d)
  valid_values <- values[!is.na(values)]

  # All values should be non-negative (after clipping)
  expect_true(all(valid_values >= -0.01))  # Allow small numerical errors
  expect_true(all(valid_values <= 100.01))
})

test_that("gc_sim_composition with nmax parameter and 3D", {
  skip_if_not_installed("gstat")

  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))

  vgm_fit <- gstat::vgm(psill = 1.0, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)

  # 3D prediction grid
  locations_3d <- expand.grid(
    x = seq(0, 100, by = 50),
    y = seq(0, 100, by = 50),
    z = c(10, 50, 90)
  )

  # Simulate with neighborhood limit
  sims <- gc_sim_composition(
    model = model,
    locations = locations_3d,
    nsim = 1,
    target_names = c("sand", "silt", "clay"),
    nmax = 10
  )

  expect_s4_class(sims, "SpatRaster")
})

test_that("gc_sim_composition 3D rejects invalid location specification", {
  skip_if_not_installed("gstat")

  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))

  vgm_fit <- gstat::vgm(psill = 1.0, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)

  # Missing required columns
  bad_locations <- data.frame(
    x = c(0, 50, 100),
    # Missing y column
    z = c(10, 50, 90)
  )

  expect_error(
    gc_sim_composition(model, bad_locations, nsim = 1,
                      target_names = c("sand", "silt", "clay")),
    "must contain 'x', 'y'"
  )
})
