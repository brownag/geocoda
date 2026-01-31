library(testthat)
library(sf)
library(terra)

# ============================================================================
# Tests for gc_prepare_zone_data()
# ============================================================================

test_that("gc_prepare_zone_data handles vector zones with spatial join", {
  # Create sample data
  set.seed(42)
  samples_df <- data.frame(
    x = c(10, 20, 30, 40, 50),
    y = c(10, 20, 30, 40, 50),
    sand = c(20, 25, 30, 22, 28),
    silt = c(60, 55, 50, 58, 52),
    clay = c(20, 20, 20, 20, 20)
  )

  # Create zone polygons
  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = c(1, 2),
      geometry = sf::st_sfc(
        sf::st_polygon(list(matrix(c(0, 50, 50, 0, 0, 0, 0, 50, 50, 0), ncol = 2))),
        sf::st_polygon(list(matrix(c(25, 75, 75, 25, 25, 25, 25, 75, 75, 25), ncol = 2)))
      )
    ),
    crs = NA
  )

  # Prepare zone data
  zone_data <- gc_prepare_zone_data(samples_df, zones_sf)

  # Verify output structure
  expect_s3_class(zone_data, "gc_zone_data")
  expect_type(zone_data, "list")
  expect_true(length(zone_data) > 0)

  # Verify zone assignments
  for (z in zone_data) {
    expect_true("zone_id" %in% names(z))
    expect_true("n_samples" %in% names(z))
    expect_true("data" %in% names(z))
    expect_true("valid" %in% names(z))
    expect_true("bounds" %in% names(z))
    expect_true("warning_msgs" %in% names(z))
  }
})

test_that("gc_prepare_zone_data handles zone_id_column parameter", {
  samples_df <- data.frame(
    x = c(10, 20, 30, 40),
    y = c(10, 20, 30, 40),
    zone_id = c(1, 1, 2, 2),
    sand = c(20, 25, 30, 22),
    silt = c(60, 55, 50, 58),
    clay = c(20, 20, 20, 20)
  )

  # Create dummy zones
  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = c(1, 2),
      geometry = sf::st_sfc(
        sf::st_point(c(0, 0)),
        sf::st_point(c(100, 100))
      )
    ),
    crs = NA
  )
  sf::st_geometry(zones_sf) <- sf::st_buffer(sf::st_geometry(zones_sf), 100)

  zone_data <- gc_prepare_zone_data(samples_df, zones_sf, zone_id_column = "zone_id")

  expect_length(zone_data, 2)
  expect_equal(zone_data[["1"]]$n_samples, 2)
  expect_equal(zone_data[["2"]]$n_samples, 2)
})

test_that("gc_prepare_zone_data warns on sparse zones", {
  samples_df <- data.frame(
    x = c(10),
    y = c(10),
    sand = c(20),
    silt = c(60),
    clay = c(20)
  )

  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = 1,
      geometry = sf::st_sfc(
        sf::st_buffer(sf::st_point(c(10, 10)), 100)
      )
    ),
    crs = NA
  )

  zone_data <- gc_prepare_zone_data(samples_df, zones_sf, min_samples = 5)

  expect_false(zone_data[["1"]]$valid)
  expect_true(length(zone_data[["1"]]$warning_msgs) > 0)
})

test_that("gc_prepare_zone_data validates compositional sums", {
  samples_df <- data.frame(
    x = c(10, 20),
    y = c(10, 20),
    sand = c(20, 50),  # Second row: 50+50+50 = 150 (invalid)
    silt = c(60, 50),
    clay = c(20, 50)
  )

  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = 1,
      geometry = sf::st_sfc(
        sf::st_buffer(sf::st_point(c(15, 15)), 100)
      )
    ),
    crs = NA
  )

  zone_data <- gc_prepare_zone_data(samples_df, zones_sf)

  # Should flag compositional sum issue
  expect_true(any(grepl("sums outside", zone_data[["1"]]$warning_msgs)))
})

# ============================================================================
# Tests for gc_fit_zone_models()
# ============================================================================

test_that("gc_fit_zone_models builds models successfully", {
  # Create minimal zone data
  set.seed(42)
  samples_df <- data.frame(
    x = runif(10, 0, 100),
    y = runif(10, 0, 100),
    sand = rnorm(10, 20, 5),
    silt = rnorm(10, 60, 5),
    clay = rnorm(10, 20, 5)
  )

  # Normalize to sum to 100
  samples_df[, c("sand", "silt", "clay")] <- samples_df[, c("sand", "silt", "clay")] / 
    rowSums(samples_df[, c("sand", "silt", "clay")]) * 100

  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = 1,
      geometry = sf::st_sfc(
        sf::st_buffer(sf::st_point(c(50, 50)), 200)
      )
    ),
    crs = NA
  )

  zone_data <- gc_prepare_zone_data(samples_df, zones_sf)
  zone_models <- gc_fit_zone_models(zone_data, verbose = FALSE)

  expect_s3_class(zone_models, "gc_zone_models")
  expect_equal(zone_models[["1"]]$fit_status, "success")
  expect_true("ilr_params" %in% names(zone_models[["1"]]))
  expect_true("variogram" %in% names(zone_models[["1"]]))
  expect_true("model" %in% names(zone_models[["1"]]))
})

test_that("gc_fit_zone_models skips sparse zones when skip_sparse=TRUE", {
  samples_df <- data.frame(
    x = c(10),
    y = c(10),
    sand = c(20),
    silt = c(60),
    clay = c(20)
  )

  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = 1,
      geometry = sf::st_sfc(
        sf::st_buffer(sf::st_point(c(10, 10)), 100)
      )
    ),
    crs = NA
  )

  zone_data <- gc_prepare_zone_data(samples_df, zones_sf, min_samples = 1)
  zone_models <- gc_fit_zone_models(zone_data, skip_sparse = TRUE, verbose = FALSE)

  expect_equal(zone_models[["1"]]$fit_status, "skipped")
})

test_that("gc_fit_zone_models handles multiple zones", {
  set.seed(42)
  samples_df <- data.frame(
    x = c(runif(5, 0, 50), runif(5, 50, 100)),
    y = c(runif(5, 0, 50), runif(5, 50, 100)),
    zone_id = c(rep(1, 5), rep(2, 5)),
    sand = rnorm(10, 20, 5),
    silt = rnorm(10, 60, 5),
    clay = rnorm(10, 20, 5)
  )

  samples_df[, c("sand", "silt", "clay")] <- samples_df[, c("sand", "silt", "clay")] / 
    rowSums(samples_df[, c("sand", "silt", "clay")]) * 100

  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = c(1, 2),
      geometry = sf::st_sfc(
        sf::st_buffer(sf::st_point(c(25, 25)), 100),
        sf::st_buffer(sf::st_point(c(75, 75)), 100)
      )
    ),
    crs = NA
  )

  zone_data <- gc_prepare_zone_data(samples_df, zones_sf, zone_id_column = "zone_id")
  zone_models <- gc_fit_zone_models(zone_data, verbose = FALSE)

  expect_length(zone_models, 2)
  expect_equal(zone_models[["1"]]$fit_status, "success")
  expect_equal(zone_models[["2"]]$fit_status, "success")
})

# ============================================================================
# Tests for gc_validate_tessellation()
# ============================================================================

test_that("gc_validate_tessellation passes valid raster", {
  # Create a simple valid raster
  set.seed(42)
  r <- terra::rast(nrows = 10, ncols = 10, resolution = 1)
  
  # Create three compositional layers (sand, silt, clay per sim)
  # Ensure sums = 100 without negatives
  sand_vals <- runif(ncell(r), 20, 30)
  silt_vals <- runif(ncell(r), 35, 45)
  clay_vals <- 100 - sand_vals - silt_vals
  # Ensure no negatives (should not happen with these ranges)
  clay_vals <- pmax(clay_vals, 1)

  r_sand <- r
  r_sand[] <- sand_vals
  names(r_sand) <- "sand.1.sim1"

  r_silt <- r
  r_silt[] <- silt_vals
  names(r_silt) <- "silt.1.sim1"

  r_clay <- r
  r_clay[] <- clay_vals
  names(r_clay) <- "clay.1.sim1"

  r_stack <- c(r_sand, r_silt, r_clay)

  validation <- gc_validate_tessellation(r_stack)

  expect_true(validation$valid)
  expect_length(validation$issues, 0)
})

test_that("gc_validate_tessellation detects sum violations", {
  set.seed(42)
  r <- terra::rast(nrows = 10, ncols = 10, resolution = 1)

  # Create invalid raster (sums don't equal 100)
  # Make sure at least some cells violate the constraint significantly
  sand_vals <- c(rep(20, 50), rep(50, 50))  # First 50 cells: 20, last 50: 50
  silt_vals <- c(rep(30, 50), rep(25, 50))
  clay_vals <- c(rep(30, 50), rep(15, 50))  # Sums: 80 in first half, invalid in second

  r_sand <- r
  r_sand[] <- sand_vals
  names(r_sand) <- "sand.1.sim1"

  r_silt <- r
  r_silt[] <- silt_vals
  names(r_silt) <- "silt.1.sim1"

  r_clay <- r
  r_clay[] <- clay_vals
  names(r_clay) <- "clay.1.sim1"

  r_stack <- c(r_sand, r_silt, r_clay)

  validation <- gc_validate_tessellation(r_stack, tolerance = 1.0)

  # At least some cells should violate (second half sums to 90)
  expect_false(validation$valid)
  expect_true(length(validation$issues) > 0)
})

test_that("gc_validate_tessellation handles multiple realizations", {
  set.seed(42)
  r <- terra::rast(nrows = 10, ncols = 10, resolution = 1)

  # Create two realizations
  r_list <- list()
  for (sim in 1:2) {
    sand_vals <- runif(ncell(r), 20, 30)
    silt_vals <- runif(ncell(r), 35, 45)
    clay_vals <- 100 - sand_vals - silt_vals
    clay_vals <- pmax(clay_vals, 1)

    r_sand <- r
    r_sand[] <- sand_vals
    names(r_sand) <- paste0("sand.1.sim", sim)

    r_silt <- r
    r_silt[] <- silt_vals
    names(r_silt) <- paste0("silt.1.sim", sim)

    r_clay <- r
    r_clay[] <- clay_vals
    names(r_clay) <- paste0("clay.1.sim", sim)

    r_list <- c(r_list, list(c(r_sand, r_silt, r_clay)))
  }

  r_stack <- do.call(c, r_list)

  validation <- gc_validate_tessellation(r_stack)

  expect_true(validation$valid)
})

# Integration tests (full workflow testing)

test_that("gc_simulate_by_zones completes full workflow", {
  skip("Gap filling workflow improvements pending")
  set.seed(42)

  # Create sample data
  samples_df <- data.frame(
    x = c(runif(5, 0, 50), runif(5, 50, 100)),
    y = c(runif(5, 0, 50), runif(5, 50, 100)),
    sand = rnorm(10, 20, 3),
    silt = rnorm(10, 60, 3),
    clay = rnorm(10, 20, 3)
  )

  samples_df[, c("sand", "silt", "clay")] <- samples_df[, c("sand", "silt", "clay")] / 
    rowSums(samples_df[, c("sand", "silt", "clay")]) * 100

  # Create zones
  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = c(1, 2),
      geometry = sf::st_sfc(
        sf::st_buffer(sf::st_point(c(25, 25)), 75),
        sf::st_buffer(sf::st_point(c(75, 75)), 75)
      )
    ),
    crs = NA
  )

  # Create extent raster
  extent_rast <- terra::rast(xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                             resolution = 20, crs = NA)

  # Run simulation
  result <- gc_simulate_by_zones(
    zones = zones_sf,
    data = samples_df,
    extent = extent_rast,
    resolution = 20,
    nsim = 2,
    verbose = FALSE
  )

  # Verify output
  expect_s4_class(result, "SpatRaster")
  expect_true(terra::nlyr(result) > 0)

  # Verify layer naming
  layer_names <- names(result)
  expect_true(any(grepl("sim1", layer_names)))
  expect_true(any(grepl("sim2", layer_names)))
})

test_that("gc_simulate_by_zones produces valid tessellation", {
  skip("Gap filling workflow improvements pending")
  set.seed(42)

  samples_df <- data.frame(
    x = runif(15, 0, 100),
    y = runif(15, 0, 100),
    sand = rnorm(15, 20, 3),
    silt = rnorm(15, 60, 3),
    clay = rnorm(15, 20, 3)
  )

  samples_df[, c("sand", "silt", "clay")] <- samples_df[, c("sand", "silt", "clay")] / 
    rowSums(samples_df[, c("sand", "silt", "clay")]) * 100

  zones_sf <- sf::st_as_sf(
    data.frame(
      zone_id = 1,
      geometry = sf::st_sfc(
        sf::st_buffer(sf::st_point(c(50, 50)), 100)
      )
    ),
    crs = NA
  )

  extent_rast <- terra::rast(xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                             resolution = 25, crs = NA)

  result <- gc_simulate_by_zones(
    zones = zones_sf,
    data = samples_df,
    extent = extent_rast,
    resolution = 25,
    nsim = 1,
    verbose = FALSE
  )

  # Validate constraints
  validation <- gc_validate_tessellation(result)
  expect_true(validation$valid)
})
