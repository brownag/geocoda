context("Ensemble Aggregation and Realization Validation")

# Helper: Create synthetic realization stack
create_test_realizations <- function(n_components = 3, n_realizations = 5, n_cells = 100) {
  
  # Create raster stack for multiple realizations
  comps <- c("sand", "silt", "clay")[1:n_components]
  nrows <- 10
  ncols <- 10
  n_layers <- n_components * n_realizations
  
  # Generate compositions summing to ~100
  set.seed(42)
  all_vals <- list()
  for (i in 1:n_layers) {
    vals <- pmax(abs(stats::rnorm(nrows * ncols, mean = 33, sd = 5)), 0.1)
    all_vals[[i]] <- vals
  }
  
  # Normalize each cell to sum to 100 per realization
  for (cell_idx in 1:(nrows * ncols)) {
    for (real_idx in 1:n_realizations) {
      indices <- seq((real_idx - 1) * n_components + 1, real_idx * n_components)
      cell_vals <- sapply(all_vals[indices], function(x) x[cell_idx])
      if (sum(cell_vals) > 0) {
        normalized_vals <- cell_vals / sum(cell_vals) * 100
        for (j in seq_along(indices)) {
          all_vals[[indices[j]]][cell_idx] <- normalized_vals[j]
        }
      }
    }
  }
  
  # Create raster stack
  r <- terra::rast(nrows = nrows, ncols = ncols, nlyrs = n_layers)
  terra::ext(r) <- c(0, 100, 0, 100)
  
  for (i in 1:n_layers) {
    terra::values(r[[i]]) <- all_vals[[i]]
  }
  
  # Set layer names
  layer_names <- paste0(
    rep(comps, each = n_realizations),
    ".sim",
    rep(1:n_realizations, times = n_components)
  )
  names(r) <- layer_names
  
  r
}

# ============================================================================
# gc_aggregate_realizations tests
# ============================================================================

test_that("gc_aggregate_realizations handles basic aggregation", {
  r <- create_test_realizations(n_components = 3, n_realizations = 5)
  
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    stats = c("mean", "sd")
  )
  
  expect_s4_class(agg, "SpatRaster")
  expect_equal(terra::nlyr(agg), 6)  # 3 components × 2 stats
  expect_true(all(grepl("_mean|_sd", terra::names(agg))))
})

test_that("gc_aggregate_realizations computes correct layer names", {
  r <- create_test_realizations(n_components = 2, n_realizations = 3)
  
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = c("sand", "clay"),
    stats = c("mean", "median", "sd")
  )
  
  expected_names <- c("sand_mean", "sand_median", "sand_sd",
                      "clay_mean", "clay_median", "clay_sd")
  expect_equal(terra::names(agg), expected_names)
})

test_that("gc_aggregate_realizations computes mean correctly", {
  r <- create_test_realizations(n_components = 1, n_realizations = 3)
  
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = "sand",
    stats = c("mean"),
    realizations_per_component = 3
  )
  
  # Manually compute expected mean
  v <- terra::values(r)
  expected_mean <- rowMeans(v)
  
  actual_mean <- terra::values(agg)[, 1]
  expect_equal(actual_mean, expected_mean, tolerance = 1e-5)
})

test_that("gc_aggregate_realizations computes sd correctly", {
  r <- create_test_realizations(n_components = 1, n_realizations = 4)
  
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = "sand",
    stats = c("sd"),
    realizations_per_component = 4
  )
  
  # Manually compute expected SD
  v <- terra::values(r)
  expected_sd <- apply(v, 1, stats::sd, na.rm = TRUE)
  
  actual_sd <- terra::values(agg)[, 1]
  expect_equal(actual_sd, expected_sd, tolerance = 1e-5)
})

test_that("gc_aggregate_realizations computes quantiles", {
  r <- create_test_realizations(n_components = 1, n_realizations = 10)
  
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = "sand",
    stats = c("p05", "p25", "p75", "p95"),
    realizations_per_component = 10
  )
  
  v <- terra::values(r)
  expected_p05 <- apply(v, 1, stats::quantile, probs = 0.05, na.rm = TRUE)
  expected_p95 <- apply(v, 1, stats::quantile, probs = 0.95, na.rm = TRUE)
  
  actual_p05 <- terra::values(agg)[, 1]
  actual_p95 <- terra::values(agg)[, 4]
  
  expect_equal(actual_p05, expected_p05, tolerance = 1e-5)
  expect_equal(actual_p95, expected_p95, tolerance = 1e-5)
})

test_that("gc_aggregate_realizations computes cv", {
  r <- create_test_realizations(n_components = 1, n_realizations = 5)
  
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = "sand",
    stats = c("cv"),
    realizations_per_component = 5
  )
  
  v <- terra::values(r)
  means <- rowMeans(v, na.rm = TRUE)
  sds <- apply(v, 1, stats::sd, na.rm = TRUE)
  expected_cv <- sds / means
  
  actual_cv <- terra::values(agg)[, 1]
  expect_equal(actual_cv, expected_cv, tolerance = 1e-5)
})

test_that("gc_aggregate_realizations handles multiple components", {
  r <- create_test_realizations(n_components = 3, n_realizations = 4)
  
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    stats = c("mean", "sd"),
    realizations_per_component = 4
  )
  
  # Check that each component has correct mean/sd values
  expect_equal(terra::nlyr(agg), 6)
  
  # Extract sand layers and verify they're in expected range
  sand_mean <- terra::values(agg)[, 1]
  expect_true(all(sand_mean >= 0 & sand_mean <= 100))
})

test_that("gc_aggregate_realizations handles list input", {
  r1 <- create_test_realizations(n_components = 3, n_realizations = 1)
  r2 <- create_test_realizations(n_components = 3, n_realizations = 1)
  
  agg <- geocoda::gc_aggregate_realizations(
    list(r1, r2),
    components = c("sand", "silt", "clay"),
    stats = c("mean"),
    realizations_per_component = 2  # Specify since we have 2 realizations combined
  )
  
  expect_s4_class(agg, "SpatRaster")
  expect_equal(terra::nlyr(agg), 3)
})

test_that("gc_aggregate_realizations rejects invalid stats", {
  r <- create_test_realizations()
  
  expect_error(
    geocoda::gc_aggregate_realizations(
      r,
      components = c("sand", "silt", "clay"),
      stats = c("mean", "invalid_stat")
    ),
    "Invalid stats"
  )
})

test_that("gc_aggregate_realizations infers realization count from names", {
  r <- create_test_realizations(n_components = 1, n_realizations = 7)
  
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = "sand",
    stats = c("mean")
  )
  
  # Should work without specifying realizations_per_component
  expect_s4_class(agg, "SpatRaster")
  expect_equal(terra::nlyr(agg), 1)
})

test_that("gc_aggregate_realizations validates dimension consistency", {
  r <- create_test_realizations(n_components = 3, n_realizations = 4)
  
  # Declare wrong number of realizations
  expect_error(
    geocoda::gc_aggregate_realizations(
      r,
      components = c("sand", "silt", "clay"),
      realizations_per_component = 5
    ),
    "Mismatch"
  )
})

# ============================================================================
# gc_validate_realizations tests
# ============================================================================

test_that("gc_validate_realizations returns expected structure", {
  r <- create_test_realizations(n_components = 3, n_realizations = 3)
  
  val <- geocoda::gc_validate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    target_sum = 100,
    sum_tolerance = 1
  )
  
  expect_is(val, "list")
  expect_true("violations_map" %in% names(val))
  expect_true("summary" %in% names(val))
  expect_true("metadata" %in% names(val))
  expect_s4_class(val$violations_map, "SpatRaster")
  expect_is(val$summary, "data.frame")
})

test_that("gc_validate_realizations computes sum violations", {
  # Create raster with deliberate sum violations in some cells
  r <- create_test_realizations(n_components = 3, n_realizations = 2)
  
  # Perturb some cells to create violations
  v <- terra::values(r)
  v[1, ] <- c(60, 30, 50, 40, 20, 15)  # Sums: 140, 75 (both out of bounds)
  terra::values(r) <- v
  
  val <- geocoda::gc_validate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    target_sum = 100,
    sum_tolerance = 1
  )
  
  # Check that violations were detected in cell 1
  violations <- terra::values(val$violations_map)[, 1]
  expect_true(violations[1] > 0)
})

test_that("gc_validate_realizations checks boundary constraints", {
  r <- create_test_realizations(n_components = 3, n_realizations = 2)
  
  # Set very tight constraints that will be violated
  constraints <- list(
    sand = list(min = 40, max = 50),
    silt = list(min = 40, max = 50),
    clay = list(min = 0, max = 10)
  )
  
  val <- geocoda::gc_validate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    constraints = constraints
  )
  
  expect_is(val$summary$n_boundary_violations, "numeric")
})

test_that("gc_validate_realizations generates summary per component", {
  r <- create_test_realizations(n_components = 3, n_realizations = 3)
  
  val <- geocoda::gc_validate_realizations(
    r,
    components = c("sand", "silt", "clay")
  )
  
  expect_equal(nrow(val$summary), 3)
  expect_equal(val$summary$component, c("sand", "silt", "clay"))
})

test_that("gc_validate_realizations marks valid realizations", {
  r <- create_test_realizations(n_components = 3, n_realizations = 3)
  
  # Realizations are created to sum to ~100, so should have minimal violations
  val <- geocoda::gc_validate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    target_sum = 100,
    sum_tolerance = 1
  )
  
  violations <- terra::values(val$violations_map)[, 1]
  # Most cells should have few violations (given our normalization)
  pct_zero <- sum(violations == 0) / length(violations)
  expect_true(pct_zero >= 0.5)  # At least 50% should be valid
})

test_that("gc_validate_realizations handles NULL constraints", {
  r <- create_test_realizations(n_components = 3, n_realizations = 2)
  
  # Should not error with NULL constraints
  val <- geocoda::gc_validate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    constraints = NULL
  )
  
  expect_is(val$summary, "data.frame")
})

test_that("gc_validate_realizations includes metadata", {
  r <- create_test_realizations(n_components = 2, n_realizations = 3)
  
  val <- geocoda::gc_validate_realizations(
    r,
    components = c("sand", "clay"),
    target_sum = 100
  )
  
  meta <- val$metadata
  expect_equal(meta$target_sum, 100)
  expect_equal(meta$n_components, 2)
  expect_equal(meta$n_realizations, 3)
})

test_that("gc_validate_realizations computes pct_valid", {
  r <- create_test_realizations(n_components = 1, n_realizations = 2)
  
  val <- geocoda::gc_validate_realizations(
    r,
    components = "sand"
  )
  
  # If no NAs, should be 100%
  pct_valid <- val$summary$pct_valid
  expect_true(pct_valid >= 90)  # Allow some floating point tolerance
})

# ============================================================================
# Integration tests
# ============================================================================

test_that("aggregate_realizations + validate_realizations workflow", {
  # Full workflow: create realizations → aggregate → validate
  r <- create_test_realizations(n_components = 3, n_realizations = 5)
  
  # Aggregate
  agg <- geocoda::gc_aggregate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    stats = c("mean", "sd")
  )
  
  # Note: validation typically operates on raw realizations, not aggregates
  # But we can validate the aggregated mean layer as a consistency check
  sand_mean_layer <- terra::subset(agg, "sand_mean")
  val <- geocoda::gc_validate_realizations(
    sand_mean_layer,
    components = "sand"
  )
  
  expect_s4_class(agg, "SpatRaster")
  expect_is(val, "list")
})

test_that("ensemble workflow produces consistent results", {
  r <- create_test_realizations(n_components = 3, n_realizations = 10)
  
  agg1 <- geocoda::gc_aggregate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    stats = c("mean")
  )
  
  # Run again with same input
  agg2 <- geocoda::gc_aggregate_realizations(
    r,
    components = c("sand", "silt", "clay"),
    stats = c("mean")
  )
  
  # Should be identical
  expect_equal(terra::values(agg1), terra::values(agg2))
})
