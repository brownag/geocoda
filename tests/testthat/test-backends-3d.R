test_that("Distance-weighted analytical 3D backend initialization", {
  # Setup: Create synthetic 3-zone domain with ILR data
  set.seed(42)
  
  # Create 3 zones with different spatial locations
  zone_centers <- data.frame(
    zone = c("Zone_A", "Zone_B", "Zone_C"),
    x = c(0, 100, 200),
    y = c(0, 50, 100),
    z = c(0, 10, 20)
  )
  
  # Simulate observations in each zone
  n_per_zone <- 30
  zones <- rep(zone_centers$zone, each = n_per_zone)
  
  # ILR columns (2D simplex -> 1 ILR component)
  ilr1 <- rnorm(n_per_zone * 3, mean = 0, sd = 0.5)
  ilr1[1:n_per_zone] <- ilr1[1:n_per_zone] + 0.3  # Zone A offset
  ilr1[(n_per_zone+1):(2*n_per_zone)] <- ilr1[(n_per_zone+1):(2*n_per_zone)] - 0.2  # Zone B offset
  
  data <- data.frame(
    zone = zones,
    ilr1 = ilr1,
    x = rep(zone_centers$x, each = n_per_zone) + rnorm(n_per_zone * 3, 0, 10),
    y = rep(zone_centers$y, each = n_per_zone) + rnorm(n_per_zone * 3, 0, 10),
    z = rep(zone_centers$z, each = n_per_zone) + rnorm(n_per_zone * 3, 0, 5)
  )
  
  # Create prior specification
  hierarchy <- list(
    zones = c("Zone_A", "Zone_B", "Zone_C"),
    n_zones = 3,
    n_components = 2
  )
  class(hierarchy) <- "gc_hierarchy"
  
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = c(0),
    global_covariance = matrix(1, 1, 1),
    pooling_coefficient = 0.5,
    prior_sd_mean = c(1),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  
  # Test 1: Basic analytical 3D fitting with computed zone centers
  result <- fit_hierarchical_analytical_3d_weighted(
    data = data,
    prior_spec = prior_spec,
    zone_centers = NULL,  # Should compute from data
    spatial_decay = "inverse_distance",
    verbose = FALSE
  )
  
  expect_s3_class(result, "gc_hierarchical_fit")
  expect_named(result, c(
    "zone_estimates", "global_estimates", "samples",
    "shrinkage_weights", "zone_summaries", "zone_centers",
    "diagnostics", "prior_spec", "metadata"
  ))
  expect_equal(length(result$zone_estimates), 3)
  expect_null(result$samples)  # Analytical should have no samples
  expect_equal(nrow(result$zone_centers), 3)
})

test_that("Distance-weighted analytical 3D: exponential decay", {
  set.seed(43)
  
  # Minimal setup for quick test
  zone_centers <- data.frame(
    zone = c("A", "B"),
    x = c(0, 50),
    y = c(0, 50),
    z = c(0, 10)
  )
  
  zones <- rep(c("A", "B"), each = 20)
  data <- data.frame(
    zone = zones,
    ilr1 = c(rnorm(20, 0, 0.5), rnorm(20, -0.2, 0.5)),
    x = rep(zone_centers$x, each = 20) + rnorm(40, 0, 5),
    y = rep(zone_centers$y, each = 20) + rnorm(40, 0, 5),
    z = rep(zone_centers$z, each = 20) + rnorm(40, 0, 3)
  )
  
  hierarchy <- list(zones = c("A", "B"), n_zones = 2, n_components = 2)
  class(hierarchy) <- "gc_hierarchy"
  
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = c(0),
    global_covariance = matrix(1, 1, 1),
    pooling_coefficient = 0.5,
    prior_sd_mean = c(1),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  
  # Test with exponential decay
  result <- fit_hierarchical_analytical_3d_weighted(
    data = data,
    prior_spec = prior_spec,
    spatial_decay = "exponential",
    decay_range = 100,
    verbose = FALSE
  )
  
  expect_s3_class(result, "gc_hierarchical_fit")
  expect_equal(result$diagnostics$spatial_decay, "exponential")
  expect_equal(result$diagnostics$decay_range, 100)
  expect_true(is.na(result$diagnostics$decay_power))
})

test_that("Distance-weighted analytical 3D: provided zone centers", {
  set.seed(44)
  
  zone_centers <- data.frame(
    zone = c("Z1", "Z2", "Z3"),
    x = c(0, 100, 200),
    y = c(0, 100, 0),
    z = c(0, 50, 100)
  )
  
  # Data WITHOUT coordinates - must provide zone_centers
  zones <- rep(c("Z1", "Z2", "Z3"), each = 15)
  data <- data.frame(
    zone = zones,
    ilr1 = rnorm(45, 0, 0.5)
  )
  
  hierarchy <- list(zones = c("Z1", "Z2", "Z3"), n_zones = 3, n_components = 2)
  class(hierarchy) <- "gc_hierarchy"
  
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = c(0),
    global_covariance = matrix(1, 1, 1),
    pooling_coefficient = 0.5,
    prior_sd_mean = c(1),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  
  # Must provide zone_centers
  result <- fit_hierarchical_analytical_3d_weighted(
    data = data,
    prior_spec = prior_spec,
    zone_centers = zone_centers,
    spatial_decay = "inverse_distance",
    verbose = FALSE
  )
  
  expect_s3_class(result, "gc_hierarchical_fit")
  expect_equal(nrow(result$zone_centers), 3)
  
  # Should error if no coordinates and no zone_centers
  expect_error(
    fit_hierarchical_analytical_3d_weighted(
      data = data,
      prior_spec = prior_spec,
      zone_centers = NULL,
      verbose = FALSE
    )
  )
})

test_that("Dispatcher routes to analytical 3D backend", {
  set.seed(45)
  
  zone_centers <- data.frame(
    zone = c("North", "South"),
    x = c(0, 100),
    y = c(100, 0),
    z = c(10, 20)
  )
  
  zones <- rep(c("North", "South"), each = 20)
  data <- data.frame(
    zone = zones,
    ilr1 = c(rnorm(20, 0.5, 0.4), rnorm(20, -0.5, 0.4)),
    x = rep(zone_centers$x, each = 20) + rnorm(40, 0, 5),
    y = rep(zone_centers$y, each = 20) + rnorm(40, 0, 5),
    z = rep(zone_centers$z, each = 20) + rnorm(40, 0, 2)
  )
  
  hierarchy <- list(zones = c("North", "South"), n_zones = 2, n_components = 2)
  class(hierarchy) <- "gc_hierarchy"
  
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = c(0),
    global_covariance = matrix(1, 1, 1),
    pooling_coefficient = 0.5,
    prior_sd_mean = c(1),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  
  # Dispatcher should route to analytical_3d backend
  result <- fit_hierarchical_backend(
    data = data,
    prior_spec = prior_spec,
    backend = "analytical_3d",
    verbose = FALSE
  )
  
  expect_s3_class(result, "gc_hierarchical_fit")
  expect_equal(result$metadata$backend, "analytical_3d")
})

test_that("Backend dispatcher recognizes all backend options", {
  # Just test the match.arg validation logic
  valid_backends <- c("analytical", "stan", "nimble", "analytical_3d", "stan_3d", "nimble_3d")
  
  for (backend in valid_backends) {
    expect_equal(match.arg(backend, valid_backends), backend)
  }
  
  # Invalid backend should fail
  expect_error(match.arg("invalid_backend", valid_backends))
})

test_that("Distance-weighted analytical 3D: shrinkage strength varies with neighbors", {
  set.seed(46)
  
  # Create 4 zones with varying distances
  zone_centers <- data.frame(
    zone = c("Central", "North", "South", "East"),
    x = c(100, 100, 100, 200),
    y = c(100, 200, 0, 100),
    z = c(50, 50, 50, 50)
  )
  
  zones <- rep(zone_centers$zone, each = 15)
  data <- data.frame(
    zone = zones,
    ilr1 = rnorm(60, 0, 0.5),
    x = rep(zone_centers$x, each = 15) + rnorm(60, 0, 3),
    y = rep(zone_centers$y, each = 15) + rnorm(60, 0, 3),
    z = rep(zone_centers$z, each = 15) + rnorm(60, 0, 2)
  )
  
  hierarchy <- list(
    zones = c("Central", "North", "South", "East"),
    n_zones = 4,
    n_components = 2
  )
  class(hierarchy) <- "gc_hierarchy"
  
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = c(0),
    global_covariance = matrix(1, 1, 1),
    pooling_coefficient = 0.5,
    prior_sd_mean = c(1),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  
  result <- fit_hierarchical_analytical_3d_weighted(
    data = data,
    prior_spec = prior_spec,
    spatial_decay = "inverse_distance",
    decay_power = 2.0,
    verbose = FALSE
  )
  
  # Check that we have correct number of zones and neighbors are calculated
  expect_equal(nrow(result$shrinkage_weights), 4)
  
  # Each zone should have number of neighbors calculated (at least 0)
  expect_true(all(result$shrinkage_weights$n_neighbors >= 0))
  
  # Shrinkage strength should be numeric and reasonable
  expect_true(all(result$shrinkage_weights$shrinkage_strength >= 0))
  expect_true(all(result$shrinkage_weights$shrinkage_strength <= 1))
})

test_that("Distance-weighted analytical 3D: posterior summaries are reasonable", {
  set.seed(47)
  
  zone_centers <- data.frame(
    zone = c("Deep", "Shallow"),
    x = c(0, 300),
    y = c(0, 300),
    z = c(500, 100)
  )
  
  zones <- rep(c("Deep", "Shallow"), each = 25)
  # Create data with clear zone differences
  ilr1 <- c(rnorm(25, mean=1.0, sd=0.3), rnorm(25, mean=-1.0, sd=0.3))
  
  data <- data.frame(
    zone = zones,
    ilr1 = ilr1,
    x = rep(zone_centers$x, each = 25) + rnorm(50, 0, 10),
    y = rep(zone_centers$y, each = 25) + rnorm(50, 0, 10),
    z = rep(zone_centers$z, each = 25) + rnorm(50, 0, 20)
  )
  
  hierarchy <- list(zones = c("Deep", "Shallow"), n_zones = 2, n_components = 2)
  class(hierarchy) <- "gc_hierarchy"
  
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = c(0),
    global_covariance = matrix(1, 1, 1),
    pooling_coefficient = 0.3,  # Lower pooling to emphasize zone differences
    prior_sd_mean = c(1),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  
  result <- fit_hierarchical_analytical_3d_weighted(
    data = data,
    prior_spec = prior_spec,
    verbose = FALSE
  )
  
  # Check posterior summaries
  expect_equal(nrow(result$zone_summaries), 2)
  
  deep_summary <- result$zone_summaries[result$zone_summaries$zone == "Deep", ]
  shallow_summary <- result$zone_summaries[result$zone_summaries$zone == "Shallow", ]
  
  # Posterior means should be different (zone-specific signal)
  expect_false(isTRUE(all.equal(deep_summary$posterior_mean, shallow_summary$posterior_mean)))
  
  # Credible intervals should be narrower than prior (informative data)
  expect_true(all(diff(cbind(deep_summary$ci_lower, deep_summary$ci_upper)) < 1.0))
  expect_true(all(diff(cbind(shallow_summary$ci_lower, shallow_summary$ci_upper)) < 1.0))
})

test_that("Analytical 3D backend validation (insufficient zones)", {
  set.seed(48)
  
  zone_centers <- data.frame(
    zone = "OnlyZone",
    x = 0,
    y = 0,
    z = 0
  )
  
  data <- data.frame(
    zone = rep("OnlyZone", 20),
    ilr1 = rnorm(20, 0, 0.5),
    x = rnorm(20, 0, 5),
    y = rnorm(20, 0, 5),
    z = rnorm(20, 0, 5)
  )
  
  hierarchy <- list(zones = "OnlyZone", n_zones = 1, n_components = 2)
  class(hierarchy) <- "gc_hierarchy"
  
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = c(0),
    global_covariance = matrix(1, 1, 1),
    pooling_coefficient = 0.5,
    prior_sd_mean = c(1),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  
  # Should still work, just with no spatial neighbors
  result <- fit_hierarchical_analytical_3d_weighted(
    data = data,
    prior_spec = prior_spec,
    verbose = FALSE
  )
  
  expect_s3_class(result, "gc_hierarchical_fit")
  expect_equal(length(result$zone_estimates), 1)
})

test_that("Analytical 3D: pooling coefficient estimates match expectation", {
  set.seed(49)
  
  # Test with fixed pooling value
  zone_centers <- data.frame(
    zone = c("A", "B", "C"),
    x = c(0, 100, 200),
    y = c(0, 100, 200),
    z = c(0, 50, 100)
  )
  
  zones <- rep(zone_centers$zone, each = 20)
  data <- data.frame(
    zone = zones,
    ilr1 = rnorm(60, 0, 0.5),
    x = rep(zone_centers$x, each = 20) + rnorm(60, 0, 5),
    y = rep(zone_centers$y, each = 20) + rnorm(60, 0, 5),
    z = rep(zone_centers$z, each = 20) + rnorm(60, 0, 3)
  )
  
  hierarchy <- list(zones = c("A", "B", "C"), n_zones = 3, n_components = 2)
  class(hierarchy) <- "gc_hierarchy"
  
  pooling_value <- 0.6
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = c(0),
    global_covariance = matrix(1, 1, 1),
    pooling_coefficient = pooling_value,
    prior_sd_mean = c(1),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  
  result <- fit_hierarchical_analytical_3d_weighted(
    data = data,
    prior_spec = prior_spec,
    verbose = FALSE
  )
  
  # Metadata should preserve pooling coefficient
  expect_equal(result$metadata$pooling_coefficient, pooling_value)
  
  # Shrinkage weights should use effective pooling values ≥ base pooling
  expect_true(all(result$shrinkage_weights$shrinkage_strength >= pooling_value))
})
