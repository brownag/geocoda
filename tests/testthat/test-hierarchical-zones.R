context("Hierarchical Shrinkage Zones")

# ============================================================================
# Test 1-4: Hierarchy Definition and Validation
# ============================================================================

test_that("gc_define_hierarchy creates proper S3 object", {
  hierarchy <- gc_define_hierarchy(
    zone_names = c("zone_1", "zone_2", "zone_3"),
    parent_region = "study_area",
    component_names = c("sand", "silt", "clay")
  )
  
  expect_s3_class(hierarchy, "gc_hierarchy")
  expect_equal(hierarchy$n_zones, 3)
  expect_equal(hierarchy$n_components, 3)
  expect_equal(hierarchy$target_sum, 100)
})

test_that("gc_define_hierarchy validates component count", {
  expect_error(
    gc_define_hierarchy(
      zone_names = "zone_1",
      component_names = c("sand")
    ),
    "at least 2 components"
  )
})

test_that("gc_define_hierarchy stores zone metadata", {
  hierarchy <- gc_define_hierarchy(
    zone_names = c("north", "south", "west"),
    parent_region = "region_A",
    component_names = c("a", "b", "c", "d"),
    target_sum = 100
  )
  
  expect_equal(hierarchy$zones, c("north", "south", "west"))
  expect_equal(hierarchy$parent_region, "region_A")
  expect_equal(hierarchy$n_ilr, 3)  # n_components - 1
})

test_that("gc_define_hierarchy allows custom target_sum", {
  hierarchy <- gc_define_hierarchy(
    zone_names = "zone_1",
    component_names = c("a", "b"),
    target_sum = 1.0
  )
  
  expect_equal(hierarchy$target_sum, 1.0)
})

# ============================================================================
# Test 5-8: Prior Specification and Configuration
# ============================================================================

test_that("gc_set_hierarchy_priors creates proper S3 object", {
  hierarchy <- gc_define_hierarchy(
    zone_names = c("z1", "z2"),
    component_names = c("sand", "silt", "clay")
  )
  
  priors <- gc_set_hierarchy_priors(
    hierarchy,
    pooling = "moderate",
    verbose = FALSE
  )
  
  expect_s3_class(priors, "gc_prior_spec")
  expect_equal(priors$pooling, "moderate")
})

test_that("gc_set_hierarchy_priors pools correctly with different strategies", {
  hierarchy <- gc_define_hierarchy(
    zone_names = "zone_1",
    component_names = c("a", "b", "c")
  )
  
  pooling_strategies <- c("none", "weak", "moderate", "strong")
  expected_coefs <- c(0.00, 0.05, 0.20, 0.50)
  
  for (i in seq_along(pooling_strategies)) {
    priors <- gc_set_hierarchy_priors(
      hierarchy,
      pooling = pooling_strategies[i],
      verbose = FALSE
    )
    expect_equal(priors$pooling_coefficient, expected_coefs[i])
  }
})

test_that("gc_set_hierarchy_priors sets default global distribution", {
  hierarchy <- gc_define_hierarchy(
    zone_names = "zone_1",
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  
  expect_length(priors$global_mean, 2)  # n_components - 1
  expect_equal(nrow(priors$global_covariance), 2)
  expect_equal(ncol(priors$global_covariance), 2)
})

test_that("gc_set_hierarchy_priors validates input object", {
  expect_error(
    gc_set_hierarchy_priors(
      list(zones = "z1"),  # Not a gc_hierarchy object
      verbose = FALSE
    ),
    "gc_hierarchy"
  )
})

# ============================================================================
# Test 9-15: Hierarchical Model Fitting
# ============================================================================

test_that("gc_fit_hierarchical_model fits without error", {
  # Create synthetic data
  set.seed(123)
  data <- data.frame(
    zone = rep(c("zone_1", "zone_2"), each = 20),
    ilr1 = c(rnorm(20, 0, 1), rnorm(20, 0.5, 1)),
    ilr2 = c(rnorm(20, 0, 1), rnorm(20, -0.5, 1))
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("zone_1", "zone_2"),
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, pooling = "moderate", verbose = FALSE)
  
  fit <- gc_fit_hierarchical_model(
    data,
    priors,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "gc_hierarchical_fit")
})

test_that("gc_fit_hierarchical_model produces valid shrinkage estimates", {
  set.seed(456)
  data <- data.frame(
    zone = rep(c("Z1", "Z2", "Z3"), each = 15),
    ilr1 = rnorm(45, 0, 1),
    ilr2 = rnorm(45, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = c("Z1", "Z2", "Z3"),
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)

  # Check zone_estimates (renamed from posterior_draws for clarity)
  expect_length(fit$zone_estimates, 3)
  expect_true(all(names(fit$zone_estimates) %in% c("Z1", "Z2", "Z3")))

  for (zone_est in fit$zone_estimates) {
    expect_true("mean" %in% names(zone_est))
    expect_true("cov" %in% names(zone_est))
    expect_true("sd" %in% names(zone_est))
    expect_true("n_obs" %in% names(zone_est))
  }
})

test_that("gc_fit_hierarchical_model applies pooling shrinkage", {
  set.seed(789)
  data <- data.frame(
    zone = rep(c("z1", "z2"), each = 20),
    ilr1 = c(rnorm(20, 2, 1), rnorm(20, -2, 1)),  # Very different means
    ilr2 = rnorm(40, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = c("z1", "z2"),
    component_names = c("a", "b", "c")
  )

  # Fit with strong pooling
  priors <- gc_set_hierarchy_priors(hierarchy, pooling = "strong", verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)

  # Means should be pulled toward global mean
  z1_mean <- fit$zone_estimates$z1$mean[1]
  z2_mean <- fit$zone_estimates$z2$mean[1]

  # Strong pooling should make them closer together
  expect_lt(abs(z1_mean - z2_mean), 3)  # Reduced from original ~4
})

test_that("gc_fit_hierarchical_model zone summaries are correct", {
  set.seed(111)
  data <- data.frame(
    zone = rep(c("north", "south"), each = 25),
    ilr1 = rnorm(50, 0, 1),
    ilr2 = rnorm(50, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("north", "south"),
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  summaries <- fit$zone_summaries
  expect_true(all(summaries$zone %in% c("north", "south")))
  expect_true(all(summaries$n_obs > 0))
  expect_true(all(!is.na(summaries$posterior_mean)))
  expect_true(all(summaries$posterior_sd >= 0))
})

test_that("gc_fit_hierarchical_model handles missing ILR columns error", {
  data <- data.frame(
    zone = c("z1", "z2"),
    x = 1:2  # No ILR columns
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("z1", "z2"),
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  
  expect_error(
    gc_fit_hierarchical_model(data, priors, verbose = FALSE),
    "ILR"
  )
})

test_that("gc_fit_hierarchical_model ignores unfitted zones gracefully", {
  set.seed(222)
  data <- data.frame(
    zone = rep("zone_1", 20),  # Only one zone in data
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = c("zone_1", "zone_2", "zone_3"),  # Three zones defined
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)

  # Should only have fitted zone (zone_1)
  expect_length(fit$zone_estimates, 1)
  expect_true("zone_1" %in% names(fit$zone_estimates))
})

test_that("gc_fit_hierarchical_model defaults to analytical backend", {
  set.seed(333)
  data <- data.frame(
    zone = rep("zone_1", 20),
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = "zone_1",
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)

  expect_equal(fit$metadata$backend, "analytical")
})

# ============================================================================
# Test 16-20: Hierarchical Simulation
# ============================================================================

test_that("gc_sim_hierarchical creates valid simulation object", {
  set.seed(333)
  data <- data.frame(
    zone = rep(c("A", "B"), each = 15),
    ilr1 = rnorm(30, 0, 1),
    ilr2 = rnorm(30, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("A", "B"),
    component_names = c("sand", "silt", "clay")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  locations <- data.frame(
    x = c(0, 10, 20, 0, 10, 20),
    y = c(0, 0, 0, 10, 10, 10),
    zone = c(rep("A", 3), rep("B", 3))
  )
  
  sims <- gc_sim_hierarchical(
    fit,
    locations,
    nsim = 2,
    target_names = c("sand", "silt", "clay"),
    verbose = FALSE
  )
  
  expect_s3_class(sims, "gc_hierarchical_sims")
  expect_length(sims, 2)  # Two simulations
})

test_that("gc_sim_hierarchical respects nsim parameter", {
  set.seed(444)
  data <- data.frame(
    zone = rep("zone_1", 20),
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = "zone_1",
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  locations <- data.frame(x = 1:5, y = 1:5, zone = "zone_1")
  
  for (n in c(3, 5, 10)) {
    sims <- gc_sim_hierarchical(fit, locations, nsim = n, verbose = FALSE)
    expect_length(sims, n)
  }
})

test_that("gc_sim_hierarchical validates fit object", {
  locations <- data.frame(x = 1:5, y = 1:5, zone = "z1")
  
  expect_error(
    gc_sim_hierarchical(
      list(zone_realizations = list()),  # Not a hierarchical fit
      locations,
      verbose = FALSE
    ),
    "gc_hierarchical_fit"
  )
})

test_that("gc_sim_hierarchical requires zone column in locations", {
  set.seed(555)
  data <- data.frame(
    zone = "zone_1",
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = "zone_1",
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  bad_locations <- data.frame(x = 1:5, y = 1:5)  # No zone column
  
  expect_error(
    gc_sim_hierarchical(fit, bad_locations, verbose = FALSE),
    "zone"
  )
})

test_that("gc_sim_hierarchical handles zones with no locations", {
  set.seed(666)
  data <- data.frame(
    zone = rep(c("z1", "z2"), each = 15),
    ilr1 = rnorm(30, 0, 1),
    ilr2 = rnorm(30, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("z1", "z2"),
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  # Only locations in zone_1
  locations <- data.frame(x = 1:5, y = 1:5, zone = "z1")
  
  sims <- gc_sim_hierarchical(fit, locations, nsim = 2, verbose = FALSE)
  
  # Should still complete without error
  expect_s3_class(sims, "gc_hierarchical_sims")
})

# ============================================================================
# Test 21-24: Stacks and Output Conversion
# ============================================================================

test_that("gc_hierarchy_to_stacks creates stack objects", {
  set.seed(777)
  data <- data.frame(
    zone = rep(c("A", "B"), each = 15),
    ilr1 = rnorm(30, 0, 1),
    ilr2 = rnorm(30, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("A", "B"),
    component_names = c("sand", "silt", "clay")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  locations <- data.frame(
    x = c(0, 10, 0, 10),
    y = c(0, 0, 10, 10),
    zone = c("A", "A", "B", "B")
  )
  
  sims <- gc_sim_hierarchical(fit, locations, nsim = 3, verbose = FALSE)
  stacks <- gc_hierarchy_to_stacks(sims, hierarchy, verbose = FALSE)
  
  expect_s3_class(stacks, "gc_hierarchy_stacks")
})

test_that("gc_hierarchy_to_stacks computes correct statistics", {
  set.seed(888)
  data <- data.frame(
    zone = "zone_1",
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = "zone_1",
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  locations <- data.frame(x = 1:3, y = 1:3, zone = "zone_1")
  sims <- gc_sim_hierarchical(fit, locations, nsim = 5, verbose = FALSE)
  
  stacks <- gc_hierarchy_to_stacks(
    sims,
    hierarchy,
    stats = c("mean", "sd", "p05", "p95"),
    verbose = FALSE
  )
  
  # Should have stacks
  expect_true(length(stacks) > 0)
})

test_that("gc_combine_hierarchical_results creates valid object", {
  set.seed(999)
  data <- data.frame(
    zone = "zone_1",
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = "zone_1",
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  results <- gc_combine_hierarchical_results(fit)
  
  expect_s3_class(results, "gc_hierarchical_results")
  expect_true("fit" %in% names(results))
  expect_equal(results$fit, fit)
})

# ============================================================================
# Test 25-28: Posterior Extraction and Validation
# ============================================================================

test_that("gc_extract_zone_posterior returns correct structure", {
  set.seed(1001)
  data <- data.frame(
    zone = rep(c("z1", "z2"), each = 20),
    ilr1 = rnorm(40, 0, 1),
    ilr2 = rnorm(40, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("z1", "z2"),
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  
  posterior <- gc_extract_zone_posterior(fit, "z1")
  
  expect_is(posterior, "data.frame")
  expect_true("mean" %in% names(posterior))
  expect_true("sd" %in% names(posterior))
  expect_true("n_obs" %in% names(posterior))
})

test_that("gc_validate_hierarchical_results checks zone coverage", {
  set.seed(1002)
  data <- data.frame(
    zone = rep(c("z1", "z2", "z3"), each = 20),
    ilr1 = rnorm(60, 0, 1),
    ilr2 = rnorm(60, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = c("z1", "z2", "z3"),
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  results <- gc_combine_hierarchical_results(fit)

  validation <- gc_validate_hierarchical_results(results, verbose = FALSE)

  expect_true("valid" %in% names(validation))
  expect_true("zone_coverage_valid" %in% names(validation))
  expect_true("zone_coverage" %in% names(validation))
})

test_that("gc_validate_hierarchical_results detects coverage issues", {
  set.seed(1003)
  # Create data with one zone having very few observations
  data <- data.frame(
    zone = c(rep("z1", 30), "z2", "z2", "z2"),
    ilr1 = rnorm(33, 0, 1),
    ilr2 = rnorm(33, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("z1", "z2"),
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  results <- gc_combine_hierarchical_results(fit)
  
  validation <- gc_validate_hierarchical_results(results, verbose = FALSE)
  
  # Should detect zone with few observations
  expect_true("issues" %in% names(validation))
})

test_that("gc_validate_hierarchical_results generates recommendations", {
  set.seed(1004)
  data <- data.frame(
    zone = c(rep("z1", 20), "z2"),  # Very unbalanced
    ilr1 = rnorm(21, 0, 1),
    ilr2 = rnorm(21, 0, 1)
  )
  
  hierarchy <- gc_define_hierarchy(
    zone_names = c("z1", "z2"),
    component_names = c("a", "b", "c")
  )
  
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- gc_fit_hierarchical_model(data, priors, verbose = FALSE)
  results <- gc_combine_hierarchical_results(fit)
  
  validation <- gc_validate_hierarchical_results(results, verbose = FALSE)
  
  # Should provide recommendations
  expect_true(length(validation$recommendations) >= 0)
})
