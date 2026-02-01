context("Hierarchical Model Backend System")

# ============================================================================
# Test 1-4: Backend Availability Checks
# ============================================================================

test_that("check_backend_available recognizes analytical backend", {
  available <- check_backend_available("analytical")
  expect_true(available)
})

test_that("check_backend_available handles unavailable Stan backend", {
  # Stan may or may not be available, but function should return logical
  available <- check_backend_available("stan")
  expect_is(available, "logical")
})

test_that("check_backend_available handles unavailable Nimble backend", {
  # Nimble may or may not be available, but function should return logical
  available <- check_backend_available("nimble")
  expect_is(available, "logical")
})

test_that("check_backend_available returns FALSE for unknown backend", {
  available <- check_backend_available("unknown_backend")
  expect_false(available)
})

# ============================================================================
# Test 5-7: Backend Installation Messages
# ============================================================================

test_that("get_backend_install_msg provides Stan installation instructions", {
  msg <- get_backend_install_msg("stan")
  expect_match(msg, "rstan")
  expect_match(msg, "install.packages")
})

test_that("get_backend_install_msg provides Nimble installation instructions", {
  msg <- get_backend_install_msg("nimble")
  expect_match(msg, "nimble")
  expect_match(msg, "install.packages")
})

test_that("get_backend_install_msg handles unknown backend", {
  msg <- get_backend_install_msg("unknown")
  expect_match(msg, "Unknown")
})

# ============================================================================
# Test 8-12: Backend Dispatcher Function
# ============================================================================

test_that("fit_hierarchical_backend with analytical backend works", {
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

  fit <- fit_hierarchical_backend(
    data,
    priors,
    backend = "analytical",
    verbose = FALSE
  )

  expect_s3_class(fit, "gc_hierarchical_fit")
})

test_that("fit_hierarchical_backend returns standardized structure", {
  set.seed(456)
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
  fit <- fit_hierarchical_backend(data, priors, backend = "analytical", verbose = FALSE)

  # Check standardized structure
  expect_true("zone_estimates" %in% names(fit))
  expect_true("global_estimates" %in% names(fit))
  expect_true("samples" %in% names(fit))
  expect_true("diagnostics" %in% names(fit))
  expect_true("prior_spec" %in% names(fit))
  expect_true("metadata" %in% names(fit))

  # Analytical backend should have NULL samples
  expect_null(fit$samples)

  # Metadata should include backend info
  expect_equal(fit$metadata$backend, "analytical")
  expect_equal(fit$metadata$method, "analytical_shrinkage")
})

test_that("fit_hierarchical_backend with invalid backend raises error", {
  data <- data.frame(
    zone = "z1",
    ilr1 = rnorm(10),
    ilr2 = rnorm(10)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = "z1",
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  expect_error(
    fit_hierarchical_backend(data, priors, backend = "invalid", verbose = FALSE),
    "should be one of"
  )
})

test_that("fit_hierarchical_backend with unavailable Stan backend raises error", {
  skip_if(requireNamespace("rstan", quietly = TRUE), "Test only runs if rstan not installed")

  data <- data.frame(
    zone = "z1",
    ilr1 = rnorm(10),
    ilr2 = rnorm(10)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = "z1",
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  expect_error(
    fit_hierarchical_backend(data, priors, backend = "stan", verbose = FALSE),
    "not available"
  )
})

# ============================================================================
# Test 13-16: Analytical Backend Direct Function
# ============================================================================

test_that("fit_hierarchical_analytical produces valid results", {
  set.seed(789)
  data <- data.frame(
    zone = rep(c("Z1", "Z2"), each = 20),
    ilr1 = rnorm(40, 0, 1),
    ilr2 = rnorm(40, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = c("Z1", "Z2"),
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- fit_hierarchical_analytical(data, priors, verbose = FALSE)

  expect_s3_class(fit, "gc_hierarchical_fit")
  expect_equal(fit$metadata$backend, "analytical")
})

test_that("fit_hierarchical_analytical matches expected zone count", {
  set.seed(234)
  data <- data.frame(
    zone = rep(c("a", "b", "c"), each = 10),
    ilr1 = rnorm(30),
    ilr2 = rnorm(30)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = c("a", "b", "c"),
    component_names = c("x", "y", "z")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- fit_hierarchical_analytical(data, priors, verbose = FALSE)

  expect_length(fit$zone_estimates, 3)
  expect_equal(fit$metadata$n_zones_fitted, 3)
})

test_that("fit_hierarchical_analytical diagnostics marked as analytical", {
  set.seed(567)
  data <- data.frame(
    zone = "z1",
    ilr1 = rnorm(20),
    ilr2 = rnorm(20)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = "z1",
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- fit_hierarchical_analytical(data, priors, verbose = FALSE)

  expect_equal(fit$diagnostics$backend, "analytical")
  expect_match(fit$diagnostics$convergence, "analytical")
})

test_that("fit_hierarchical_analytical handles single zone", {
  set.seed(345)
  data <- data.frame(
    zone = rep("single_zone", 25),
    ilr1 = rnorm(25),
    ilr2 = rnorm(25)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = "single_zone",
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)
  fit <- fit_hierarchical_analytical(data, priors, verbose = FALSE)

  expect_length(fit$zone_estimates, 1)
  expect_true("single_zone" %in% names(fit$zone_estimates))
})

# ============================================================================
# Test 17-19: Stan Skeleton (should raise informative error)
# ============================================================================

# Note: Removed obsolete "not yet implemented" tests since backends are now fully implemented

# ============================================================================
# Test 20: Default backend selection
# ============================================================================

test_that("fit_hierarchical_backend defaults to analytical backend", {
  set.seed(890)
  data <- data.frame(
    zone = "z1",
    ilr1 = rnorm(20),
    ilr2 = rnorm(20)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = "z1",
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  # Call without specifying backend
  fit <- fit_hierarchical_backend(data, priors, verbose = FALSE)

  # Should default to analytical
  expect_equal(fit$metadata$backend, "analytical")
})

# ============================================================================
# Test 21-25: Stan Backend (full implementation)
# ============================================================================

test_that("Stan backend fits successfully with rstan available", {
  skip_if_not(requireNamespace("rstan", quietly = TRUE), "rstan not available")

  set.seed(456)
  data <- data.frame(
    zone = rep(c("Z1", "Z2"), each = 15),
    ilr1 = rnorm(30, 0, 1),
    ilr2 = rnorm(30, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = c("Z1", "Z2"),
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, pooling = "moderate", verbose = FALSE)

  fit_stan <- fit_hierarchical_stan(
    data,
    priors,
    n_iter = 500,
    n_warmup = 250,
    n_chains = 1,
    verbose = FALSE
  )

  expect_s3_class(fit_stan, "gc_hierarchical_fit")
  expect_equal(fit_stan$metadata$backend, "stan")
})

test_that("Stan backend returns standardized structure", {
  skip_if_not(requireNamespace("rstan", quietly = TRUE), "rstan not available")

  set.seed(789)
  data <- data.frame(
    zone = rep("Z1", 20),
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = "Z1",
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  fit_stan <- fit_hierarchical_stan(data, priors, n_iter = 400, n_warmup = 200, n_chains = 1, verbose = FALSE)

  # Check standardized structure
  expect_true("zone_estimates" %in% names(fit_stan))
  expect_true("global_estimates" %in% names(fit_stan))
  expect_true("samples" %in% names(fit_stan))
  expect_true("diagnostics" %in% names(fit_stan))
  expect_true("prior_spec" %in% names(fit_stan))
  expect_true("metadata" %in% names(fit_stan))

  # Stan should have samples (not NULL like analytical)
  expect_false(is.null(fit_stan$samples))
  expect_true(is.matrix(fit_stan$samples) || is.data.frame(fit_stan$samples))
})

test_that("Stan backend computes convergence diagnostics", {
  skip_if_not(requireNamespace("rstan", quietly = TRUE), "rstan not available")

  set.seed(321)
  data <- data.frame(
    zone = "Z1",
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(zone_names = "Z1", component_names = c("a", "b", "c"))
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  fit_stan <- fit_hierarchical_stan(data, priors, n_iter = 400, n_warmup = 200, n_chains = 1, verbose = FALSE)

  # Check diagnostics
  expect_true("rhat" %in% names(fit_stan$diagnostics))
  expect_true("ess_bulk" %in% names(fit_stan$diagnostics))
  expect_true("n_divergent" %in% names(fit_stan$diagnostics))
  expect_true("adapt_delta" %in% names(fit_stan$diagnostics))
})

test_that("Stan backend handles estimate_pooling option", {
  skip_if_not(requireNamespace("rstan", quietly = TRUE), "rstan not available")

  set.seed(654)
  data <- data.frame(
    zone = rep("Z1", 20),
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(zone_names = "Z1", component_names = c("a", "b", "c"))
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  # Test with estimate_pooling = TRUE
  fit_stan_est <- fit_hierarchical_stan(
    data, priors,
    n_iter = 300,
    n_warmup = 150,
    n_chains = 1,
    estimate_pooling = TRUE,
    verbose = FALSE
  )

  expect_equal(fit_stan_est$diagnostics$estimate_pooling, TRUE)
})

# ============================================================================
# Test 26-30: Nimble Backend (full implementation)
# ============================================================================

test_that("Nimble backend fits successfully with nimble available", {
  skip_if_not(requireNamespace("nimble", quietly = TRUE), "nimble not available")

  set.seed(111)
  data <- data.frame(
    zone = rep(c("Z1", "Z2"), each = 15),
    ilr1 = rnorm(30, 0, 1),
    ilr2 = rnorm(30, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(
    zone_names = c("Z1", "Z2"),
    component_names = c("a", "b", "c")
  )

  priors <- gc_set_hierarchy_priors(hierarchy, pooling = "moderate", verbose = FALSE)

  fit_nimble <- fit_hierarchical_nimble(
    data,
    priors,
    n_iter = 500,
    n_warmup = 250,
    n_chains = 1,
    verbose = FALSE
  )

  expect_s3_class(fit_nimble, "gc_hierarchical_fit")
  expect_equal(fit_nimble$metadata$backend, "nimble")
})

test_that("Nimble backend returns standardized structure", {
  skip_if_not(requireNamespace("nimble", quietly = TRUE), "nimble not available")

  set.seed(222)
  data <- data.frame(
    zone = rep("Z1", 20),
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(zone_names = "Z1", component_names = c("a", "b", "c"))
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  fit_nimble <- fit_hierarchical_nimble(
    data, priors,
    n_iter = 400,
    n_warmup = 200,
    n_chains = 1,
    verbose = FALSE
  )

  # Check standardized structure
  expect_true("zone_estimates" %in% names(fit_nimble))
  expect_true("global_estimates" %in% names(fit_nimble))
  expect_true("samples" %in% names(fit_nimble))
  expect_true("diagnostics" %in% names(fit_nimble))

  # Nimble should have samples
  expect_false(is.null(fit_nimble$samples))
  expect_true(is.matrix(fit_nimble$samples))
})

test_that("Nimble backend computes convergence diagnostics", {
  skip_if_not(requireNamespace("nimble", quietly = TRUE), "nimble not available")

  set.seed(333)
  data <- data.frame(
    zone = "Z1",
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(zone_names = "Z1", component_names = c("a", "b", "c"))
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  fit_nimble <- fit_hierarchical_nimble(
    data, priors,
    n_iter = 400,
    n_warmup = 200,
    n_chains = 1,
    verbose = FALSE
  )

  # Check diagnostics
  expect_true("gelman_rubin" %in% names(fit_nimble$diagnostics))
  expect_true("estimate_pooling" %in% names(fit_nimble$diagnostics))
})

test_that("MCMC backends produce posterior samples", {
  skip_if_not(requireNamespace("rstan", quietly = TRUE) || requireNamespace("nimble", quietly = TRUE),
              "Neither rstan nor nimble available")

  set.seed(444)
  data <- data.frame(
    zone = "Z1",
    ilr1 = rnorm(20, 0, 1),
    ilr2 = rnorm(20, 0, 1)
  )

  hierarchy <- gc_define_hierarchy(zone_names = "Z1", component_names = c("a", "b", "c"))
  priors <- gc_set_hierarchy_priors(hierarchy, verbose = FALSE)

  if (requireNamespace("rstan", quietly = TRUE)) {
    fit_stan <- fit_hierarchical_stan(
      data, priors,
      n_iter = 300,
      n_warmup = 150,
      n_chains = 1,
      verbose = FALSE
    )
    expect_true(nrow(fit_stan$samples) > 0)
    expect_true(ncol(fit_stan$samples) > 0)
  }

  if (requireNamespace("nimble", quietly = TRUE)) {
    fit_nimble <- fit_hierarchical_nimble(
      data, priors,
      n_iter = 300,
      n_warmup = 150,
      n_chains = 1,
      verbose = FALSE
    )
    expect_true(nrow(fit_nimble$samples) > 0)
    expect_true(ncol(fit_nimble$samples) > 0)
  }
})
