context("Per-zone diagnostics")

# Helper to create test data
create_diagnostic_test_data <- function(n_per_zone = 50) {
  set.seed(42)
  zones <- rep(c("z1", "z2", "z3"), each = n_per_zone)

  sand <- c(
    rnorm(n_per_zone, mean = 70, sd = 8),
    rnorm(n_per_zone, mean = 50, sd = 10),
    rnorm(n_per_zone, mean = 30, sd = 12)
  )
  silt <- c(
    rnorm(n_per_zone, mean = 20, sd = 5),
    rnorm(n_per_zone, mean = 35, sd = 7),
    rnorm(n_per_zone, mean = 40, sd = 8)
  )
  clay <- c(
    rnorm(n_per_zone, mean = 10, sd = 3),
    rnorm(n_per_zone, mean = 15, sd = 5),
    rnorm(n_per_zone, mean = 30, sd = 7)
  )

  # Closure
  denom <- sand + silt + clay
  sand <- (sand / denom) * 100
  silt <- (silt / denom) * 100
  clay <- (clay / denom) * 100

  data.frame(
    zone = zones,
    SAND = sand,
    SILT = silt,
    CLAY = clay
  )
}

# ============================================================================
# gc_cross_validate_per_zone tests
# ============================================================================

test_that("gc_cross_validate_per_zone requires data frame", {
  expect_error(
    gc_cross_validate_per_zone(
      data = list(a = 1),
      zone_vector = c("z1", "z2")
    ),
    "data must be a data frame"
  )
})

test_that("gc_cross_validate_per_zone validates zone vector length", {
  data <- create_diagnostic_test_data()
  expect_error(
    gc_cross_validate_per_zone(
      data = data,
      zone_vector = rep("z1", 5),
      comp_cols = c("SAND", "SILT", "CLAY")
    ),
    "zone_vector length must equal nrow"
  )
})

test_that("gc_cross_validate_per_zone returns data frame with expected columns", {
  data <- create_diagnostic_test_data()
  result <- gc_cross_validate_per_zone(
    data = data,
    zone_vector = data$zone,
    comp_cols = c("SAND", "SILT", "CLAY"),
    k = 3,
    verbose = FALSE
  )

  expect_is(result, "data.frame")
  expect_true("zone" %in% names(result))
  expect_true("n_samples" %in% names(result))
  expect_true("mean_rmse" %in% names(result))
})

test_that("gc_cross_validate_per_zone handles all zones", {
  data <- create_diagnostic_test_data()
  result <- gc_cross_validate_per_zone(
    data = data,
    zone_vector = data$zone,
    comp_cols = c("SAND", "SILT", "CLAY"),
    k = 3,
    verbose = FALSE
  )

  zones_in_result <- sort(unique(result$zone))
  zones_in_data <- sort(unique(data$zone))
  expect_equal(zones_in_result, zones_in_data)
})

test_that("gc_cross_validate_per_zone RMSE values are positive", {
  data <- create_diagnostic_test_data()
  result <- gc_cross_validate_per_zone(
    data = data,
    zone_vector = data$zone,
    comp_cols = c("SAND", "SILT", "CLAY"),
    k = 3,
    verbose = FALSE
  )

  expect_true(all(result$mean_rmse > 0, na.rm = TRUE))
})

test_that("gc_cross_validate_per_zone respects k parameter", {
  data <- create_diagnostic_test_data(n_per_zone = 30)
  result_k3 <- gc_cross_validate_per_zone(
    data = data,
    zone_vector = data$zone,
    comp_cols = c("SAND", "SILT", "CLAY"),
    k = 3,
    verbose = FALSE
  )
  result_k5 <- gc_cross_validate_per_zone(
    data = data,
    zone_vector = data$zone,
    comp_cols = c("SAND", "SILT", "CLAY"),
    k = 5,
    verbose = FALSE
  )

  # Different k should give (potentially) different results
  expect_true(nrow(result_k3) > 0)
  expect_true(nrow(result_k5) > 0)
})

test_that("gc_cross_validate_per_zone handles LOO method", {
  data <- create_diagnostic_test_data(n_per_zone = 20)
  result <- gc_cross_validate_per_zone(
    data = data,
    zone_vector = data$zone,
    comp_cols = c("SAND", "SILT", "CLAY"),
    method = "loo",
    verbose = FALSE
  )

  expect_true(nrow(result) > 0)
  expect_true(all(result$n_samples > 0))
})

# ============================================================================
# gc_compute_entropy_per_zone tests
# ============================================================================

test_that("gc_compute_entropy_per_zone requires ensemble list", {
  expect_error(
    gc_compute_entropy_per_zone(
      ensemble = list(),
      zone_col = "zone"
    ),
    "ensemble must be from gc_sim_hierarchical"
  )
})

test_that("gc_compute_entropy_per_zone requires zone specification", {
  ensemble <- list(data = data.frame(SAND = c(50, 60), SILT = c(30, 25), CLAY = c(20, 15)))

  expect_error(
    gc_compute_entropy_per_zone(
      ensemble = ensemble
    ),
    "Either zone_col or zone_definition required"
  )
})

test_that("gc_compute_entropy_per_zone returns data frame with entropy column", {
  # Create mock ensemble
  ensemble <- list(
    data = data.frame(
      SAND = rep(c(60, 50, 40), 10),
      SILT = rep(c(25, 35, 40), 10),
      CLAY = rep(c(15, 15, 20), 10),
      zone_id = rep(c("z1", "z2", "z3"), 10)
    )
  )

  result <- gc_compute_entropy_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    method = "shannon",
    verbose = FALSE
  )

  expect_is(result, "data.frame")
  expect_true("shannon_entropy" %in% names(result))
  expect_true("entropy_normalized" %in% names(result))
  expect_true("dominant_component" %in% names(result))
})

test_that("gc_compute_entropy_per_zone entropy values positive and bounded", {
  ensemble <- list(
    data = data.frame(
      SAND = rep(c(60, 50, 40), 10),
      SILT = rep(c(25, 35, 40), 10),
      CLAY = rep(c(15, 15, 20), 10),
      zone_id = rep(c("z1", "z2", "z3"), 10)
    )
  )

  result <- gc_compute_entropy_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    method = "shannon",
    normalize = TRUE,
    verbose = FALSE
  )

  expect_true(all(result$shannon_entropy >= 0))
  expect_true(all(result$entropy_normalized >= 0))
  expect_true(all(result$entropy_normalized <= 1))
})

test_that("gc_compute_entropy_per_zone handles all zones", {
  ensemble <- list(
    data = data.frame(
      SAND = rep(c(60, 50, 40), 10),
      SILT = rep(c(25, 35, 40), 10),
      CLAY = rep(c(15, 15, 20), 10),
      zone_id = rep(c("z1", "z2", "z3"), 10)
    )
  )

  result <- gc_compute_entropy_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    verbose = FALSE
  )

  zones_in_result <- sort(unique(result$zone))
  zones_in_ensemble <- sort(unique(ensemble$data$zone_id))
  expect_equal(zones_in_result, zones_in_ensemble)
})

test_that("gc_compute_entropy_per_zone classifies uncertainty", {
  ensemble <- list(
    data = data.frame(
      SAND = rep(c(60, 50, 40), 10),
      SILT = rep(c(25, 35, 40), 10),
      CLAY = rep(c(15, 15, 20), 10),
      zone_id = rep(c("z1", "z2", "z3"), 10)
    )
  )

  result <- gc_compute_entropy_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    verbose = FALSE
  )

  expect_true("uncertainty_class" %in% names(result))
  expect_true(all(result$uncertainty_class %in% c("Low", "Medium", "High")))
})

test_that("gc_compute_entropy_per_zone identifies dominant component", {
  ensemble <- list(
    data = data.frame(
      SAND = c(rep(80, 10), rep(50, 10), rep(40, 10)),
      SILT = c(rep(15, 10), rep(35, 10), rep(40, 10)),
      CLAY = c(rep(5, 10), rep(15, 10), rep(20, 10)),
      zone_id = rep(c("z1", "z2", "z3"), each = 10)
    )
  )

  result <- gc_compute_entropy_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    verbose = FALSE
  )

  # Z1 should have SAND as dominant
  z1_row <- result[result$zone == "z1", ]
  expect_equal(z1_row$dominant_component, "SAND")
})

test_that("gc_compute_entropy_per_zone handles simpson method", {
  ensemble <- list(
    data = data.frame(
      SAND = rep(c(60, 50, 40), 10),
      SILT = rep(c(25, 35, 40), 10),
      CLAY = rep(c(15, 15, 20), 10),
      zone_id = rep(c("z1", "z2", "z3"), 10)
    )
  )

  result <- gc_compute_entropy_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    method = "simpson",
    verbose = FALSE
  )

  expect_is(result, "data.frame")
  expect_true("shannon_entropy" %in% names(result))
})

# ============================================================================
# gc_bootstrap_uncertainty_per_zone tests
# ============================================================================

test_that("gc_bootstrap_uncertainty_per_zone requires data frame", {
  expect_error(
    gc_bootstrap_uncertainty_per_zone(
      data = list(a = 1),
      zone_vector = c("z1", "z2")
    ),
    "data must be a data frame"
  )
})

test_that("gc_bootstrap_uncertainty_per_zone validates zone vector length", {
  data <- create_diagnostic_test_data()
  expect_error(
    gc_bootstrap_uncertainty_per_zone(
      data = data,
      zone_vector = rep("z1", 5)
    ),
    "zone_vector length must equal nrow"
  )
})

test_that("gc_bootstrap_uncertainty_per_zone returns list with zone_estimates", {
  data <- create_diagnostic_test_data()
  result <- gc_bootstrap_uncertainty_per_zone(
    data = data,
    zone_vector = data$zone,
    n_bootstrap = 100,
    verbose = FALSE
  )

  expect_is(result, "list")
  expect_true("zone_estimates" %in% names(result))
  expect_is(result$zone_estimates, "data.frame")
})

test_that("gc_bootstrap_uncertainty_per_zone contains all zones", {
  data <- create_diagnostic_test_data()
  result <- gc_bootstrap_uncertainty_per_zone(
    data = data,
    zone_vector = data$zone,
    n_bootstrap = 100,
    verbose = FALSE
  )

  zones_in_result <- sort(unique(result$zone_estimates$zone))
  zones_in_data <- sort(unique(data$zone))
  expect_equal(zones_in_result, zones_in_data)
})

test_that("gc_bootstrap_uncertainty_per_zone CIs are sensible", {
  data <- create_diagnostic_test_data()
  result <- gc_bootstrap_uncertainty_per_zone(
    data = data,
    zone_vector = data$zone,
    n_bootstrap = 100,
    conf_level = 0.95,
    verbose = FALSE
  )

  # Lower bound < estimate < upper bound
  for (i in seq_len(nrow(result$zone_estimates))) {
    est <- result$zone_estimates$param_1_est[i]
    lower <- result$zone_estimates$param_1_lower[i]
    upper <- result$zone_estimates$param_1_upper[i]

    if (!is.na(est) && !is.na(lower) && !is.na(upper)) {
      expect_true(lower <= est)
      expect_true(est <= upper)
    }
  }
})

test_that("gc_bootstrap_uncertainty_per_zone respects confidence level", {
  data <- create_diagnostic_test_data(n_per_zone = 100)
  result_90 <- gc_bootstrap_uncertainty_per_zone(
    data = data,
    zone_vector = data$zone,
    n_bootstrap = 500,
    conf_level = 0.90,
    verbose = FALSE
  )
  result_99 <- gc_bootstrap_uncertainty_per_zone(
    data = data,
    zone_vector = data$zone,
    n_bootstrap = 500,
    conf_level = 0.99,
    verbose = FALSE
  )

  # 99% CI should be wider than 90% CI
  ci_width_90 <- mean(result_90$zone_estimates$param_1_upper - 
                       result_90$zone_estimates$param_1_lower, na.rm = TRUE)
  ci_width_99 <- mean(result_99$zone_estimates$param_1_upper - 
                       result_99$zone_estimates$param_1_lower, na.rm = TRUE)

  expect_true(ci_width_99 > ci_width_90)
})

test_that("gc_bootstrap_uncertainty_per_zone stores metadata", {
  data <- create_diagnostic_test_data()
  result <- gc_bootstrap_uncertainty_per_zone(
    data = data,
    zone_vector = data$zone,
    n_bootstrap = 250,
    conf_level = 0.95,
    verbose = FALSE
  )

  expect_equal(result$n_bootstrap, 250)
  expect_equal(result$conf_level, 0.95)
  expect_true(length(result$zones) > 0)
})

test_that("gc_bootstrap_uncertainty_per_zone handles small zones", {
  data <- create_diagnostic_test_data(n_per_zone = 10)
  result <- gc_bootstrap_uncertainty_per_zone(
    data = data,
    zone_vector = data$zone,
    n_bootstrap = 100,
    verbose = FALSE
  )

  expect_is(result, "list")
  expect_true(nrow(result$zone_estimates) > 0)
})

test_that("gc_bootstrap_uncertainty_per_zone parameter estimates reasonable", {
  data <- create_diagnostic_test_data()
  result <- gc_bootstrap_uncertainty_per_zone(
    data = data,
    zone_vector = data$zone,
    n_bootstrap = 200,
    verbose = FALSE
  )

  # Estimates should be non-zero for typical soil data
  expect_true(any(abs(result$zone_estimates$param_1_est) > 0.01, na.rm = TRUE))
})
