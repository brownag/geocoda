context("Risk Assessment Functions - Probability & Percentile Mapping")

# Setup: Create synthetic rasters for testing
create_gaussian_raster <- function(n_realization = 100, n_cell_x = 10, n_cell_y = 10,
                                   mean = 0, sd = 1) {
  set.seed(42)
  rast <- terra::rast(nrow = n_cell_y, ncol = n_cell_x, nlyr = n_realization)

  # Fill with Gaussian random values
  for (i in 1:n_realization) {
    terra::values(rast[[i]]) <- rnorm(n_cell_x * n_cell_y, mean = mean, sd = sd)
  }

  return(rast)
}

# Test: gc_probability_map basic functionality
test_that("gc_probability_map returns SpatRaster with correct structure", {
  sims <- create_gaussian_raster(n_realization = 50, mean = 0, sd = 1)
  prob <- gc_probability_map(sims, threshold = 0, operator = ">")

  expect_s4_class(prob, "SpatRaster")
  expect_equal(terra::nlyr(prob), 1)
  expect_equal(terra::ncol(prob), 10)
  expect_equal(terra::nrow(prob), 10)
})

# Test: gc_probability_map values in [0, 1]
test_that("gc_probability_map returns values in [0, 1]", {
  sims <- create_gaussian_raster(n_realization = 100)
  prob <- gc_probability_map(sims, threshold = 0)

  values <- terra::values(prob)
  expect_true(all(values >= 0 & values <= 1, na.rm = TRUE))
})

# Test: gc_probability_map with Gaussian data matches theory
test_that("gc_probability_map with Gaussian data matches theoretical probability", {
  # Generate large number of realizations from N(0, 1)
  sims <- create_gaussian_raster(n_realization = 1000, n_cell_x = 50, n_cell_y = 50,
                                 mean = 0, sd = 1)

  # P(Z > 0) should be approximately 0.5 for N(0,1)
  prob_0 <- gc_probability_map(sims, threshold = 0, operator = ">")
  mean_prob_0 <- mean(terra::values(prob_0), na.rm = TRUE)
  expect_true(abs(mean_prob_0 - 0.5) < 0.05)  # Within 5%

  # P(Z > 1) should be approximately 0.16 for N(0,1)
  prob_1 <- gc_probability_map(sims, threshold = 1, operator = ">")
  mean_prob_1 <- mean(terra::values(prob_1), na.rm = TRUE)
  expect_true(abs(mean_prob_1 - 0.16) < 0.05)
})

# Test: gc_probability_map operators
test_that("gc_probability_map operators work correctly", {
  sims <- create_gaussian_raster(n_realization = 100)

  prob_gt <- gc_probability_map(sims, threshold = 0, operator = ">")
  prob_gte <- gc_probability_map(sims, threshold = 0, operator = ">=")
  prob_lt <- gc_probability_map(sims, threshold = 0, operator = "<")
  prob_lte <- gc_probability_map(sims, threshold = 0, operator = "<=")

  # Complementary probabilities should sum to ~1
  prob_sum_gt_lte <- terra::values(prob_gt) + terra::values(prob_lte)
  expect_true(all(abs(prob_sum_gt_lte - 1.0) < 0.02, na.rm = TRUE))
})

# Test: gc_probability_map NA handling
test_that("gc_probability_map handles NA values correctly", {
  sims <- create_gaussian_raster(n_realization = 50)

  # Add some NAs
  for (i in 1:5) {
    values <- terra::values(sims[[i]])
    values[1:5] <- NA
    terra::values(sims[[i]]) <- values
  }

  # With na.rm = TRUE
  prob_rm <- gc_probability_map(sims, threshold = 0, na.rm = TRUE)
  expect_true(!all(is.na(terra::values(prob_rm))))

  # Values should still be in [0, 1]
  values_rm <- terra::values(prob_rm)
  expect_true(all(values_rm >= 0 & values_rm <= 1, na.rm = TRUE))
})

# Test: gc_percentile_map basic functionality
test_that("gc_percentile_map returns correct number of layers", {
  sims <- create_gaussian_raster(n_realization = 100)
  perc <- gc_percentile_map(sims, percentiles = c(0.1, 0.5, 0.9))

  expect_s4_class(perc, "SpatRaster")
  expect_equal(terra::nlyr(perc), 3)
  expect_equal(names(perc), c("P0.1", "P0.5", "P0.9"))
})

# Test: gc_percentile_map ordering
test_that("gc_percentile_map percentiles are in correct order", {
  sims <- create_gaussian_raster(n_realization = 100)
  perc <- gc_percentile_map(sims, percentiles = c(0.1, 0.5, 0.9))

  p10_values <- terra::values(perc[[1]])
  p50_values <- terra::values(perc[[2]])
  p90_values <- terra::values(perc[[3]])

  # At each cell: p10 <= p50 <= p90
  expect_true(all(p10_values <= p50_values, na.rm = TRUE))
  expect_true(all(p50_values <= p90_values, na.rm = TRUE))
})

# Test: gc_percentile_map custom names
test_that("gc_percentile_map accepts custom names_format", {
  sims <- create_gaussian_raster(n_realization = 100)
  perc <- gc_percentile_map(
    sims,
    percentiles = c(0.25, 0.5, 0.75),
    names_format = "Q{p*100}"
  )

  expect_equal(names(perc), c("Q25", "Q50", "Q75"))
})

# Test: gc_percentile_map with Gaussian data
test_that("gc_percentile_map matches empirical quantiles", {
  sims <- create_gaussian_raster(n_realization = 1000)
  perc <- gc_percentile_map(sims, percentiles = c(0.25, 0.5, 0.75))

  perc_values <- terra::as.data.frame(perc)
  sims_values <- terra::as.data.frame(sims)

  # For each cell, check percentile calculation
  for (cell in 1:nrow(sims_values)) {
    cell_sims <- as.numeric(sims_values[cell, ])
    cell_perc <- as.numeric(perc_values[cell, ])

    # Compute empirical quantiles
    empirical_q <- stats::quantile(cell_sims, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

    # Should be very close
    expect_equal(cell_perc, as.numeric(empirical_q), tolerance = 0.1)
  }
})

# Test: gc_risk_assessment basic functionality
test_that("gc_risk_assessment returns list with expected structure", {
  sims <- create_gaussian_raster(n_realization = 50)
  risk <- gc_risk_assessment(
    sims,
    threshold = 0,
    loss_function = "linear",
    cost_params = list(a = 1, b = 1)
  )

  expect_is(risk, "list")
  expect_true(all(c("expected_loss", "decision", "cost_params", "metadata") %in% names(risk)))
  expect_s4_class(risk$expected_loss, "SpatRaster")
  expect_s4_class(risk$decision, "SpatRaster")
})

# Test: gc_risk_assessment symmetric loss
test_that("gc_risk_assessment symmetric loss is symmetric around threshold", {
  sims <- create_gaussian_raster(n_realization = 1000, mean = 50, sd = 10)

  loss_below <- gc_risk_assessment(
    sims,
    threshold = 50,
    loss_function = "linear",
    cost_params = list(a = 1, b = 1),
    operator = "<"
  )

  loss_above <- gc_risk_assessment(
    sims,
    threshold = 50,
    loss_function = "linear",
    cost_params = list(a = 1, b = 1),
    operator = ">"
  )

  # Expected losses should be similar (symmetric)
  el_below <- mean(terra::values(loss_below$expected_loss), na.rm = TRUE)
  el_above <- mean(terra::values(loss_above$expected_loss), na.rm = TRUE)

  expect_true(abs(el_below - el_above) / mean(el_below, el_above) < 0.1)  # Within 10%
})

# Test: gc_risk_assessment asymmetric loss
test_that("gc_risk_assessment asymmetric loss reflects cost structure", {
  sims <- create_gaussian_raster(n_realization = 500, mean = 50, sd = 10)

  # Higher cost for under-estimate
  loss_asym <- gc_risk_assessment(
    sims,
    threshold = 50,
    loss_function = "asymmetric",
    cost_params = list(a = 2, b = 1),
    operator = "<"
  )

  el <- mean(terra::values(loss_asym$expected_loss), na.rm = TRUE)
  expect_true(el > 0)  # Non-zero loss
})

# Test: gc_carbon_audit basic functionality
test_that("gc_carbon_audit returns expected output structure", {
  sims <- create_gaussian_raster(n_realization = 100, mean = 55, sd = 15)
  audit <- gc_carbon_audit(
    sims,
    credit_threshold = 50,
    confidence = 0.9,
    audit_type = "conservative"
  )

  expect_is(audit, "list")
  expect_true(all(c("probability_map", "conservative_estimate", "issued_credit",
                     "summary", "audit_report", "metadata") %in% names(audit)))
  expect_is(audit$audit_report, "character")
  expect_true(nchar(audit$audit_report) > 100)
})

# Test: gc_carbon_audit audit types ordering
test_that("gc_carbon_audit estimates follow: conservative < expected < optimistic", {
  sims <- create_gaussian_raster(n_realization = 500, mean = 60, sd = 20)

  audit_cons <- gc_carbon_audit(sims, credit_threshold = 50, confidence = 0.9,
                                 audit_type = "conservative")
  audit_exp <- gc_carbon_audit(sims, credit_threshold = 50, confidence = 0.9,
                                audit_type = "expected")
  audit_opt <- gc_carbon_audit(sims, credit_threshold = 50, confidence = 0.9,
                                audit_type = "optimistic")

  # Extract mean estimates
  cons_mean <- mean(terra::values(audit_cons$conservative_estimate), na.rm = TRUE)
  exp_mean <- mean(terra::values(audit_exp$conservative_estimate), na.rm = TRUE)
  opt_mean <- mean(terra::values(audit_opt$conservative_estimate), na.rm = TRUE)

  expect_true(cons_mean < exp_mean)
  expect_true(exp_mean < opt_mean)
})

# Test: gc_carbon_audit generates report
test_that("gc_carbon_audit generates human-readable report", {
  sims <- create_gaussian_raster(n_realization = 100)
  audit <- gc_carbon_audit(sims, credit_threshold = 50)

  report <- audit$audit_report

  expect_is(report, "character")
  expect_true(grepl("Credit Threshold", report))
  expect_true(grepl("Confidence Level", report))
  expect_true(grepl("Audit Type", report))
  expect_true(nchar(report) > 200)
})

# Test: Input validation
test_that("Functions validate inputs appropriately", {
  sims <- create_gaussian_raster()

  # Invalid operator
  expect_error(
    gc_probability_map(sims, threshold = 0, operator = "INVALID"),
    "operator must be"
  )

  # Invalid loss function
  expect_error(
    gc_risk_assessment(sims, threshold = 0, loss_function = "INVALID"),
    "Unknown loss function"
  )

  # Invalid threshold
  expect_error(
    gc_carbon_audit(sims, credit_threshold = "not_numeric"),
    "must be a single numeric"
  )
})

# Test: Edge cases
test_that("Functions handle edge cases gracefully", {
  sims <- create_gaussian_raster(n_realization = 10)

  # Very low confidence
  audit_low <- gc_carbon_audit(sims, credit_threshold = 50, confidence = 0.05)
  expect_is(audit_low, "list")

  # Threshold outside data range
  prob_extreme <- gc_probability_map(sims, threshold = 1000)
  prob_values <- terra::values(prob_extreme)
  expect_true(all(prob_values == 0 | is.na(prob_values)))

  # Step loss function
  risk_step <- gc_risk_assessment(
    sims,
    threshold = 0,
    loss_function = "step",
    cost_params = list(a = 1, b = 2)
  )
  expect_s4_class(risk_step$expected_loss, "SpatRaster")
})

# Integration test with realistic data
test_that("Functions work with realistic carbon stock data", {
  # Simulate realistic SOC data: Mean ~50 Mg/ha, SD ~15
  carbon_sims <- create_gaussian_raster(
    n_realization = 200,
    n_cell_x = 20,
    n_cell_y = 20,
    mean = 50,
    sd = 15
  )

  # Full workflow
  prob_map <- gc_probability_map(carbon_sims, threshold = 50)
  perc_map <- gc_percentile_map(carbon_sims, percentiles = c(0.1, 0.5, 0.9))
  audit <- gc_carbon_audit(carbon_sims, credit_threshold = 50, confidence = 0.9)

  # Verify outputs
  expect_s4_class(prob_map, "SpatRaster")
  expect_s4_class(perc_map, "SpatRaster")
  expect_is(audit, "list")

  # P(SOC > 50) should be close to 0.5 for mean=50 data
  mean_prob <- mean(terra::values(prob_map), na.rm = TRUE)
  expect_true(abs(mean_prob - 0.5) < 0.15)
})

# Test: Multiple operators with complementary logic
test_that("Complementary operators sum to 1", {
  sims <- create_gaussian_raster(n_realization = 100)
  threshold <- 0

  prob_gt <- gc_probability_map(sims, threshold = threshold, operator = ">")
  prob_lte <- gc_probability_map(sims, threshold = threshold, operator = "<=")

  sum_probs <- terra::values(prob_gt) + terra::values(prob_lte)

  # Should sum to ~1 (allowing for small numerical errors)
  expect_true(all(abs(sum_probs - 1.0) < 0.01, na.rm = TRUE))
})
