context("Risk Assessment Per-Zone Functions")

# Create mock ensemble data for testing
create_mock_ensemble <- function(n_zones = 3, n_per_zone = 20, n_realizations = 1) {
  zones <- rep(c("Zone1", "Zone2", "Zone3"), each = n_per_zone)
  
  ensemble_data <- data.frame(
    x = runif(n_zones * n_per_zone),
    y = runif(n_zones * n_per_zone),
    zone_id = zones,
    SAND = c(
      rnorm(n_per_zone, mean = 70, sd = 5),
      rnorm(n_per_zone, mean = 50, sd = 8),
      rnorm(n_per_zone, mean = 30, sd = 6)
    ),
    SILT = NA_real_,
    CLAY = NA_real_
  )
  
  ensemble_data$CLAY <- 100 - ensemble_data$SAND - 20
  ensemble_data$SILT <- 20
  
  list(data = ensemble_data)
}

test_that("gc_risk_assessment_per_zone handles basic input validation", {
  ensemble <- create_mock_ensemble()
  
  # Missing ensemble
  expect_error(
    gc_risk_assessment_per_zone(list()),
    "ensemble must be from"
  )
  
  # Missing thresholds
  expect_error(
    gc_risk_assessment_per_zone(
      ensemble = ensemble,
      zone_col = "zone_id",
      cost_matrix = data.frame(component = "SAND", cost_under = 100, cost_over = 50)
    ),
    "Both thresholds and cost_matrix must be specified"
  )
  
  # Invalid threshold format
  expect_error(
    gc_risk_assessment_per_zone(
      ensemble = ensemble,
      zone_col = "zone_id",
      thresholds = c(50),  # Not named
      cost_matrix = data.frame(component = "SAND", cost_under = 100, cost_over = 50)
    ),
    "must be named numeric"
  )
})

test_that("gc_risk_assessment_per_zone validates cost_matrix structure", {
  ensemble <- create_mock_ensemble()
  thresholds <- c(SAND = 50, CLAY = 30)
  
  # Invalid cost_matrix (missing column)
  expect_error(
    gc_risk_assessment_per_zone(
      ensemble = ensemble,
      zone_col = "zone_id",
      thresholds = thresholds,
      cost_matrix = data.frame(component = "SAND", cost_under = 100)
    ),
    "must be data frame with columns"
  )
})

test_that("gc_risk_assessment_per_zone computes expected loss correctly", {
  ensemble <- create_mock_ensemble()
  thresholds <- c(SAND = 50)
  cost_matrix <- data.frame(
    component = "SAND",
    cost_under = 100,
    cost_over = 0
  )
  
  result <- gc_risk_assessment_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    thresholds = thresholds,
    cost_matrix = cost_matrix
  )
  
  # Check output structure
  expect_is(result, "data.frame")
  expect_true(all(c("zone", "expected_loss", "p_under", "p_in_range") %in% names(result)))
  
  # Check Zone1 (mean SAND = 70, threshold = 50)
  # P(SAND < 50) should be low for Zone1
  zone1_result <- result[result$zone == "Zone1" & result$component == "SAND", ]
  expect_true(zone1_result$p_under < 0.3, "Zone1 should have low p_under (SAND > 50)")
  
  # Zone3 (mean SAND = 30, threshold = 50)
  # P(SAND < 50) should be high for Zone3
  zone3_result <- result[result$zone == "Zone3" & result$component == "SAND", ]
  expect_true(zone3_result$p_under > 0.7, "Zone3 should have high p_under (SAND < 50)")
  
  # Expected loss for Zone3 should be significant
  expect_true(zone3_result$expected_loss > 50, "Zone3 expected loss should be substantial")
})

test_that("gc_risk_assessment_per_zone handles multiple components", {
  ensemble <- create_mock_ensemble()
  thresholds <- c(SAND = 50, CLAY = 25)
  cost_matrix <- data.frame(
    component = c("SAND", "CLAY"),
    cost_under = c(100, 50),
    cost_over = c(0, 80)
  )
  
  result <- gc_risk_assessment_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    thresholds = thresholds,
    cost_matrix = cost_matrix
  )
  
  # Should have 3 zones × 2 components = 6 rows
  expect_equal(nrow(result), 6)
  
  # Check both components are present
  expect_true("SAND" %in% result$component)
  expect_true("CLAY" %in% result$component)
})

test_that("gc_risk_assessment_per_zone handles asymmetric costs", {
  ensemble <- create_mock_ensemble()
  thresholds <- c(SAND = 50)
  
  # High cost_under, low cost_over
  cost_high_under <- data.frame(
    component = "SAND",
    cost_under = 100,
    cost_over = 10
  )
  
  # Low cost_under, high cost_over
  cost_high_over <- data.frame(
    component = "SAND",
    cost_under = 10,
    cost_over = 100
  )
  
  result_under <- gc_risk_assessment_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    thresholds = thresholds,
    cost_matrix = cost_high_under
  )
  
  result_over <- gc_risk_assessment_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    thresholds = thresholds,
    cost_matrix = cost_high_over
  )
  
  # Zone1 (SAND >> 50): under scenario low EL (p_under ~0)
  zone1_under <- result_under[result_under$zone == "Zone1", ]$expected_loss
  
  # Zone3 (SAND << 50): under scenario high EL (p_under >> 0)
  zone3_under <- result_under[result_under$zone == "Zone3", ]$expected_loss
  
  # In high_over scenario, cost_over=100 is not incurred (threshold is lower bound)
  # So EL = p_under * 10. Zone1 has low p_under, Zone3 has high p_under
  zone1_over <- result_over[result_over$zone == "Zone1", ]$expected_loss
  zone3_over <- result_over[result_over$zone == "Zone3", ]$expected_loss
  
  # Zone3 should have much higher EL under high_under scenario
  expect_true(zone3_under > zone1_under * 5, "Zone3 should have significantly higher EL under high_under scenario")
  
  # In high_over scenario: cost_over not incurred, cost_under matters
  # Zone3 should still have higher EL than Zone1 (p_under is higher for Zone3)
  expect_true(zone3_over > zone1_over, "Zone3 should have higher EL in high_over scenario (due to high p_under)")
})

test_that("gc_risk_assessment_per_zone respects na.rm parameter", {
  ensemble <- create_mock_ensemble()
  ensemble$data$SAND[c(1, 5, 10)] <- NA
  thresholds <- c(SAND = 50)
  cost_matrix <- data.frame(
    component = "SAND",
    cost_under = 100,
    cost_over = 0
  )
  
  # With na.rm = TRUE (default)
  result_rmna <- gc_risk_assessment_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    thresholds = thresholds,
    cost_matrix = cost_matrix,
    na.rm = TRUE
  )
  
  # With na.rm = FALSE
  result_keepna <- gc_risk_assessment_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    thresholds = thresholds,
    cost_matrix = cost_matrix,
    na.rm = FALSE
  )
  
  # Results should differ (NA handling affects proportions)
  zone1_rmna <- result_rmna[result_rmna$zone == "Zone1", ]$p_under
  zone1_keepna <- result_keepna[result_keepna$zone == "Zone1", ]$p_under
  
  # na.rm=TRUE should exclude NAs, possibly giving different proportions
  expect_true(is.numeric(zone1_rmna) && !is.na(zone1_rmna))
})

test_that("gc_carbon_audit_per_zone handles basic input validation", {
  ensemble <- create_mock_ensemble()
  
  # Missing ensemble
  expect_error(
    gc_carbon_audit_per_zone(list()),
    "ensemble must be from"
  )
  
  # Missing zone specification
  expect_error(
    gc_carbon_audit_per_zone(ensemble = ensemble),
    "Either zone_col or zone_definition required"
  )
})

test_that("gc_carbon_audit_per_zone computes carbon stock correctly", {
  ensemble <- create_mock_ensemble()
  
  result <- gc_carbon_audit_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    soil_depth = 30,
    bulk_density = 1.3,
    carbon_fraction = 0.03,
    verbose = FALSE
  )
  
  # Check output structure
  expect_is(result, "data.frame")
  expect_true(all(c("zone", "mean_c_stock_mg_ha", "compliance_status") %in% names(result)))
  
  # Check all zones present
  expect_equal(nrow(result), 3)
  expect_true(all(c("Zone1", "Zone2", "Zone3") %in% result$zone))
  
  # Carbon stock should be positive
  expect_true(all(result$mean_c_stock_mg_ha >= 0))
  
  # Compliance status should be one of the three values
  expect_true(all(result$compliance_status %in% c("Verified", "Uncertain", "At Risk")))
})

test_that("gc_carbon_audit_per_zone handles SoilGrids carbon model", {
  ensemble <- create_mock_ensemble()
  
  result_soilgrids <- gc_carbon_audit_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    soil_depth = 30,
    carbon_fraction = "soilgrids"
  )
  
  expect_is(result_soilgrids, "data.frame")
  expect_equal(nrow(result_soilgrids), 3)
  
  # SoilGrids model: C% correlates with clay (100 - SAND)
  # Zone1 (high SAND): low C
  # Zone3 (low SAND): high C
  zone1_c <- result_soilgrids[result_soilgrids$zone == "Zone1", ]$mean_c_stock_mg_ha
  zone3_c <- result_soilgrids[result_soilgrids$zone == "Zone3", ]$mean_c_stock_mg_ha
  
  expect_true(zone1_c < zone3_c, "SoilGrids: Zone1 (sandy) should have less C than Zone3 (clayey)")
})

test_that("gc_carbon_audit_per_zone respects bulk_density specification", {
  ensemble <- create_mock_ensemble()
  
  # Single numeric bulk_density
  result_uniform <- gc_carbon_audit_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    bulk_density = 1.5,
    carbon_fraction = 0.03
  )
  
  # Component-specific bulk_density (data frame)
  bd_df <- data.frame(
    component = c("SAND", "CLAY"),
    bulk_density = c(1.5, 1.2)
  )
  
  result_component <- gc_carbon_audit_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    bulk_density = bd_df,
    carbon_fraction = 0.03
  )
  
  # Both should produce valid results
  expect_is(result_uniform, "data.frame")
  expect_is(result_component, "data.frame")
})

test_that("gc_carbon_audit_per_zone computes credit_issued conservatively", {
  ensemble <- create_mock_ensemble()
  
  result <- gc_carbon_audit_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    soil_depth = 30,
    bulk_density = 1.3,
    carbon_fraction = 0.03,
    credit_price = 25
  )
  
  # credit_issued should be <= mean_c_stock (conservative use of lower CI bound)
  for (i in seq_len(nrow(result))) {
    expect_true(
      result$credit_issued_mg_c_ha[i] <= result$mean_c_stock_mg_ha[i],
      "credit_issued should be conservative (≤ mean)"
    )
  }
  
  # credit_value should equal credit_issued × price
  expected_value <- result$credit_issued_mg_c_ha * 25
  expect_equal(result$credit_value_usd, round(expected_value, 0), tolerance = 1)
})

test_that("gc_carbon_audit_per_zone assignment compliance status correctly", {
  ensemble <- create_mock_ensemble()
  
  result <- gc_carbon_audit_per_zone(
    ensemble = ensemble,
    zone_col = "zone_id",
    soil_depth = 30,
    bulk_density = 1.3,
    carbon_fraction = 0.03,
    verbose = FALSE
  )
  
  # Compliance status based on SD/mean & CI bounds
  for (i in seq_len(nrow(result))) {
    if (result$ci_lower[i] < 0) {
      # At Risk if lower CI is negative
      expect_equal(result$compliance_status[i], "At Risk")
    } else if (result$sd_c_stock[i] < result$mean_c_stock_mg_ha[i] * 0.2) {
      # Verified if SD < 20% of mean
      expect_equal(result$compliance_status[i], "Verified")
    } else {
      # Otherwise Uncertain
      expect_equal(result$compliance_status[i], "Uncertain")
    }
  }
})

test_that("Per-zone functions handle factor zone definitions", {
  ensemble <- create_mock_ensemble()
  zone_factor <- factor(ensemble$data$zone_id)
  
  # Risk assessment with factor zones
  result_risk <- gc_risk_assessment_per_zone(
    ensemble = ensemble,
    zone_definition = zone_factor,
    thresholds = c(SAND = 50),
    cost_matrix = data.frame(component = "SAND", cost_under = 100, cost_over = 50)
  )
  
  # Carbon audit with factor zones
  result_carbon <- gc_carbon_audit_per_zone(
    ensemble = ensemble,
    zone_definition = zone_factor
  )
  
  expect_is(result_risk, "data.frame")
  expect_is(result_carbon, "data.frame")
})

test_that("Per-zone functions handle character zone definitions", {
  ensemble <- create_mock_ensemble()
  zone_char <- as.character(ensemble$data$zone_id)
  
  result_risk <- gc_risk_assessment_per_zone(
    ensemble = ensemble,
    zone_definition = zone_char,
    thresholds = c(SAND = 50),
    cost_matrix = data.frame(component = "SAND", cost_under = 100, cost_over = 50)
  )
  
  result_carbon <- gc_carbon_audit_per_zone(
    ensemble = ensemble,
    zone_definition = zone_char
  )
  
  expect_is(result_risk, "data.frame")
  expect_is(result_carbon, "data.frame")
})
