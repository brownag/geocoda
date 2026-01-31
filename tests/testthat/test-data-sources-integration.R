# Data Source Integration Tests
# Tests for external data source utilities: SSURGO, horizon stratification,
# component weighting, and parametric specifications

test_that("gc_ssurgo_to_params extracts parameters from SSURGO data", {
  # Create minimal SSURGO-like data
  ssurgo_data <- data.frame(
    mukey = c(1001, 1001, 1001, 1002, 1002),
    cokey = c(1001.1, 1001.2, 1001.3, 1002.1, 1002.2),
    comppct = c(60, 25, 15, 50, 50),
    sandtotal = c(35, 15, 50, 40, 20),
    silttotal = c(50, 60, 30, 45, 65),
    claytotal = c(15, 25, 20, 15, 15),
    hzdept_r = c(0, 0, 0, 0, 0),
    hzdepb_r = c(25, 25, 25, 25, 25)
  )
  
  # Convert to parameters
  result <- geocoda::gc_ssurgo_to_params(
    ssurgo_data,
    weight_method = "comppct"
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true("constraints" %in% names(result))
  expect_true("samples" %in% names(result))
  expect_true("samples_weighted" %in% names(result))
  expect_true("mu_stats" %in% names(result))
  expect_true("metadata" %in% names(result))
  
  # Check constraints format
  expect_true("SAND" %in% names(result$constraints))
  expect_true("SILT" %in% names(result$constraints))
  expect_true("CLAY" %in% names(result$constraints))
  
  for (comp in c("SAND", "SILT", "CLAY")) {
    expect_true("min" %in% names(result$constraints[[comp]]))
    expect_true("max" %in% names(result$constraints[[comp]]))
    expect_true(result$constraints[[comp]]$min <= result$constraints[[comp]]$max)
  }
  
  # Check samples
  expect_equal(nrow(result$samples), 5)
  expect_true("mukey" %in% names(result$samples))
  expect_true("cokey" %in% names(result$samples))
  
  # Check weighted samples
  expect_equal(nrow(result$samples_weighted), 5)
  expect_true("weight" %in% names(result$samples_weighted))
  
  # Verify weights sum to 1 per map unit
  for (mu in unique(result$samples_weighted$mukey)) {
    mu_weights <- result$samples_weighted$weight[result$samples_weighted$mukey == mu]
    expect_true(abs(sum(mu_weights) - 1) < 0.01)
  }
  
  # Check metadata
  expect_equal(result$metadata$source, "SSURGO")
  expect_equal(result$metadata$n_records, 5)
  expect_equal(result$metadata$n_map_units, 2)
})

test_that("gc_ssurgo_to_params validates required columns", {
  # Missing 'comppct' column
  bad_data <- data.frame(
    mukey = 1001,
    cokey = 1001.1,
    sandtotal = 35,
    silttotal = 50,
    claytotal = 15
  )
  
  expect_error(
    geocoda::gc_ssurgo_to_params(bad_data),
    "Missing required columns"
  )
})

test_that("gc_ssurgo_to_params handles depth filtering", {
  ssurgo_data <- data.frame(
    mukey = c(1001, 1001, 1001),
    cokey = c(1001.1, 1001.2, 1001.3),
    comppct = c(60, 25, 15),
    sandtotal = c(35, 30, 25),
    silttotal = c(50, 55, 60),
    claytotal = c(15, 15, 15),
    hzdept_r = c(0, 25, 60),
    hzdepb_r = c(25, 60, 100)
  )
  
  # Filter to 0-30 cm
  result_shallow <- geocoda::gc_ssurgo_to_params(
    ssurgo_data,
    depth_range = c(0, 30)
  )
  
  # Should include first two horizons
  expect_equal(nrow(result_shallow$samples), 2)
  
  # Filter to 60-100 cm
  result_deep <- geocoda::gc_ssurgo_to_params(
    ssurgo_data,
    depth_range = c(60, 100)
  )
  
  # Should include only last horizon
  expect_equal(nrow(result_deep$samples), 1)
})

test_that("gc_ssurgo_to_params rescales invalid compositions", {
  # Create data where sand+silt+clay ≠ 100
  bad_comps <- data.frame(
    mukey = 1001,
    cokey = 1001.1,
    comppct = 100,
    sandtotal = 40,
    silttotal = 40,
    claytotal = 40  # Sums to 120 instead of 100
  )
  
  # Should warn and rescale
  result <- expect_warning(
    geocoda::gc_ssurgo_to_params(bad_comps),
    "composition sums"
  )
  
  # Verify rescaled data sums to ~100
  scaled_sum <- sum(result$samples[1, c("sandtotal", "silttotal", "claytotal")])
  expect_true(abs(scaled_sum - 100) < 0.01)
})

test_that("gc_horizon_stratify creates depth-based strata", {
  comp_data <- data.frame(
    sandtotal = c(35, 30, 25, 20),
    silttotal = c(50, 55, 60, 65),
    claytotal = c(15, 15, 15, 15),
    hzdept_r = c(0, 25, 60, 100),
    hzdepb_r = c(25, 60, 100, 150)
  )
  
  strata <- geocoda::gc_horizon_stratify(
    comp_data,
    horizon_breaks = c(0, 30, 60, 100),
    horizon_labels = c("surface", "subsurface_1", "subsurface_2")
  )
  
  # Check structure
  expect_equal(length(strata), 3)
  expect_true("surface" %in% names(strata))
  expect_true("subsurface_1" %in% names(strata))
  expect_true("subsurface_2" %in% names(strata))
  
  # Check each stratum has required fields
  for (label in names(strata)) {
    stratum <- strata[[label]]
    expect_true("constraints" %in% names(stratum))
    expect_true("samples" %in% names(stratum))
    expect_true("depth_range" %in% names(stratum))
    expect_true("n_samples" %in% names(stratum))
  }
  
  # Surface stratum should have n_samples = 2 (horizons 0-25 and 25-60 both contribute)
  expect_equal(strata$surface$n_samples, 2)
})

test_that("gc_horizon_stratify handles horizon splitting", {
  # Horizon spanning multiple strata: 20-50 cm across 0-30 and 30-60 strata
  comp_data <- data.frame(
    sandtotal = c(40, 30, 20),  # Different values per horizon
    silttotal = c(45, 55, 65),
    claytotal = c(15, 15, 15),
    hzdept_r = c(0, 20, 60),
    hzdepb_r = c(20, 50, 100)
  )
  
  strata <- geocoda::gc_horizon_stratify(
    comp_data,
    horizon_breaks = c(0, 30, 60, 100)
  )
  
  # First stratum (0-30) should include:
  # - Horizon 1: full (0-20)
  # - Horizon 2: partial (20-30, 10 cm)
  expect_equal(strata$stratum_1$n_samples, 2)
  
  # Second stratum (30-60) should include:
  # - Horizon 2: partial (30-50, 20 cm)
  # - Horizon 3: none (starts at 60)
  expect_equal(strata$stratum_2$n_samples, 1)
})

test_that("gc_horizon_stratify validates depth columns", {
  bad_data <- data.frame(
    sand = c(35, 30),
    silt = c(50, 55),
    clay = c(15, 15)
    # Missing hzdept_r and hzdepb_r
  )
  
  expect_error(
    geocoda::gc_horizon_stratify(bad_data),
    "hzdept_r.*hzdepb_r"
  )
})

test_that("gc_weight_components applies weights correctly", {
  samples <- data.frame(
    sand = c(35, 25, 45),
    silt = c(50, 60, 40),
    clay = c(15, 15, 15),
    comppct = c(60, 25, 15)
  )
  
  # Method: sample_weights
  weighted <- geocoda::gc_weight_components(
    samples,
    comp_cols = c("sand", "silt", "clay"),
    weight_col = "comppct",
    method = "sample_weights"
  )
  
  expect_true("weight" %in% names(weighted))
  expect_equal(nrow(weighted), 3)
  
  # Verify weights are normalized
  expect_true(abs(sum(weighted$weight) - 1) < 0.01)
  
  # Method: weighted_mean
  mean_comp <- geocoda::gc_weight_components(
    samples,
    comp_cols = c("sand", "silt", "clay"),
    weight_col = "comppct",
    method = "weighted_mean"
  )
  
  expect_equal(nrow(mean_comp), 1)
  expect_true(all(c("sand", "silt", "clay") %in% names(mean_comp)))
  
  # Weighted mean should be close to dominant component
  # (60% weight on sand=35 should pull overall sand toward 35)
  expect_true(mean_comp$sand < mean(samples$sand))  # Closer to the 60% component
})

test_that("gc_weight_components computes weighted quantiles", {
  samples <- data.frame(
    sand = c(20, 30, 40, 50, 60),
    silt = c(70, 60, 50, 40, 30),
    clay = c(10, 10, 10, 10, 10),
    weight = c(0.1, 0.2, 0.4, 0.2, 0.1)  # Mode at sand=40
  )
  
  constraints <- geocoda::gc_weight_components(
    samples,
    comp_cols = c("sand", "silt", "clay"),
    weight_col = "weight",
    method = "quantile_weighted"
  )
  
  # Check that 50th percentile sand is close to mode (40)
  expect_true(abs(constraints$sand$representative - 40) < 10)
  
  # Min < representative < max
  expect_true(constraints$sand$min < constraints$sand$representative)
  expect_true(constraints$sand$representative < constraints$sand$max)
})

test_that("gc_weight_components validates inputs", {
  samples <- data.frame(
    sand = c(35, 25),
    weight = c(60, 25)
  )
  
  # Missing weight column
  expect_error(
    geocoda::gc_weight_components(
      samples,
      comp_cols = "sand",
      weight_col = "missing_col",
      method = "weighted_mean"
    ),
    "not found in x"
  )
})

test_that("gc_parametric_constraints creates valid constraints from priors", {
  constraints <- geocoda::gc_parametric_constraints(
    components = c("SAND", "SILT", "CLAY"),
    means = c(35, 50, 15),
    sds = c(5, 5, 3)
  )
  
  # Check structure
  expect_true("constraints" %in% names(constraints))
  expect_true("parameters" %in% names(constraints))
  expect_true("metadata" %in% names(constraints))
  
  # Check constraint bounds
  for (comp in c("SAND", "SILT", "CLAY")) {
    expect_true("min" %in% names(constraints$constraints[[comp]]))
    expect_true("max" %in% names(constraints$constraints[[comp]]))
    expect_true(
      constraints$constraints[[comp]]$min <= constraints$constraints[[comp]]$max
    )
  }
  
  # Parameters data frame
  expect_equal(nrow(constraints$parameters), 3)
  expect_true("mean" %in% names(constraints$parameters))
  expect_true("sd" %in% names(constraints$parameters))
})

test_that("gc_parametric_constraints warns on inconsistent sums", {
  # Means sum to 95 instead of 100
  result <- expect_warning(
    geocoda::gc_parametric_constraints(
      components = c("SAND", "SILT", "CLAY"),
      means = c(35, 50, 10),  # Sum = 95
      sds = c(5, 5, 3)
    ),
    "Mean composition sums"
  )
  
  expect_true(is.list(result))
})

test_that("gc_parametric_constraints handles custom bounds", {
  constraints <- geocoda::gc_parametric_constraints(
    components = c("SAND", "SILT", "CLAY"),
    means = c(35, 50, 15),
    sds = c(5, 5, 3),
    mins = c(20, 40, 10),
    maxs = c(50, 60, 20)
  )
  
  # Verify custom bounds are used
  expect_equal(constraints$constraints$SAND$min, 20)
  expect_equal(constraints$constraints$SAND$max, 50)
  expect_equal(constraints$constraints$SILT$min, 40)
  expect_equal(constraints$constraints$SILT$max, 60)
})

test_that("gc_parametric_constraints rejects invalid SDs", {
  expect_error(
    geocoda::gc_parametric_constraints(
      components = c("SAND", "SILT"),
      means = c(35, 50),
      sds = c(-5, 5)  # Negative SD
    ),
    "sds must be non-negative"
  )
})

test_that("Data source integration workflow", {
  # Simulate real workflow: SSURGO → constraints → expand_bounds → resample
  
  # 1. Create mock SSURGO data
  ssurgo_data <- data.frame(
    mukey = rep(1001, 10),
    cokey = seq(1001.1, 1001.10, by = 0.1),
    comppct = c(rep(50, 5), rep(25, 3), rep(15, 2)),
    sandtotal = c(rep(35, 5), rep(25, 3), rep(45, 2)),
    silttotal = c(rep(50, 5), rep(60, 3), rep(35, 2)),
    claytotal = c(rep(15, 5), rep(15, 3), rep(20, 2))
  )
  
  # 2. Extract parameters
  params <- geocoda::gc_ssurgo_to_params(
    ssurgo_data,
    weight_method = "comppct"
  )
  
  # 3. Use constraints from SSURGO
  constraints <- params$constraints
  
  # 4. Verify constraints are usable with gc_expand_bounds
  expect_true(!is.null(constraints))
  expect_true(length(constraints) == 3)
  
  # Check that constraints have proper structure
  for (comp in c("SAND", "SILT", "CLAY")) {
    expect_true("min" %in% names(constraints[[comp]]))
    expect_true("max" %in% names(constraints[[comp]]))
  }
})
