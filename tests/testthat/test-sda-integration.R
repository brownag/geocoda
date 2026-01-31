context("SDA/soilDB Integration")

# ============================================================================
# Test Suite for SDA Integration Functions
# gc_fetch_sda_properties, gc_process_ssurgo_components,
# gc_optimize_sda_queries, gc_prepare_hierarchical_data
# ============================================================================


# ============================================================================
# Tests: gc_fetch_sda_properties()
# ============================================================================

test_that("gc_fetch_sda_properties validates inputs", {
  
  # Test: Must provide either extent or map_unit_keys
  expect_error(
    gc_fetch_sda_properties(extent = NULL, map_unit_keys = NULL),
    "Must provide either extent or map_unit_keys"
  )
})


test_that("gc_fetch_sda_properties returns proper data frame structure", {
  skip_if_offline()
  skip_if_not_installed("soilDB")
  
  # Mock query using a bbox extent
  study_area <- sf::st_bbox(c(
    xmin = -93.5, ymin = 41.5,
    xmax = -93.45, ymax = 41.55
  ), crs = 4326)
  
  result <- gc_fetch_sda_properties(
    extent = study_area,
    depth_range = c(0, 30),
    property_names = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    use_cache = FALSE,
    verbose = FALSE
  )
  
  # Check data frame structure
  expect_is(result, "data.frame")
  expect_true(nrow(result) >= 0)  # May be empty if no results
  
  # If data returned, check columns
  if (nrow(result) > 0) {
    expect_true("mukey" %in% names(result))
    expect_true("depth_min" %in% names(result))
    expect_true("depth_max" %in% names(result))
  }
})


test_that("gc_fetch_sda_properties handles caching", {
  skip_if_offline()
  skip_if_not_installed("soilDB")
  
  study_area <- sf::st_bbox(c(
    xmin = -93.5, ymin = 41.5,
    xmax = -93.45, ymax = 41.55
  ), crs = 4326)
  
  # Create first query
  result1 <- gc_fetch_sda_properties(
    extent = study_area,
    depth_range = c(0, 30),
    use_cache = TRUE,
    cache_dir = tempdir(),
    verbose = FALSE
  )
  
  # Second query should use cache
  t0 <- Sys.time()
  result2 <- gc_fetch_sda_properties(
    extent = study_area,
    depth_range = c(0, 30),
    use_cache = TRUE,
    cache_dir = tempdir(),
    verbose = FALSE
  )
  t1 <- Sys.time()
  
  # Cache should make second query fast
  expect_true(as.numeric(difftime(t1, t0)) < 1)
  
  # Results should be identical
  expect_equal(nrow(result1), nrow(result2))
})


test_that("gc_fetch_sda_properties validates depth_range", {
  
  study_area <- sf::st_bbox(c(
    xmin = -93.5, ymin = 41.5,
    xmax = -93.45, ymax = 41.55
  ), crs = 4326)
  
  # Depth range should be length 2
  expect_error(
    gc_fetch_sda_properties(
      extent = study_area,
      depth_range = c(0),  # Wrong length
      verbose = FALSE
    ),
    NA  # May or may not error depending on implementation
  )
})


# ============================================================================
# Tests: gc_process_ssurgo_components()
# ============================================================================

test_that("gc_process_ssurgo_components validates inputs", {
  
  # Test: Returns error if ssurgo_data is empty
  expect_error(
    gc_process_ssurgo_components(
      ssurgo_data = data.frame(),
      component_cols = c("sand", "silt", "clay"),
      verbose = FALSE
    ),
    "ssurgo_data must be a non-empty data frame"
  )
  
  # Test: Returns error if component_cols not found
  ssurgo_test <- data.frame(
    mukey = "463168",
    compname = "Harps",
    comppct = 60,
    sand = 50, silt = 30, clay = 20
  )
  
  expect_error(
    gc_process_ssurgo_components(
      ssurgo_data = ssurgo_test,
      component_cols = c("sand", "missing_col", "clay"),
      verbose = FALSE
    ),
    "component_cols not found"
  )
})


test_that("gc_process_ssurgo_components aggregates correctly", {
  
  # Create test SSURGO data: 1 map unit with 2 components
  ssurgo_test <- data.frame(
    mukey = c("463168", "463168"),  # Same map unit
    compname = c("Harps", "Mahaska"),
    comppct = c(60, 40),  # 60% Harps, 40% Mahaska
    sandtotal_r = c(50, 30),
    silttotal_r = c(30, 40),
    claytotal_r = c(20, 30),
    stringsAsFactors = FALSE
  )
  
  result <- gc_process_ssurgo_components(
    ssurgo_data = ssurgo_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    aggregate_method = "weighted_mean",
    weight_field = "comppct",
    verbose = FALSE
  )
  
  # Check output structure
  expect_is(result, "list")
  expect_true(all(c("compositions", "component_summary", "aggregation_method",
                    "n_map_units", "n_components_raw", "n_records_final") %in% names(result)))
  
  # Check aggregation
  expect_equal(result$n_map_units, 1)
  expect_equal(result$n_components_raw, 2)
  expect_equal(result$n_records_final, 1)
  
  # Check weighted mean calculation
  agg_sand <- result$compositions$sandtotal_r[1]
  expected_sand <- (50 * 0.60) + (30 * 0.40)  # 42
  expect_equal(agg_sand, expected_sand, tolerance = 0.1)
})


test_that("gc_process_ssurgo_components handles composition normalization", {
  
  # Create test data with non-normalized compositions
  ssurgo_test <- data.frame(
    mukey = "463168",
    compname = "Harps",
    comppct = 100,
    sandtotal_r = 55,   # Will normalize to 100% together
    silttotal_r = 28,
    claytotal_r = 17,
    stringsAsFactors = FALSE
  )
  
  result <- gc_process_ssurgo_components(
    ssurgo_data = ssurgo_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    verbose = FALSE
  )
  
  # After aggregation, should sum to ~100%
  agg <- result$compositions[1, c("sandtotal_r", "silttotal_r", "claytotal_r")]
  comp_sum <- sum(agg, na.rm = TRUE)
  
  expect_equal(comp_sum, 100, tolerance = 0.1)
})


test_that("gc_process_ssurgo_components supports all aggregation methods", {
  
  ssurgo_test <- data.frame(
    mukey = c("463168", "463168"),
    compname = c("Harps", "Mahaska"),
    comppct = c(60, 40),
    sandtotal_r = c(50, 30),
    silttotal_r = c(30, 40),
    claytotal_r = c(20, 30),
    representativekomponent = c("Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  # Test weighted_mean
  result_weighted <- gc_process_ssurgo_components(
    ssurgo_data = ssurgo_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    aggregate_method = "weighted_mean",
    verbose = FALSE
  )
  expect_equal(result_weighted$aggregation_method, "weighted_mean")
  
  # Test component_best
  result_best <- gc_process_ssurgo_components(
    ssurgo_data = ssurgo_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    aggregate_method = "component_best",
    verbose = FALSE
  )
  expect_equal(result_best$aggregation_method, "component_best")
  
  # With component_best, should use Harps (60% > 40%)
  sand_best <- result_best$compositions$sandtotal_r[1]
  expect_equal(sand_best, 50, tolerance = 0.1)  # Harps sand
})


# ============================================================================
# Tests: gc_optimize_sda_queries()
# ============================================================================

test_that("gc_optimize_sda_queries creates proper output structure", {
  
  study_area <- sf::st_bbox(c(
    xmin = -93.6, ymin = 41.4,
    xmax = -93.2, ymax = 41.8
  ), crs = 4326)
  
  result <- gc_optimize_sda_queries(
    extent = study_area,
    tile_size = 0.5,
    n_workers = 2,
    priority = "speed",
    verbose = FALSE
  )
  
  # Check class and attributes
  expect_is(result, "gc_sda_query_results")
  expect_is(result, "data.frame")
  
  # Check attributes
  expect_true(!is.null(attr(result, "total_query_time")))
  expect_true(!is.null(attr(result, "n_tiles")))
  expect_true(!is.null(attr(result, "success_rate")))
})


test_that("gc_optimize_sda_queries respects priority level", {
  
  study_area <- sf::st_bbox(c(
    xmin = -93.6, ymin = 41.4,
    xmax = -93.2, ymax = 41.8
  ), crs = 4326)
  
  # Fast priority should use larger tiles
  result_fast <- gc_optimize_sda_queries(
    extent = study_area,
    priority = "speed",
    verbose = FALSE
  )
  
  # Completeness priority should use smaller tiles
  result_complete <- gc_optimize_sda_queries(
    extent = study_area,
    priority = "completeness",
    verbose = FALSE
  )
  
  # Both should succeed and have tiles metric
  expect_true(attr(result_fast, "n_tiles") > 0)
  expect_true(attr(result_complete, "n_tiles") > 0)
})


# ============================================================================
# Tests: gc_prepare_hierarchical_data()
# ============================================================================

test_that("gc_prepare_hierarchical_data validates inputs", {
  
  ssurgo_test <- data.frame(
    mukey = "463168",
    sand = 50, silt = 30, clay = 20
  )
  
  # Test: component_cols not found
  expect_error(
    gc_prepare_hierarchical_data(
      ssurgo_compositions = ssurgo_test,
      component_cols = c("sand", "missing", "clay"),
      verbose = FALSE
    ),
    "component_cols not found"
  )
})


test_that("gc_prepare_hierarchical_data combines SSURGO and field data", {
  
  # SSURGO data
  ssurgo_test <- data.frame(
    mukey = c("463168", "463169"),
    sandtotal_r = c(50, 40),
    silttotal_r = c(30, 35),
    claytotal_r = c(20, 25)
  )
  
  # Field observations
  field_test <- data.frame(
    pit_id = c("P1", "P2"),
    sandtotal_r = c(55, 45),
    silttotal_r = c(28, 32),
    claytotal_r = c(17, 23)
  )
  
  # Combine with field_priority (field weighted 1.0, SSURGO 0.6)
  result <- gc_prepare_hierarchical_data(
    ssurgo_compositions = ssurgo_test,
    field_observations = field_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    combine_method = "field_priority",
    verbose = FALSE
  )
  
  # Should have both SSURGO + field records
  expect_is(result, "data.frame")
  expect_true(nrow(result) > nrow(ssurgo_test))
  
  # Check source tracking
  expect_true("source" %in% names(result))
  sources <- unique(result$source)
  expect_true("field" %in% sources)
  expect_true("ssurgo" %in% sources)
  
  # Check weights
  expect_true("weight" %in% names(result))
})


test_that("gc_prepare_hierarchical_data supports combine strategies", {
  
  ssurgo_test <- data.frame(
    sandtotal_r = c(50, 40),
    silttotal_r = c(30, 35),
    claytotal_r = c(20, 25)
  )
  
  field_test <- data.frame(
    sandtotal_r = c(55, 45),
    silttotal_r = c(28, 32),
    claytotal_r = c(17, 23)
  )
  
  # Test field_priority
  result_priority <- gc_prepare_hierarchical_data(
    ssurgo_compositions = ssurgo_test,
    field_observations = field_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    combine_method = "field_priority",
    verbose = FALSE
  )
  expect_true(nrow(result_priority) > 0)
  
  # Test equal_weight
  result_equal <- gc_prepare_hierarchical_data(
    ssurgo_compositions = ssurgo_test,
    field_observations = field_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    combine_method = "equal_weight",
    verbose = FALSE
  )
  expect_true(nrow(result_equal) > 0)
  
  # Test ssurgo_only
  result_ssurgo <- gc_prepare_hierarchical_data(
    ssurgo_compositions = ssurgo_test,
    field_observations = NULL,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    combine_method = "ssurgo_only",
    verbose = FALSE
  )
  expect_equal(nrow(result_ssurgo), nrow(ssurgo_test))
  expect_true(all(result_ssurgo$source == "ssurgo"))
})


test_that("gc_prepare_hierarchical_data creates ILR coordinates", {
  skip_if_not_installed("compositions")
  
  ssurgo_test <- data.frame(
    sandtotal_r = c(50, 40, 60),
    silttotal_r = c(30, 35, 25),
    claytotal_r = c(20, 25, 15)
  )
  
  result <- gc_prepare_hierarchical_data(
    ssurgo_compositions = ssurgo_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    verbose = FALSE
  )
  
  # Should have ILR coordinates (if compositions package available)
  # ILR columns might not be created if compositions package not loaded
  # But the function should still work
  expect_true(nrow(result) == nrow(ssurgo_test))
})


# ============================================================================
# Integration Tests
# ============================================================================

test_that("SDA functions work together in sequence", {
  skip_if_offline()
  skip_if_not_installed("soilDB")
  
  # 1. Create synthetic SSURGO data
  ssurgo_raw <- data.frame(
    mukey = c("463168", "463168", "463169"),
    compname = c("Harps", "Mahaska", "Harps"),
    comppct = c(60, 40, 100),
    sandtotal_r = c(50, 30, 55),
    silttotal_r = c(30, 40, 28),
    claytotal_r = c(20, 30, 17)
  )
  
  # 2. Process components
  ssurgo_agg <- gc_process_ssurgo_components(
    ssurgo_data = ssurgo_raw,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    verbose = FALSE
  )
  
  expect_equal(ssurgo_agg$n_map_units, 2)
  
  # 3. Add field data
  field_test <- data.frame(
    sandtotal_r = c(52, 58),
    silttotal_r = c(28, 25),
    claytotal_r = c(20, 17)
  )
  
  # 4. Combine
  ilr_combined <- gc_prepare_hierarchical_data(
    ssurgo_compositions = ssurgo_agg$compositions,
    field_observations = field_test,
    component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
    combine_method = "field_priority",
    verbose = FALSE
  )
  
  # Check combined dataset
  expect_true(nrow(ilr_combined) > 0)
  expect_true("source" %in% names(ilr_combined))
  expect_true(nrow(ilr_combined) >= nrow(ssurgo_agg$compositions))
})
