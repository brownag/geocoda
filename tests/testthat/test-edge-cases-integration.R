# Tests for edge cases, boundary conditions, and integration workflows

test_that("gc_handle_zeros works with all-zero rows", {
  # Create data with a completely zero row
  comp_data <- data.frame(
    sand = c(40, 0),
    silt = c(35, 0),
    clay = c(25, 0)
  )

  # Should handle without crashing
  result <- gc_handle_zeros(comp_data, method = "mzero")
  
  expect_is(result, "list")
  expect_true("imputed_data" %in% names(result))
  expect_true("n_zeros_imputed" %in% names(result))
  
  # Imputed data should all be positive
  expect_true(all(result$imputed_data > 0))
  
  # Rows should sum to 100 (or close)
  row_sums <- rowSums(result$imputed_data)
  expect_true(all(abs(row_sums - 100) < 1))
})

test_that("gc_expand_bounds handles tight constraints", {
  # Constraints where min+max = target_sum exactly
  constraints <- list(
    SAND = list(min = 30, max = 30),
    SILT = list(min = 50, max = 50),
    CLAY = list(min = 20, max = 20)
  )

  grid <- gc_expand_bounds(constraints, step = 1.0, target_sum = 100)
  
  expect_is(grid, "data.frame")
  expect_equal(nrow(grid), 1) # Only one valid composition
  expect_true(abs(sum(grid[1, ]) - 100) < 1e-6)
})

test_that("gc_ilr_params handles data with single observation", {
  # Single row of compositional data
  samples <- data.frame(
    sand = 30,
    silt = 50,
    clay = 20
  )

  # This should work (mean is the value, covariance might be NA/0)
  params <- gc_ilr_params(samples)
  
  expect_is(params, "list")
  expect_true("mean" %in% names(params))
  expect_true("cov" %in% names(params))
  expect_equal(length(params$mean), 2) # 3 components - 1
})

test_that("gc_resample_compositions with n = 0 returns empty data", {
  constraints <- list(
    SAND = list(min = 20, max = 30),
    SILT = list(min = 50, max = 60),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  
  result <- gc_resample_compositions(grid, n = 0, method = "uniform")
  
  expect_is(result$samples, "data.frame")
  expect_equal(nrow(result$samples), 0)
  expect_equal(result$n_sampled, 0)
})

test_that("gc_assess_stationarity handles uniform spatial patterns", {
  # Create uniform ILR data (no spatial trend)
  data_df <- data.frame(
    x = rep(seq(0, 100, by = 20), 2),
    y = rep(c(0, 50), each = 6),
    ilr1 = rnorm(12, mean = 0, sd = 0.1),
    ilr2 = rnorm(12, mean = 0, sd = 0.1)
  )

  # Should detect stationarity (low spatial variance)
  result <- gc_assess_stationarity(data_df, method = "pca")
  
  expect_is(result, "list")
  expect_true("is_stationary" %in% names(result))
  expect_true("spatial_variance" %in% names(result))
  expect_true(result$spatial_variance < 0.5)
})

test_that("Integration: Full workflow with small dataset", {
  # Minimal workflow test
  constraints <- list(
    SAND = list(min = 25, max = 35),
    SILT = list(min = 50, max = 60),
    CLAY = list(min = 10, max = 20)
  )

  # Step 1: Expand bounds
  grid <- gc_expand_bounds(constraints, step = 2.5, target_sum = 100)
  expect_true(nrow(grid) > 0)

  # Step 2: Resample
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 15, seed = 55)
  samples <- samples_res$samples
  expect_equal(nrow(samples), 15)

  # Step 3: Estimate parameters
  ilr_params <- gc_ilr_params(samples)
  expect_equal(length(ilr_params$mean), 2)

  # Step 4: Get variogram suggestions
  extent <- c(0, 0, 100, 100)
  vgm_suggestions <- gc_vgm_defaults(ilr_params, extent)
  expect_true(vgm_suggestions$range > 0)

  # Step 5: Build model
  vgm_model <- gstat::vgm(
    psill = vgm_suggestions$mean_sill,
    model = "Exp",
    range = vgm_suggestions$range,
    nugget = vgm_suggestions$nugget
  )
  model <- gc_ilr_model(ilr_params, vgm_model)
  expect_s3_class(model, "gstat")

  # Step 6: Create simulation grid
  sim_grid_df <- data.frame(
    x = c(25, 75),
    y = c(25, 75)
  )
  sim_grid <- sf::st_as_sf(sim_grid_df, coords = c("x", "y"), crs = "local")

  # Step 7: Simulate
  sims <- gc_sim_composition(
    model, sim_grid,
    nsim = 2,
    target_names = c("sand", "silt", "clay")
  )
  
  expect_s4_class(sims, "SpatRaster")
  
  # Step 8: Validate
  sims_matrix <- terra::as.matrix(sims, wide = TRUE)
  row_sums <- rowSums(sims_matrix, na.rm = TRUE)
  expect_true(all(abs(row_sums[1:2] - 100) < 1e-6)) # First 2 rows (first realization at 2 locations)
  expect_true(all(abs(row_sums[3:4] - 100) < 1e-6)) # Next 2 rows (second realization)
})

test_that("Integration: Assessment functions work on real ensemble", {
  # Create ensemble data
  constraints <- list(
    SAND = list(min = 20, max = 35),
    SILT = list(min = 45, max = 65),
    CLAY = list(min = 10, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 50, seed = 56)
  samples <- samples_res$samples

  # Test diagnostics on ensemble
  stationarity_results <- gc_assess_stationarity(
    cbind(data.frame(x = 1:50, y = 1:50), samples),
    method = "pca"
  )
  expect_is(stationarity_results, "list")

  gaussianity_results <- gc_assess_gaussianity(
    samples,
    methods = c("ad", "ks", "sw")
  )
  expect_is(gaussianity_results, "data.frame")
  expect_equal(nrow(gaussianity_results), 3)
})

test_that("gc_identify_strata works with small clusters", {
  # Create data with clear clusters
  cluster1 <- data.frame(
    x = c(10, 12, 11),
    y = c(10, 12, 11),
    sand = c(25, 26, 24),
    silt = c(55, 54, 56),
    clay = c(20, 20, 20)
  )
  
  cluster2 <- data.frame(
    x = c(85, 87, 86),
    y = c(85, 87, 86),
    sand = c(35, 34, 36),
    silt = c(45, 46, 44),
    clay = c(20, 20, 20)
  )
  
  data_combined <- rbind(cluster1, cluster2)
  
  # Should identify 2 clusters
  result <- gc_identify_strata(data_combined, method = "hierarchical", k = 2)
  
  expect_is(result, "list")
  expect_true("clusters" %in% names(result))
  expect_equal(length(unique(result$clusters)), 2)
})

test_that("gc_fit_vgm with correct.diagonal parameter affects admissibility", {
  constraints <- list(
    SAND = list(min = 20, max = 30),
    SILT = list(min = 50, max = 60),
    CLAY = list(min = 15, max = 25)
  )

  grid <- gc_expand_bounds(constraints, step = 5, target_sum = 100)
  samples_res <- gc_resample_compositions(grid, method = "uniform", n = 40, seed = 57)
  samples <- samples_res$samples
  ilr_params <- gc_ilr_params(samples)

  # Create sample locations
  set.seed(42)
  sample_locs <- expand.grid(
    x = seq(10, 90, by = 20),
    y = seq(10, 90, by = 20)
  )
  
  n_samples <- nrow(sample_locs)
  samples_subset <- samples[sample(nrow(samples), n_samples, replace = TRUE), ]
  ilr_samples <- compositions::ilr(compositions::acomp(samples_subset))
  colnames(ilr_samples) <- paste0("ilr", seq_len(ncol(ilr_samples)))
  
  sample_data_df <- cbind(sample_locs, as.data.frame(ilr_samples))

  # Test with different correction factors
  result_low <- gc_fit_vgm(ilr_params, sample_data_df, correct.diagonal = 1.00)
  result_high <- gc_fit_vgm(ilr_params, sample_data_df, correct.diagonal = 1.05)

  expect_is(result_low, "data.frame")
  expect_is(result_high, "data.frame")
  
  # Both should have psill columns reflecting the corrected diagonal
  expect_true(all(result_high$psill >= result_low$psill))
})
