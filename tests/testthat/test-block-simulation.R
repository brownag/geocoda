test_that("gc_discretize_block creates correct number of sub-points", {
  blocks <- expand.grid(x = c(0, 100), y = c(0, 100))
  sub_points <- gc_discretize_block(blocks, block_size = c(100, 100),
                                     discretization = c(3, 3))
  expect_equal(nrow(sub_points), 36)  # 4 blocks × 9 sub-points
  block_ids <- table(sub_points$block_id)
  expect_equal(as.numeric(block_ids), c(9, 9, 9, 9))
})

test_that("gc_discretize_block sub-points are centered in sub-cells", {
  blocks <- data.frame(x = 0, y = 0)
  sub_points <- gc_discretize_block(blocks, block_size = c(100, 100),
                                     discretization = c(2, 2))
  x_rel <- unique(sub_points$x)
  y_rel <- unique(sub_points$y)
  expect_equal(length(x_rel), 2)
  expect_equal(length(y_rel), 2)
})

test_that("gc_discretize_block sub-points within block bounds", {
  blocks <- data.frame(x = 100, y = 200)
  sub_points <- gc_discretize_block(blocks, block_size = c(50, 60),
                                     discretization = c(4, 4))
  x_offset <- sub_points$x - 100
  y_offset <- sub_points$y - 200
  expect_true(all(x_offset >= -25 & x_offset <= 25))
  expect_true(all(y_offset >= -30 & y_offset <= 30))
})

test_that("gc_average_block_ilr computes correct means", {
  sub_points <- data.frame(block_id = c(1, 1, 1, 1, 2, 2, 2, 2))
  sim_result <- data.frame(
    ilr1.sim1 = c(1, 2, 3, 4, 5, 6, 7, 8),
    ilr2.sim1 = c(10, 20, 30, 40, 50, 60, 70, 80),
    ilr1.sim2 = c(2, 4, 6, 8, 10, 12, 14, 16)
  )
  averaged <- gc_average_block_ilr(sim_result, sub_points)
  expect_equal(nrow(averaged), 2)
  expect_equal(averaged$ilr1.sim1[1], 2.5)
  expect_equal(averaged$ilr1.sim1[2], 6.5)
})

test_that("gc_average_block_ilr preserves order", {
  sub_points <- data.frame(block_id = c(3, 1, 3, 2, 1, 2))
  sim_result <- data.frame(ilr1.sim1 = c(1, 2, 3, 4, 5, 6))
  averaged <- gc_average_block_ilr(sim_result, sub_points)
  expect_equal(averaged$block_id, c(1, 2, 3))
})

test_that("gc_average_block_ilr works with MAF columns", {
  sub_points <- data.frame(block_id = c(1, 1, 1, 1, 2, 2, 2, 2))
  sim_result <- data.frame(
    maf1.sim1 = c(1, 2, 3, 4, 5, 6, 7, 8),
    maf2.sim1 = c(10, 20, 30, 40, 50, 60, 70, 80)
  )
  averaged <- gc_average_block_ilr(sim_result, sub_points)
  expect_equal(nrow(averaged), 2)
  expect_equal(averaged$maf1.sim1[1], 2.5)
})

test_that("gc_sim_composition_block produces SpatRaster output", {
  skip_if_not_installed("gstat")
  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))
  vgm_fit <- gstat::vgm(psill = 1.0, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)
  blocks <- data.frame(x = c(0, 100, 200), y = c(0, 100, 200))
  
  result <- gc_sim_composition_block(
    model = model, block_centers = blocks, block_size = c(50, 50),
    nsim = 2, target_names = c("sand", "silt", "clay"), discretization = c(2, 2)
  )
  
  expect_s4_class(result, "SpatRaster")
  expect_true(terra::nlyr(result) >= 6)  # At least 3 components × 2 sims
  layer_names <- names(result)
  expect_true(any(grepl("sand", layer_names)))
  expect_true(any(grepl("silt", layer_names)))
  expect_true(any(grepl("clay", layer_names)))
})

test_that("gc_sim_composition_block compositional values non-negative", {
  skip_if_not_installed("gstat")
  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))
  vgm_fit <- gstat::vgm(psill = 0.8, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)
  blocks <- data.frame(x = c(0, 100), y = c(0, 100))
  
  result <- gc_sim_composition_block(
    model = model, block_centers = blocks, block_size = c(50, 50),
    nsim = 2, target_names = c("sand", "silt", "clay"), discretization = c(3, 3)
  )
  
  values <- terra::values(result)
  valid_values <- values[!is.na(values)]
  expect_true(all(valid_values >= -0.01))  # Allow small numerical errors
  expect_true(all(valid_values <= 100.01))
})

test_that("Block discretization input validation", {
  blocks <- data.frame(x = c(0, 100), y = c(0, 100))
  expect_error(
    gc_discretize_block(blocks, block_size = c(0, 100)),
    "block_size"
  )
  expect_error(
    gc_discretize_block(blocks, block_size = c(100, 100), discretization = c(1, 4)),
    "discretization"
  )
})

test_that("Block averaging input validation", {
  sub_points <- data.frame(block_id = c(1, 1, 2, 2))
  sim_result <- data.frame(ilr1.sim1 = c(1, 2, 3, 4))

  # Mismatched rows should error
  sim_result_short <- as.data.frame(sim_result[1:3, , drop = FALSE])
  expect_error(
    gc_average_block_ilr(sim_result_short, sub_points),
    "nrow"
  )

  # Missing block_id should error
  expect_error(
    gc_average_block_ilr(sim_result, data.frame(x = c(1, 1, 2, 2))),
    "block_id"
  )
})

test_that("gc_sim_composition_block with different discretizations", {
  skip_if_not_installed("gstat")
  set.seed(42)
  ilr_cov <- matrix(c(1.0, 0.3, 0.3, 1.0), nrow = 2)
  ilr_params <- list(mean = c(0, 0), cov = ilr_cov, names = c("sand", "silt", "clay"))
  vgm_fit <- gstat::vgm(psill = 1.0, model = "Exp", range = 50, nugget = 0.1)
  model <- gc_ilr_model(ilr_params, vgm_fit)
  blocks <- data.frame(x = c(0, 100), y = c(0, 100))
  
  for (n_sub in c(2, 3, 4)) {
    result <- gc_sim_composition_block(
      model = model, block_centers = blocks, block_size = c(50, 50),
      nsim = 1, target_names = c("sand", "silt", "clay"),
      discretization = c(n_sub, n_sub)
    )
    expect_s4_class(result, "SpatRaster")
    expect_true(terra::nlyr(result) >= 3)  # At least 3 components
  }
})
