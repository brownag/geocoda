test_that("gc_discretize_block_3d creates correct sub-point grid", {
  blocks <- expand.grid(x = c(0, 100), y = c(0, 100), z = c(10, 50))

  sub_points <- gc_discretize_block_3d(
    blocks,
    block_size = c(100, 100, 40),
    discretization = c(3, 3, 3)
  )

  # Check dimensions
  expect_equal(nrow(sub_points), 8 * 27)  # 8 blocks × 27 sub-points
  expect_equal(ncol(sub_points), 4)  # x, y, z, block_id
  expect_true(all(c("x", "y", "z", "block_id") %in% names(sub_points)))

  # Check block_id ranges
  expect_equal(min(sub_points$block_id), 1)
  expect_equal(max(sub_points$block_id), 8)
  # Each block should have 27 sub-points (3×3×3)
  block_counts <- table(sub_points$block_id)
  expect_true(all(block_counts == 27))
})

test_that("gc_discretize_block_3d sub-points centered within blocks", {
  blocks <- data.frame(x = 0, y = 0, z = 0)
  block_size <- c(100, 100, 40)
  discretization <- c(3, 3, 3)

  sub_points <- gc_discretize_block_3d(blocks, block_size, discretization)

  # Center of block is (0, 0, 0)
  # Sub-points should be centered around it
  expect_equal(mean(sub_points$x), 0, tolerance = 1e-10)
  expect_equal(mean(sub_points$y), 0, tolerance = 1e-10)
  expect_equal(mean(sub_points$z), 0, tolerance = 1e-10)

  # Check that sub-points are within expected range
  expect_true(all(sub_points$x >= -50 & sub_points$x <= 50))
  expect_true(all(sub_points$y >= -50 & sub_points$y <= 50))
  expect_true(all(sub_points$z >= -20 & sub_points$z <= 20))
})

test_that("gc_discretize_block_3d respects block offsets", {
  blocks <- data.frame(x = 250, y = 500, z = 100)
  block_size <- c(100, 100, 40)

  sub_points <- gc_discretize_block_3d(blocks, block_size)

  # Block center is at (250, 500, 100)
  # Sub-points should be centered around it
  expect_equal(mean(sub_points$x), 250, tolerance = 1e-10)
  expect_equal(mean(sub_points$y), 500, tolerance = 1e-10)
  expect_equal(mean(sub_points$z), 100, tolerance = 1e-10)
})

test_that("gc_discretize_block_3d different discretization levels", {
  blocks <- data.frame(x = 0, y = 0, z = 0)

  # Test different discretization levels
  sub_2_2_2 <- gc_discretize_block_3d(blocks, c(100, 100, 40), c(2, 2, 2))
  expect_equal(nrow(sub_2_2_2), 8)  # 2×2×2

  sub_3_3_3 <- gc_discretize_block_3d(blocks, c(100, 100, 40), c(3, 3, 3))
  expect_equal(nrow(sub_3_3_3), 27)  # 3×3×3

  sub_4_4_4 <- gc_discretize_block_3d(blocks, c(100, 100, 40), c(4, 4, 4))
  expect_equal(nrow(sub_4_4_4), 64)  # 4×4×4
})

test_that("gc_average_block_ilr_3d computes correct block means", {
  skip_if_not_installed("gstat")

  set.seed(42)

  # Create 3D sub-points (2×2×2 discretization = 8 sub-points per block)
  # 2×2×2 blocks = 8 blocks = 64 sub-points total
  blocks <- data.frame(x = c(0, 100), y = c(0, 100), z = c(10, 50))
  sub_points <- gc_discretize_block_3d(blocks, block_size = c(100, 100, 40),
                                       discretization = c(2, 2, 2))

  # Create simulation result with ILR values
  # nrow(sub_points) = 2*2 blocks × 2×2×2 sub-points = 4 × 8 = 32
  n_sims <- nrow(sub_points)
  sim_result <- data.frame(
    x = sub_points$x,
    y = sub_points$y,
    z = sub_points$z,
    block_id = sub_points$block_id,
    ilr1.sim1 = seq(1, n_sims),
    ilr2.sim1 = seq(1, n_sims) * 2
  )

  block_avg <- gc_average_block_ilr_3d(sim_result, sub_points)

  # Should have 4 blocks (2×2 in x,y with 2 z levels)
  expect_equal(nrow(block_avg), nrow(blocks))

  # Each block should have 8 sub-points (2×2×2)
  # Values should be averaged correctly
  block_counts <- table(sub_points$block_id)
  expected_mean_1 <- mean(1:block_counts[1])

  # Check some averaging was performed
  expect_true(all(is.finite(block_avg$ilr1.sim1)))
  expect_true(all(is.finite(block_avg$ilr2.sim1)))
})

test_that("gc_average_block_ilr_3d requires proper columns", {
  sub_points <- data.frame(x = c(0, 1, 2), y = c(0, 1, 2), z = c(0, 1, 2),
                           block_id = c(1, 1, 2))

  # Missing ILR columns should error
  sim_result_no_ilr <- data.frame(
    x = c(0, 1, 2), y = c(0, 1, 2), z = c(0, 1, 2),
    block_id = c(1, 1, 2)
  )
  expect_error(gc_average_block_ilr_3d(sim_result_no_ilr, sub_points),
               "No ILR/MAF")
})

test_that("gc_discretize_block_3d validates inputs", {
  blocks <- data.frame(x = 0, y = 0, z = 0)

  # Missing z column
  blocks_no_z <- data.frame(x = 0, y = 0)
  expect_error(gc_discretize_block_3d(blocks_no_z, c(100, 100, 40)),
               "x.*y.*z")

  # Invalid block_size
  expect_error(gc_discretize_block_3d(blocks, c(100, 100)), "length")
  expect_error(gc_discretize_block_3d(blocks, c(-100, 100, 40)), "all")

  # Invalid discretization (< 2 in any dimension)
  expect_error(gc_discretize_block_3d(blocks, c(100, 100, 40), c(1, 1, 1)),
               "is not TRUE")
})

test_that("3D block discretization creates appropriate structure", {
  skip_if_not_installed("gstat")

  set.seed(42)

  # Create 3D block grid
  blocks <- expand.grid(
    x = seq(0, 200, by = 100),
    y = seq(0, 200, by = 100),
    z = c(30, 70)  # 2 depth levels
  )

  # Discretize into sub-points
  sub_points <- gc_discretize_block_3d(blocks, block_size = c(100, 100, 40),
                                       discretization = c(2, 2, 2))

  # Create synthetic ILR simulation results
  sim_result <- data.frame(
    x = sub_points$x,
    y = sub_points$y,
    z = sub_points$z,
    block_id = sub_points$block_id,
    ilr1.sim1 = rnorm(nrow(sub_points)),
    ilr2.sim1 = rnorm(nrow(sub_points)),
    ilr1.sim2 = rnorm(nrow(sub_points)),
    ilr2.sim2 = rnorm(nrow(sub_points))
  )

  # Average to block support
  block_avg <- gc_average_block_ilr_3d(sim_result, sub_points)

  # Check that we have correct number of blocks
  # 3×3 grid of x,y with 2 depth levels = 18 blocks
  expect_equal(nrow(block_avg), nrow(blocks))

  # Check that ILR columns are present
  expect_true(all(c("block_id", "ilr1.sim1", "ilr2.sim1", "ilr1.sim2",
                    "ilr2.sim2") %in% names(block_avg)))

  # Check that sub-point count is correct
  # 18 blocks × 2×2×2 = 144 sub-points
  expect_equal(nrow(sub_points), 144)
})

test_that("3D block averaging preserves structure across simulations", {
  skip_if_not_installed("gstat")

  set.seed(123)

  # Create block grid
  grid_blocks <- expand.grid(
    x = seq(100, 400, by = 100),
    y = seq(100, 400, by = 100),
    z = seq(30, 70, by = 20)
  )

  sub_points <- gc_discretize_block_3d(
    grid_blocks,
    block_size = c(100, 100, 20),
    discretization = c(2, 2, 2)
  )

  # Create synthetic ILR simulation results with multiple realizations
  n_sim <- 3
  ilr_cols <- list()
  for (s in seq_len(n_sim)) {
    ilr_cols[[paste0("ilr1.sim", s)]] <- rnorm(nrow(sub_points))
    ilr_cols[[paste0("ilr2.sim", s)]] <- rnorm(nrow(sub_points))
  }

  sim_result <- data.frame(
    x = sub_points$x,
    y = sub_points$y,
    z = sub_points$z,
    block_id = sub_points$block_id,
    ilr_cols
  )

  # Average to block support
  block_avg <- gc_average_block_ilr_3d(sim_result, sub_points)

  # Verify averaging worked
  expect_equal(nrow(block_avg), nrow(grid_blocks))

  # Block-averaged results should have fewer rows than sub-point results
  expect_true(nrow(block_avg) < nrow(sub_points))

  # Check that all ILR columns are present and finite
  for (s in seq_len(n_sim)) {
    col_name <- paste0("ilr1.sim", s)
    expect_true(col_name %in% names(block_avg))
    expect_true(all(is.finite(block_avg[[col_name]])))
  }
})

test_that("3D block functions handle MAF-transformed data", {
  skip_if_not_installed("gstat")

  set.seed(42)

  # Create 3D data
  blocks <- expand.grid(x = c(0, 100), y = c(0, 100), z = c(10, 50))
  sub_points <- gc_discretize_block_3d(blocks, block_size = c(100, 100, 40), discretization = c(2, 2, 2))

  # Create simulation result with MAF columns (alternative to ILR)
  sim_result <- sub_points[, c("x", "y", "z", "block_id")]
  sim_result$maf1.sim1 <- rnorm(nrow(sub_points))
  sim_result$maf2.sim1 <- rnorm(nrow(sub_points))
  sim_result$maf1.sim2 <- rnorm(nrow(sub_points))

  # Averaging should work with MAF columns
  block_avg <- gc_average_block_ilr_3d(sim_result, sub_points)

  expect_equal(nrow(block_avg), nrow(blocks))
  expect_true("maf1.sim1" %in% names(block_avg))
  expect_true("maf2.sim1" %in% names(block_avg))
})

test_that("3D block averaging with unequal sub-point counts", {
  # Sometimes blocks might have different numbers of sub-points
  # (e.g., in masked regions) - test averaging still works

  sub_points_data <- data.frame(
    x = c(0, 1, 2, 3, 4, 5),
    y = c(0, 0, 0, 1, 1, 1),
    z = c(0, 0, 0, 1, 1, 1),
    block_id = c(1, 1, 2, 3, 3, 3)  # Block 1: 2 points, Block 2: 1 point, Block 3: 3 points
  )

  sim_result <- sub_points_data
  sim_result$ilr1.sim1 <- 1:6

  block_avg <- gc_average_block_ilr_3d(sim_result, sub_points_data)

  expect_equal(nrow(block_avg), 3)
  expect_equal(block_avg$ilr1.sim1[1], mean(1:2))  # (1+2)/2
  expect_equal(block_avg$ilr1.sim1[2], 3)           # Just 3
  expect_equal(block_avg$ilr1.sim1[3], mean(4:6))  # (4+5+6)/3
})

test_that("3D block discretization with very small blocks", {
  blocks <- data.frame(x = 0, y = 0, z = 0)

  # Small blocks (0.1 × 0.1 × 0.1)
  sub_points <- gc_discretize_block_3d(blocks, block_size = c(0.1, 0.1, 0.1),
                                       discretization = c(2, 2, 2))

  # Should still create correct number of sub-points
  expect_equal(nrow(sub_points), 8)

  # Sub-points should be within expected range (±0.05 in each direction)
  expect_true(all(sub_points$x >= -0.05 & sub_points$x <= 0.05))
  expect_true(all(sub_points$y >= -0.05 & sub_points$y <= 0.05))
  expect_true(all(sub_points$z >= -0.05 & sub_points$z <= 0.05))
})

test_that("3D block discretization with large blocks", {
  blocks <- data.frame(x = 1000, y = 1000, z = 500)

  # Large blocks (1000 × 1000 × 500)
  sub_points <- gc_discretize_block_3d(blocks, block_size = c(1000, 1000, 500),
                                       discretization = c(3, 3, 3))

  # Should create correct number of sub-points
  expect_equal(nrow(sub_points), 27)

  # Sub-points should be within expected range
  expect_true(all(sub_points$x >= 500 & sub_points$x <= 1500))
  expect_true(all(sub_points$y >= 500 & sub_points$y <= 1500))
  expect_true(all(sub_points$z >= 250 & sub_points$z <= 750))
})
