test_that("gc_standardize_depths input validation", {
  profiles <- data.frame(
    mukey = c(1, 1, 2, 2),
    hzdept_r = c(0, 20, 0, 30),
    hzdepb_r = c(20, 50, 30, 100),
    clay = c(15, 25, 20, 35)
  )

  # Missing required columns
  expect_error(
    gc_standardize_depths(profiles[, -3]),
    "hzdepb_r"
  )

  # Invalid method
  expect_error(
    gc_standardize_depths(profiles, method = "invalid"),
    "method"
  )

  # Valid call should work
  result <- gc_standardize_depths(
    profiles,
    method = "linear",
    comp_cols = "clay",
    id_col = "mukey"
  )
  expect_s3_class(result, "data.frame")
})

test_that("gc_standardize_depths with linear method", {
  profiles <- data.frame(
    profile_id = c(1, 1, 1),
    hzdept_r = c(0, 20, 50),
    hzdepb_r = c(20, 50, 100),
    clay = c(10, 20, 30)
  )

  result <- gc_standardize_depths(
    profiles,
    method = "linear",
    intervals = c(20, 50, 100),
    comp_cols = "clay",
    id_col = "profile_id"
  )

  # Check structure
  expect_equal(nrow(result), 3)  # 3 intervals
  expect_true(all(c("depth_top", "depth_bot", "clay", "profile_id") %in% names(result)))

  # Check depth intervals
  expect_equal(result$depth_top, c(0, 20, 50))
  expect_equal(result$depth_bot, c(20, 50, 100))

  # Check profile IDs
  expect_equal(result$profile_id, c(1, 1, 1))
})

test_that("gc_standardize_depths with multiple components", {
  profiles <- data.frame(
    mukey = c(1, 1),
    hzdept_r = c(0, 30),
    hzdepb_r = c(30, 100),
    sand = c(60, 40),
    silt = c(20, 30),
    clay = c(20, 30)
  )

  result <- gc_standardize_depths(
    profiles,
    method = "linear",
    intervals = c(30, 100),
    comp_cols = c("sand", "silt", "clay"),
    id_col = "mukey"
  )

  # Check all components are present
  expect_true(all(c("sand", "silt", "clay") %in% names(result)))
  expect_equal(nrow(result), 2)  # 2 intervals
})

test_that("gc_standardize_depths with multiple profiles", {
  profiles <- data.frame(
    mukey = c(1, 1, 2, 2),
    hzdept_r = c(0, 30, 0, 50),
    hzdepb_r = c(30, 100, 50, 100),
    clay = c(15, 25, 20, 35)
  )

  result <- gc_standardize_depths(
    profiles,
    method = "linear",
    intervals = c(30, 100),
    comp_cols = "clay",
    id_col = "mukey"
  )

  # Should have 2 profiles × 2 intervals = 4 rows
  expect_equal(nrow(result), 4)

  # Check profile grouping
  prof1_rows <- result[result$mukey == 1, ]
  prof2_rows <- result[result$mukey == 2, ]
  expect_equal(nrow(prof1_rows), 2)
  expect_equal(nrow(prof2_rows), 2)
})

test_that("gc_standardize_depths GlobalSoilMap intervals", {
  profiles <- data.frame(
    pedon_id = 1,
    hzdept_r = c(0, 5, 15, 30, 60, 100),
    hzdepb_r = c(5, 15, 30, 60, 100, 150),
    clay = c(10, 12, 18, 25, 30, 32)
  )

  result <- gc_standardize_depths(
    profiles,
    method = "linear",
    intervals = c(5, 15, 30, 60, 100),  # GlobalSoilMap standard
    comp_cols = "clay",
    id_col = "pedon_id"
  )

  expect_equal(nrow(result), 5)
  expect_equal(result$depth_bot, c(5, 15, 30, 60, 100))
})

test_that("gc_standardize_depths custom intervals", {
  profiles <- data.frame(
    id = 1,
    hzdept_r = c(0, 10, 20),
    hzdepb_r = c(10, 20, 50),
    value = c(1, 2, 3)
  )

  # Custom 10 cm intervals
  result <- gc_standardize_depths(
    profiles,
    method = "linear",
    intervals = c(10, 20, 30, 40, 50),
    comp_cols = "value",
    id_col = "id"
  )

  expect_equal(nrow(result), 5)
  expect_equal(result$depth_bot, c(10, 20, 30, 40, 50))
})

test_that("gc_standardize_depths returns attributes", {
  profiles <- data.frame(
    id = 1,
    hzdept_r = 0,
    hzdepb_r = 50,
    clay = 20
  )

  result <- gc_standardize_depths(
    profiles,
    method = "linear",
    intervals = c(50),
    comp_cols = "clay",
    id_col = "id"
  )

  expect_equal(attr(result, "method"), "linear")
  expect_equal(attr(result, "intervals"), c(50))
  expect_equal(attr(result, "comp_cols"), "clay")
})

test_that("gc_validate_depth_standardization input validation", {
  profiles_original <- data.frame(
    id = c(1, 1),
    hzdept_r = c(0, 30),
    hzdepb_r = c(30, 100),
    clay = c(20, 30)
  )

  profiles_standardized <- data.frame(
    id = c(1, 1),
    depth_top = c(0, 30),
    depth_bot = c(30, 100),
    clay = c(20, 30)
  )
  attr(profiles_standardized, "comp_cols") <- "clay"

  # Valid call
  result <- gc_validate_depth_standardization(
    profiles_original, profiles_standardized,
    comp_cols = "clay", id_col = "id"
  )
  expect_s3_class(result, "data.frame")
})

test_that("gc_validate_depth_standardization detects mass conservation", {
  # Profile with exact mass conservation (linear interpolation of linear function)
  profiles_original <- data.frame(
    id = 1,
    hzdept_r = c(0, 50),
    hzdepb_r = c(50, 100),
    value = c(10, 20)
  )

  # Standardized perfectly
  profiles_standardized <- data.frame(
    id = 1,
    depth_top = c(0, 50),
    depth_bot = c(50, 100),
    value = c(10, 20)
  )
  attr(profiles_standardized, "comp_cols") <- "value"

  result <- gc_validate_depth_standardization(
    profiles_original, profiles_standardized,
    comp_cols = "value", id_col = "id"
  )

  expect_equal(result$mean_error_pct, 0, tolerance = 0.01)
  expect_equal(result$n_conserved, 1)
})

test_that("gc_validate_depth_standardization handles multiple profiles", {
  profiles_original <- data.frame(
    id = c(1, 1, 2, 2),
    hzdept_r = c(0, 50, 0, 50),
    hzdepb_r = c(50, 100, 50, 100),
    clay = c(15, 25, 20, 30)
  )

  profiles_standardized <- data.frame(
    id = c(1, 1, 2, 2),
    depth_top = c(0, 50, 0, 50),
    depth_bot = c(50, 100, 50, 100),
    clay = c(15, 25, 20, 30)
  )
  attr(profiles_standardized, "comp_cols") <- "clay"

  result <- gc_validate_depth_standardization(
    profiles_original, profiles_standardized,
    comp_cols = "clay", id_col = "id"
  )

  expect_equal(result$n_profiles, 2)
  expect_equal(result$n_conserved, 2)
})

test_that("gc_validate_depth_standardization with tolerance", {
  profiles_original <- data.frame(
    id = 1,
    hzdept_r = c(0, 50),
    hzdepb_r = c(50, 100),
    value = c(100, 200)
  )

  # Slightly different values (small error)
  profiles_standardized <- data.frame(
    id = 1,
    depth_top = c(0, 50),
    depth_bot = c(50, 100),
    value = c(100.5, 200.5)
  )
  attr(profiles_standardized, "comp_cols") <- "value"

  result_tight <- gc_validate_depth_standardization(
    profiles_original, profiles_standardized,
    comp_cols = "value", id_col = "id", tolerance = 0.001
  )
  expect_equal(result_tight$n_conserved, 0)  # Fails strict tolerance

  result_loose <- gc_validate_depth_standardization(
    profiles_original, profiles_standardized,
    comp_cols = "value", id_col = "id", tolerance = 0.01
  )
  expect_equal(result_loose$n_conserved, 1)  # Passes loose tolerance
})

test_that("gc_standardize_depths with spline method requires aqp", {
  skip_if_not_installed("aqp")

  profiles <- data.frame(
    id = c(1, 1, 1),
    hzdept_r = c(0, 20, 50),
    hzdepb_r = c(20, 50, 100),
    clay = c(10, 20, 30)
  )

  # Should work with aqp installed
  result <- gc_standardize_depths(
    profiles,
    method = "spline",
    intervals = c(50, 100),
    comp_cols = "clay",
    id_col = "id"
  )

  expect_s3_class(result, "data.frame")
  expect_equal(attr(result, "method"), "spline")
})

test_that("gc_standardize_depths spline method error without aqp", {
  if (requireNamespace("aqp", quietly = TRUE)) {
    skip("aqp is installed, skipping fallback test")
  }

  profiles <- data.frame(
    id = 1,
    hzdept_r = 0,
    hzdepb_r = 50,
    clay = 20
  )

  expect_error(
    gc_standardize_depths(profiles, method = "spline"),
    "aqp"
  )
})
