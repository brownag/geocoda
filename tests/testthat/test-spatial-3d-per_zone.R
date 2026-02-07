context("3D spatial stratification")

# Helper to create 3D test data
create_3d_test_data <- function(n_per_zone = 40) {
  set.seed(42)
  zones <- rep(c("z1", "z2", "z3"), each = n_per_zone)

  x <- rep(1:3, each = n_per_zone) * 30 + rnorm(3 * n_per_zone, sd = 5)
  y <- rnorm(3 * n_per_zone, mean = 50, sd = 10)
  z <- rnorm(3 * n_per_zone, mean = 40, sd = 20)

  sand <- c(
    rnorm(n_per_zone, mean = 70, sd = 8),
    rnorm(n_per_zone, mean = 50, sd = 10),
    rnorm(n_per_zone, mean = 30, sd = 12)
  )
  silt <- c(
    rnorm(n_per_zone, mean = 20, sd = 5),
    rnorm(n_per_zone, mean = 30, sd = 7),
    rnorm(n_per_zone, mean = 45, sd = 8)
  )
  clay <- c(
    rnorm(n_per_zone, mean = 10, sd = 3),
    rnorm(n_per_zone, mean = 20, sd = 5),
    rnorm(n_per_zone, mean = 25, sd = 7)
  )

  # Closure
  denom <- sand + silt + clay
  sand <- (sand / denom) * 100
  silt <- (silt / denom) * 100
  clay <- (clay / denom) * 100

  data.frame(
    zone = zones,
    x = x,
    y = y,
    z = z,
    SAND = sand,
    SILT = silt,
    CLAY = clay
  )
}

# ============================================================================
# gc_fit_vgm_3d_per_zone tests
# ============================================================================

test_that("gc_fit_vgm_3d_per_zone validates input data", {
  expect_error(
    gc_fit_vgm_3d_per_zone(
      data = list(),
      zone_col = "zone",
      value_col = "SAND"
    ),
    "data must be a data frame"
  )
})

test_that("gc_fit_vgm_3d_per_zone checks zone_col exists", {
  data <- create_3d_test_data()
  expect_error(
    gc_fit_vgm_3d_per_zone(
      data = data,
      zone_col = "nonexistent",
      value_col = "SAND"
    ),
    "zone_col.*not found"
  )
})

test_that("gc_fit_vgm_3d_per_zone checks value_col exists", {
  data <- create_3d_test_data()
  expect_error(
    gc_fit_vgm_3d_per_zone(
      data = data,
      zone_col = "zone",
      value_col = "nonexistent"
    ),
    "value_col.*not found"
  )
})

test_that("gc_fit_vgm_3d_per_zone requires xyz coordinates", {
  data <- create_3d_test_data()
  data_no_z <- data[, !(names(data) == "z")]

  expect_error(
    gc_fit_vgm_3d_per_zone(
      data = data_no_z,
      zone_col = "zone",
      value_col = "SAND"
    ),
    "must contain x, y, z"
  )
})

test_that("gc_fit_vgm_3d_per_zone returns list with required elements", {
  data <- create_3d_test_data()
  result <- gc_fit_vgm_3d_per_zone(
    data = data,
    zone_col = "zone",
    value_col = "SAND"
  )

  expect_is(result, "list")
  expect_true("zone_vgm" %in% names(result))
  expect_true("zone_ranges" %in% names(result))
  expect_true("zone_sill" %in% names(result))
  expect_true("data_summary" %in% names(result))
})

test_that("gc_fit_vgm_3d_per_zone creates vgm per zone", {
  data <- create_3d_test_data()
  result <- gc_fit_vgm_3d_per_zone(
    data = data,
    zone_col = "zone",
    value_col = "SAND"
  )

  n_zones <- length(unique(data$zone))
  expect_equal(length(result$zone_vgm), n_zones)
})

test_that("gc_fit_vgm_3d_per_zone computes anisotropy ratios", {
  data <- create_3d_test_data()
  result <- gc_fit_vgm_3d_per_zone(
    data = data,
    zone_col = "zone",
    value_col = "SAND"
  )

  expect_true("anisotropy_ratio" %in% names(result$zone_ranges))
  expect_true(all(result$zone_ranges$anisotropy_ratio > 0))
})

test_that("gc_fit_vgm_3d_per_zone lateral range > vertical range", {
  data <- create_3d_test_data()
  result <- gc_fit_vgm_3d_per_zone(
    data = data,
    zone_col = "zone",
    value_col = "SAND"
  )

  expect_true(all(result$zone_ranges$lateral_range > 
                  result$zone_ranges$vertical_range))
})

test_that("gc_fit_vgm_3d_per_zone respects range_limiter", {
  data <- create_3d_test_data()
  result_1 <- gc_fit_vgm_3d_per_zone(
    data = data,
    zone_col = "zone",
    value_col = "SAND",
    range_limiter = 1.0
  )
  result_05 <- gc_fit_vgm_3d_per_zone(
    data = data,
    zone_col = "zone",
    value_col = "SAND",
    range_limiter = 0.5
  )

  # With smaller range_limiter, ranges should be smaller
  expect_true(mean(result_05$zone_ranges$lateral_range) <
              mean(result_1$zone_ranges$lateral_range))
})

test_that("gc_fit_vgm_3d_per_zone data_summary has correct counts", {
  data <- create_3d_test_data(n_per_zone = 30)
  result <- gc_fit_vgm_3d_per_zone(
    data = data,
    zone_col = "zone",
    value_col = "SAND"
  )

  for (i in seq_len(nrow(result$data_summary))) {
    zone_name <- result$data_summary$zone[i]
    expected_n <- sum(data$zone == zone_name)
    actual_n <- result$data_summary$n_obs[i]
    expect_equal(actual_n, expected_n)
  }
})

# ============================================================================
# gc_define_nested_hierarchy tests
# ============================================================================

test_that("gc_define_nested_hierarchy validates domain_zones", {
  expect_error(
    gc_define_nested_hierarchy(
      domain_zones = NULL,
      depth_strata = c("shallow", "deep"),
      depth_breaks = c(0, 50, 100)
    ),
    "domain_zones must be non-empty character"
  )
})

test_that("gc_define_nested_hierarchy validates depth_strata", {
  expect_error(
    gc_define_nested_hierarchy(
      domain_zones = c("z1", "z2"),
      depth_strata = NULL,
      depth_breaks = c(0, 50, 100)
    ),
    "depth_strata must be non-empty character"
  )
})

test_that("gc_define_nested_hierarchy validates depth_breaks length", {
  expect_error(
    gc_define_nested_hierarchy(
      domain_zones = c("z1", "z2"),
      depth_strata = c("shallow", "deep"),
      depth_breaks = c(0, 100)  # wrong length
    ),
    "depth_breaks must have length"
  )
})

test_that("gc_define_nested_hierarchy validates breaks are increasing", {
  expect_error(
    gc_define_nested_hierarchy(
      domain_zones = c("z1", "z2"),
      depth_strata = c("shallow", "deep"),
      depth_breaks = c(0, 100, 50)  # not increasing
    ),
    "must be strictly increasing"
  )
})

test_that("gc_define_nested_hierarchy creates correct hierarchy size", {
  result <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2", "z3"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  expect_equal(result$n_units, 3 * 2)  # 3 domains × 2 depths
})

test_that("gc_define_nested_hierarchy returns data frame with nested_units", {
  result <- gc_define_nested_hierarchy(
    domain_zones = c("upland", "lowland"),
    depth_strata = c("surface", "deep"),
    depth_breaks = c(0, 30, 100)
  )

  expect_true("nested_units" %in% names(result))
  expect_is(result$nested_units, "data.frame")
  expect_equal(nrow(result$nested_units), 4)
})

test_that("gc_define_nested_hierarchy hierarchy_id is unique", {
  result <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2", "z3"),
    depth_strata = c("shallow", "mid", "deep"),
    depth_breaks = c(0, 30, 60, 100)
  )

  expect_equal(length(unique(result$nested_units$hierarchy_id)), 
               result$n_units)
})

# ============================================================================
# gc_sim_hierarchical_3d_per_zone tests
# ============================================================================

test_that("gc_sim_hierarchical_3d_per_zone validates object type", {
  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  expect_error(
    gc_sim_hierarchical_3d_per_zone(
      object = NULL,
      n = 10,
      nested_hierarchy = h_nested,
      coords = data.frame(x = 1:3, y = 1:3, z = c(10, 30, 70)),
      zone_vector = rep("z1", 3)
    ),
    "object must be a fitted"
  )
})

test_that("gc_sim_hierarchical_3d_per_zone validates n parameter", {
  data <- create_3d_test_data()
  hz_def <- gc_define_hierarchy(
    data = data,
    comp_cols = c("SAND", "SILT", "CLAY"),
    group_col = "zone"
  )
  hz_fit <- gc_fit_hierarchical_model(data = data, object = hz_def)

  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  coords <- data.frame(x = 1:2, y = 1:2, z = c(10, 30))

  expect_error(
    gc_sim_hierarchical_3d_per_zone(
      object = hz_fit,
      n = -1,
      nested_hierarchy = h_nested,
      coords = coords,
      zone_vector = rep("z1", 2)
    ),
    "n must be positive"
  )
})

test_that("gc_sim_hierarchical_3d_per_zone requires xyz coordinates", {
  data <- create_3d_test_data()
  hz_def <- gc_define_hierarchy(
    data = data,
    comp_cols = c("SAND", "SILT", "CLAY"),
    group_col = "zone"
  )
  hz_fit <- gc_fit_hierarchical_model(data = data, object = hz_def)

  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  coords <- data.frame(x = 1:2, y = 1:2)  # missing z

  expect_error(
    gc_sim_hierarchical_3d_per_zone(
      object = hz_fit,
      n = 10,
      nested_hierarchy = h_nested,
      coords = coords,
      zone_vector = rep("z1", 2)
    ),
    "must contain x, y, z"
  )
})

test_that("gc_sim_hierarchical_3d_per_zone validates zone vector length", {
  data <- create_3d_test_data()
  hz_def <- gc_define_hierarchy(
    data = data,
    comp_cols = c("SAND", "SILT", "CLAY"),
    group_col = "zone"
  )
  hz_fit <- gc_fit_hierarchical_model(data = data, object = hz_def)

  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2", "z3"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  coords <- data.frame(x = 1:5, y = 1:5, z = c(10, 30, 50, 70, 90))

  expect_error(
    gc_sim_hierarchical_3d_per_zone(
      object = hz_fit,
      n = 10,
      nested_hierarchy = h_nested,
      coords = coords,
      zone_vector = rep("z1", 3)  # wrong length
    ),
    "zone_vector length must equal"
  )
})

test_that("gc_sim_hierarchical_3d_per_zone returns list with data", {
  data <- create_3d_test_data()
  hz_def <- gc_define_hierarchy(
    data = data,
    comp_cols = c("SAND", "SILT", "CLAY"),
    group_col = "zone"
  )
  hz_fit <- gc_fit_hierarchical_model(data = data, object = hz_def)

  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2", "z3"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  coords <- data.frame(x = 1:3, y = 1:3, z = c(10, 50, 90))
  zone_vec <- data$zone[1:3]

  ensemble <- gc_sim_hierarchical_3d_per_zone(
    object = hz_fit,
    n = 5,
    nested_hierarchy = h_nested,
    coords = coords,
    zone_vector = zone_vec
  )

  expect_is(ensemble, "list")
  expect_true("data" %in% names(ensemble))
  expect_is(ensemble$data, "data.frame")
})

test_that("gc_sim_hierarchical_3d_per_zone respects n parameter", {
  data <- create_3d_test_data()
  hz_def <- gc_define_hierarchy(
    data = data,
    comp_cols = c("SAND", "SILT", "CLAY"),
    group_col = "zone"
  )
  hz_fit <- gc_fit_hierarchical_model(data = data, object = hz_def)

  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2", "z3"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  coords <- data.frame(x = 1:2, y = 1:2, z = c(10, 50))
  zone_vec <- c("z1", "z2")

  ensemble <- gc_sim_hierarchical_3d_per_zone(
    object = hz_fit,
    n = 20,
    nested_hierarchy = h_nested,
    coords = coords,
    zone_vector = zone_vec
  )

  expect_equal(ensemble$n_realizations, 20)
})

test_that("gc_sim_hierarchical_3d_per_zone SAND+SILT+CLAY close to 100", {
  data <- create_3d_test_data()
  hz_def <- gc_define_hierarchy(
    data = data,
    comp_cols = c("SAND", "SILT", "CLAY"),
    group_col = "zone"
  )
  hz_fit <- gc_fit_hierarchical_model(data = data, object = hz_def)

  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2", "z3"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  coords <- data.frame(x = 1:3, y = 1:3, z = c(10, 50, 90))
  zone_vec <- data$zone[1:3]

  ensemble <- gc_sim_hierarchical_3d_per_zone(
    object = hz_fit,
    n = 10,
    nested_hierarchy = h_nested,
    coords = coords,
    zone_vector = zone_vec
  )

  closure <- ensemble$data$SAND + ensemble$data$SILT + ensemble$data$CLAY
  expect_true(all(closure > 99, na.rm = TRUE))
  expect_true(all(closure < 101, na.rm = TRUE))
})

test_that("gc_sim_hierarchical_3d_per_zone depth_stratum assigned correctly", {
  data <- create_3d_test_data()
  hz_def <- gc_define_hierarchy(
    data = data,
    comp_cols = c("SAND", "SILT", "CLAY"),
    group_col = "zone"
  )
  hz_fit <- gc_fit_hierarchical_model(data = data, object = hz_def)

  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2", "z3"),
    depth_strata = c("shallow", "mid", "deep"),
    depth_breaks = c(0, 30, 60, 100)
  )

  coords <- data.frame(x = 1:3, y = 1:3, z = c(10, 45, 80))
  zone_vec <- data$zone[1:3]

  ensemble <- gc_sim_hierarchical_3d_per_zone(
    object = hz_fit,
    n = 5,
    nested_hierarchy = h_nested,
    coords = coords,
    zone_vector = zone_vec
  )

  # Check that depth strata are assigned
  expect_true("depth_stratum" %in% names(ensemble$data))
  expect_true(all(!is.na(ensemble$data$depth_stratum)))
})

test_that("gc_sim_hierarchical_3d_per_zone handles multiple zones", {
  data <- create_3d_test_data()
  hz_def <- gc_define_hierarchy(
    data = data,
    comp_cols = c("SAND", "SILT", "CLAY"),
    group_col = "zone"
  )
  hz_fit <- gc_fit_hierarchical_model(data = data, object = hz_def)

  h_nested <- gc_define_nested_hierarchy(
    domain_zones = c("z1", "z2", "z3"),
    depth_strata = c("shallow", "deep"),
    depth_breaks = c(0, 50, 100)
  )

  coords <- data.frame(x = 1:5, y = 1:5, z = c(10, 30, 50, 70, 90))
  zone_vec <- rep(c("z1", "z2", "z3", "z1", "z2"), 1)

  ensemble <- gc_sim_hierarchical_3d_per_zone(
    object = hz_fit,
    n = 10,
    nested_hierarchy = h_nested,
    coords = coords,
    zone_vector = zone_vec
  )

  expect_is(ensemble, "list")
  expect_equal(ensemble$n_locations, 5)
})
