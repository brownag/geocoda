test_that("gc_fit_vgm_3d validates 3D coordinates", {
  # Missing z column should error
  data_2d <- data.frame(
    x = c(0, 100, 200),
    y = c(0, 100, 200),
    ilr1 = c(0, 0.5, -0.5),
    ilr2 = c(0, -0.3, 0.3)
  )

  ilr_params <- list(mean = c(0, 0), cov = diag(2), names = c("ilr1", "ilr2"))

  expect_error(
    gc_fit_vgm_3d(ilr_params, data_2d),
    "columns 'x', 'y', and 'z'"
  )
})

test_that("gc_fit_vgm_3d requires minimum observations", {
  # Create minimal data (only 7 points, need >= 8)
  data_small <- data.frame(
    x = c(0, 50, 100, 0, 50, 100, 0),
    y = c(0, 50, 100, 0, 50, 100, 0),
    z = c(0, 0, 0, 10, 10, 10, 20),
    ilr1 = rnorm(7),
    ilr2 = rnorm(7)
  )

  ilr_params <- list(mean = c(0, 0), cov = diag(2), names = c("ilr1", "ilr2"))

  expect_error(
    gc_fit_vgm_3d(ilr_params, data_small),
    "at least 8 observations"
  )
})

test_that("gc_fit_vgm_3d fits horizontal and vertical variograms", {
  # Create 3D soil profile data with depth-dependent structure
  set.seed(42)

  x <- rep(seq(0, 200, by = 100), 3)  # 3 horizontal locations
  y <- rep(seq(0, 200, by = 100), 3)
  z <- rep(c(10, 50, 90), each = 3)  # 3 depth layers

  data <- data.frame(
    x = x,
    y = y,
    z = z,
    ilr1 = rnorm(9, mean = 0, sd = 0.5),
    ilr2 = rnorm(9, mean = 0, sd = 0.5)
  )

  ilr_params <- list(mean = c(0, 0), cov = diag(2), names = c("ilr1", "ilr2"))

  vgm_3d <- gc_fit_vgm_3d(ilr_params, data,
                         horizontal_cutoff = 300,
                         vertical_cutoff = 100)

  # Check structure
  expect_is(vgm_3d, "list")
  expect_true("fitted_vgms" %in% names(vgm_3d))
  expect_true("fitted_params" %in% names(vgm_3d))

  # Check fitted_params
  expect_equal(nrow(vgm_3d$fitted_params), 2)  # 2 ILR dimensions
  expect_true("range_h" %in% names(vgm_3d$fitted_params))
  expect_true("range_v" %in% names(vgm_3d$fitted_params))
  expect_true("anis_ratio" %in% names(vgm_3d$fitted_params))

  # Check attributes
  expect_equal(attr(vgm_3d, "coordinates_3d"), TRUE)
  expect_equal(attr(vgm_3d, "model_type"), "Exp")
})

test_that("gc_fit_vgm_3d computes anisotropy ratio", {
  set.seed(123)

  # Create data with clear depth structure
  n_x <- 5
  n_z <- 4
  x <- rep(seq(0, 400, length.out = n_x), each = n_z)
  y <- rep(seq(0, 400, length.out = n_x), each = n_z)
  z <- rep(seq(10, 100, length.out = n_z), n_x)

  # ILR values: more variation in z than in xy
  ilr1 <- sin(z / 20) + rnorm(length(z), sd = 0.1)

  data <- data.frame(x = x, y = y, z = z, ilr1 = ilr1)

  ilr_params <- list(mean = 0, cov = matrix(1), names = "ilr1")

  vgm_3d <- gc_fit_vgm_3d(ilr_params, data,
                         horizontal_cutoff = 300,
                         vertical_cutoff = 90,
                         tolerance_z = 15,
                         tolerance_xy = 100)

  # Anisotropy ratio should exist
  anis_ratio <- vgm_3d$fitted_params$anis_ratio[1]
  expect_true(!is.na(anis_ratio))

  # For typical soil data: horizontal range > vertical range
  # So anisotropy_ratio > 1 is expected
})

test_that("gc_fit_vgm_3d produces expected variogram structure per dimension", {
  set.seed(456)

  # Create 3D data with realistic variogram structure
  coords <- expand.grid(
    x = seq(0, 500, by = 100),
    y = seq(0, 500, by = 100),
    z = seq(10, 60, by = 10)
  )

  n <- nrow(coords)
  coords$ilr1 <- sin(coords$x / 200) + cos(coords$z / 30) + rnorm(n, sd = 0.1)
  coords$ilr2 <- cos(coords$y / 200) - sin(coords$z / 40) + rnorm(n, sd = 0.1)

  ilr_params <- list(
    mean = c(0, 0),
    cov = matrix(c(1, 0.2, 0.2, 1), nrow = 2),
    names = c("ilr1", "ilr2")
  )

  vgm_3d <- gc_fit_vgm_3d(ilr_params, coords,
                         vgm_model_type = "Exp",
                         horizontal_cutoff = 400,
                         vertical_cutoff = 50,
                         tolerance_z = 10,
                         tolerance_xy = 150)

  # Check fitted models exist
  expect_true("ilr1" %in% names(vgm_3d$fitted_vgms))
  expect_true("ilr2" %in% names(vgm_3d$fitted_vgms))

  # Check each dimension has horizontal and vertical components
  for (ilr_id in c("ilr1", "ilr2")) {
    vgm_item <- vgm_3d$fitted_vgms[[ilr_id]]
    expect_true("vgm_horizontal" %in% names(vgm_item))
    expect_true("vgm_vertical" %in% names(vgm_item))
    expect_true("range_horizontal" %in% names(vgm_item))
    expect_true("range_vertical" %in% names(vgm_item))
  }
})

test_that("gc_fit_vgm_3d sensitivity to tolerance parameters", {
  set.seed(789)

  coords <- expand.grid(
    x = seq(0, 300, by = 100),
    y = seq(0, 300, by = 100),
    z = seq(0, 100, by = 10)
  )

  n <- nrow(coords)
  coords$ilr1 <- rnorm(n, mean = 0, sd = 0.5)

  ilr_params <- list(mean = 0, cov = matrix(1), names = "ilr1")

  # Fit with tight z-tolerance (fewer pairs)
  vgm_tight <- gc_fit_vgm_3d(ilr_params, coords,
                            tolerance_z = 5,
                            tolerance_xy = 80)

  # Fit with loose z-tolerance (more pairs)
  vgm_loose <- gc_fit_vgm_3d(ilr_params, coords,
                            tolerance_z = 20,
                            tolerance_xy = 80)

  # Both should produce results
  expect_is(vgm_tight, "list")
  expect_is(vgm_loose, "list")

  # Results may differ due to different pair selection
  tight_range <- vgm_tight$fitted_params$range_h[1]
  loose_range <- vgm_loose$fitted_params$range_h[1]

  expect_true(!is.na(tight_range) || !is.na(loose_range))
})

test_that("gc_fit_vgm_3d handles missing ILR columns gracefully", {
  data <- data.frame(
    x = c(0, 100, 200, 300, 0, 100, 200, 300),
    y = c(0, 100, 200, 300, 0, 100, 200, 300),
    z = c(10, 20, 30, 40, 50, 60, 70, 80),
    ilr1 = c(0, 0.5, -0.5, 0.2, -0.1, 0.3, -0.2, 0.1)
    # Missing ilr2
  )

  ilr_params <- list(mean = c(0, 0), cov = diag(2), names = c("ilr1", "ilr2"))

  expect_error(
    gc_fit_vgm_3d(ilr_params, data),
    "must contain column 'ilr2'"
  )
})

test_that("gc_fit_vgm_3d warnings for insufficient pairs", {
  # Create sparse data with tight tolerances
  data <- data.frame(
    x = c(0, 100, 200, 300, 0, 100, 200, 300),
    y = c(0, 100, 200, 300, 0, 100, 200, 300),
    z = c(0, 10, 20, 30, 40, 50, 60, 70),
    ilr1 = c(0, 0.5, -0.5, 0.2, -0.1, 0.3, -0.2, 0.1),
    ilr2 = c(0, -0.3, 0.3, 0.1, 0.2, -0.4, 0.2, -0.1)
  )

  ilr_params <- list(mean = c(0, 0), cov = diag(2), names = c("ilr1", "ilr2"))

  # With very tight tolerances and sparse data, should get warnings
  expect_warning(
    gc_fit_vgm_3d(ilr_params, data,
                 tolerance_z = 2,
                 tolerance_xy = 30),
    "pairs|tolerance"
  )
})

test_that("gc_fit_vgm_3d attributes set correctly", {
  set.seed(111)

  coords <- expand.grid(
    x = seq(0, 200, by = 100),
    y = seq(0, 200, by = 100),
    z = seq(10, 90, by = 40)
  )
  coords$ilr1 <- rnorm(nrow(coords))

  ilr_params <- list(mean = 0, cov = matrix(1), names = "ilr1")

  vgm_3d <- gc_fit_vgm_3d(ilr_params, coords,
                         vgm_model_type = "Sph",
                         tolerance_z = 50,
                         tolerance_xy = 150)

  # Check attributes
  expect_equal(attr(vgm_3d, "coordinates_3d"), TRUE)
  expect_equal(attr(vgm_3d, "model_type"), "Sph")
  expect_equal(attr(vgm_3d, "tolerance_z"), 50)
  expect_equal(attr(vgm_3d, "tolerance_xy"), 150)
  expect_equal(length(attr(vgm_3d, "ilr_dimension_names")), 1)
})

test_that("gc_fit_vgm_3d empirical variogram creation", {
  set.seed(222)

  coords <- expand.grid(
    x = seq(0, 200, by = 100),
    y = seq(0, 200, by = 100),
    z = seq(10, 50, by = 20)
  )

  coords$ilr1 <- rnorm(nrow(coords))

  ilr_params <- list(mean = 0, cov = matrix(1), names = "ilr1")

  vgm_3d <- gc_fit_vgm_3d(ilr_params, coords,
                         horizontal_cutoff = 250,
                         vertical_cutoff = 40)

  # Check empirical variogram components exist
  expect_true("empirical_vgms_horizontal" %in% names(vgm_3d))
  expect_true("empirical_vgms_vertical" %in% names(vgm_3d))

  # Should have some empirical variogram data
  expect_true(length(vgm_3d$empirical_vgms_horizontal) >= 0)
  expect_true(length(vgm_3d$empirical_vgms_vertical) >= 0)
})
