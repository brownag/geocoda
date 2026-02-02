test_that("gc_define_surfaces with constant surfaces", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  expect_s4_class(surfaces$top, "SpatRaster")
  expect_s4_class(surfaces$bottom, "SpatRaster")
  expect_s4_class(surfaces$thickness, "SpatRaster")
  expect_equal(surfaces$type, "constant")

  # Check values
  expect_equal(mean(terra::values(surfaces$top), na.rm = TRUE), 0)
  expect_equal(mean(terra::values(surfaces$bottom), na.rm = TRUE), -100)
  expect_equal(mean(terra::values(surfaces$thickness), na.rm = TRUE), 100)
})

test_that("gc_define_surfaces validates top > bottom", {
  expect_error(
    gc_define_surfaces(
      top_surface = -100,
      bottom_surface = 0,
      type = "constant",
      extent = terra::ext(c(0, 1000, 0, 1000)),
      resolution = c(10, 10)
    ),
    "top must be > bottom"
  )
})

test_that("gc_define_surfaces with raster surfaces", {
  # Create test rasters
  extent <- terra::ext(c(0, 100, 0, 100))
  top_rast <- terra::rast(extent, res = 10, vals = 0)
  bottom_rast <- terra::rast(extent, res = 10, vals = -50)

  surfaces <- gc_define_surfaces(
    top_surface = top_rast,
    bottom_surface = bottom_rast,
    type = "raster"
  )

  expect_s4_class(surfaces$top, "SpatRaster")
  expect_s4_class(surfaces$thickness, "SpatRaster")
  expect_equal(mean(terra::values(surfaces$thickness), na.rm = TRUE), 50)
})

test_that("gc_stratigraphic_transform proportional method", {
  # Create test surfaces
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  # Create test data at different depths
  data <- data.frame(
    x = c(100, 200, 300),
    y = c(100, 200, 300),
    z = c(0, -50, -100)  # top, middle, bottom
  )

  # Transform using proportional method
  data_strat <- gc_stratigraphic_transform(
    data, surfaces,
    transform_type = "proportional",
    h_const = 100
  )

  # Check results
  expect_equal(data_strat$z_strat[1], 0)     # At top
  expect_equal(data_strat$z_strat[2], 50)    # At middle
  expect_equal(data_strat$z_strat[3], 100)   # At bottom
  expect_true(all(c("x", "y", "z", "z_strat") %in% names(data_strat)))

  # Check metadata
  expect_equal(attr(data_strat, "transform_type"), "proportional")
  expect_equal(attr(data_strat, "h_const"), 100)
})

test_that("gc_stratigraphic_transform top_onlap method", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data <- data.frame(
    x = c(100, 200, 300),
    y = c(100, 200, 300),
    z = c(0, -50, -100)
  )

  data_strat <- gc_stratigraphic_transform(
    data, surfaces,
    transform_type = "top_onlap"
  )

  expect_equal(data_strat$z_strat[1], 0)    # At surface
  expect_equal(data_strat$z_strat[2], 50)   # 50 cm below surface
  expect_equal(data_strat$z_strat[3], 100)  # 100 cm below surface
})

test_that("gc_stratigraphic_transform bottom_erosional method", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data <- data.frame(
    x = c(100, 200, 300),
    y = c(100, 200, 300),
    z = c(0, -50, -100)
  )

  data_strat <- gc_stratigraphic_transform(
    data, surfaces,
    transform_type = "bottom_erosional"
  )

  expect_equal(data_strat$z_strat[1], 100)  # 100 cm above base
  expect_equal(data_strat$z_strat[2], 50)   # 50 cm above base
  expect_equal(data_strat$z_strat[3], 0)    # At base
})

test_that("gc_stratigraphic_transform rotational method", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data <- data.frame(
    x = c(100, 100, 100),
    y = c(100, 100, 100),
    z = c(0, -50, -100)
  )

  data_strat <- gc_stratigraphic_transform(
    data, surfaces,
    transform_type = "rotational",
    dip = 0,
    azimuth = 0
  )

  # With 0 dip, should have x_strat, y_strat, z_strat
  expect_true("x_strat" %in% names(data_strat))
  expect_true("y_strat" %in% names(data_strat))
  expect_true("z_strat" %in% names(data_strat))
})

test_that("gc_inverse_stratigraphic_transform proportional round-trip", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data_orig <- data.frame(
    x = c(100, 200, 300),
    y = c(100, 200, 300),
    z = c(-10, -50, -90)
  )

  # Forward transformation
  data_strat <- gc_stratigraphic_transform(
    data_orig, surfaces,
    transform_type = "proportional"
  )

  # Inverse transformation
  data_back <- gc_inverse_stratigraphic_transform(data_strat)

  # Check round-trip
  expect_equal(data_back$z, data_orig$z, tolerance = 1e-6)
  expect_equal(data_back$x, data_orig$x, tolerance = 1e-6)
  expect_equal(data_back$y, data_orig$y, tolerance = 1e-6)
})

test_that("gc_inverse_stratigraphic_transform top_onlap round-trip", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data_orig <- data.frame(
    x = c(100, 200, 300),
    y = c(100, 200, 300),
    z = c(-10, -50, -90)
  )

  data_strat <- gc_stratigraphic_transform(
    data_orig, surfaces,
    transform_type = "top_onlap"
  )

  data_back <- gc_inverse_stratigraphic_transform(data_strat)

  expect_equal(data_back$z, data_orig$z, tolerance = 1e-6)
})

test_that("gc_inverse_stratigraphic_transform bottom_erosional round-trip", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data_orig <- data.frame(
    x = c(100, 200, 300),
    y = c(100, 200, 300),
    z = c(-10, -50, -90)
  )

  data_strat <- gc_stratigraphic_transform(
    data_orig, surfaces,
    transform_type = "bottom_erosional"
  )

  data_back <- gc_inverse_stratigraphic_transform(data_strat)

  expect_equal(data_back$z, data_orig$z, tolerance = 1e-6)
})

test_that("gc_stratigraphic_transform warns for points outside unit", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data <- data.frame(
    x = c(100, 200),
    y = c(100, 200),
    z = c(50, -150)  # Second point is outside unit
  )

  expect_warning(
    gc_stratigraphic_transform(data, surfaces, transform_type = "proportional", validate = TRUE),
    "points outside stratigraphic unit"
  )
})

test_that("gc_stratigraphic_transform with no validation", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data <- data.frame(
    x = c(100, 200),
    y = c(100, 200),
    z = c(50, -150)
  )

  # Should not warn when validate = FALSE
  result <- gc_stratigraphic_transform(
    data, surfaces,
    transform_type = "proportional",
    validate = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true("z_strat" %in% names(result))
})

test_that("gc_stratigraphic_transform preserves other columns", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data <- data.frame(
    x = c(100, 200),
    y = c(100, 200),
    z = c(-50, -50),
    sample_id = c("A", "B"),
    value = c(1.5, 2.3)
  )

  data_strat <- gc_stratigraphic_transform(data, surfaces, transform_type = "proportional")

  expect_true("sample_id" %in% names(data_strat))
  expect_true("value" %in% names(data_strat))
  expect_equal(data_strat$sample_id, c("A", "B"))
  expect_equal(data_strat$value, c(1.5, 2.3))
})

test_that("gc_inverse_stratigraphic_transform requires metadata", {
  data_no_meta <- data.frame(
    x = c(100, 200),
    y = c(100, 200),
    z_strat = c(50, 50)
  )

  expect_error(
    gc_inverse_stratigraphic_transform(data_no_meta),
    "must have transformation metadata"
  )
})

test_that("gc_stratigraphic_transform with variable thickness surfaces", {
  # Create rasters with varying surfaces
  extent <- terra::ext(c(0, 100, 0, 100))
  top_rast <- terra::rast(extent, res = 10)
  bottom_rast <- terra::rast(extent, res = 10)

  # Top surface varies from 0 to 10
  top_vals <- seq(0, 10, length.out = 100)
  terra::values(top_rast) <- top_vals

  # Bottom surface varies from -90 to -110
  bottom_vals <- seq(-90, -110, length.out = 100)
  terra::values(bottom_rast) <- bottom_vals

  surfaces <- gc_define_surfaces(
    top_surface = top_rast,
    bottom_surface = bottom_rast,
    type = "raster"
  )

  # Create points at different locations
  data <- data.frame(
    x = c(5, 15, 25),
    y = c(5, 15, 25),
    z = c(-50, -50, -50)  # Same relative depth
  )

  data_strat <- gc_stratigraphic_transform(
    data, surfaces,
    transform_type = "proportional"
  )

  # All should have similar z_strat values (around 50)
  expect_true(all(data_strat$z_strat > 40 & data_strat$z_strat < 60))
})

test_that("gc_define_surfaces input validation", {
  # Missing required parameters for constant type
  expect_error(
    gc_define_surfaces(0, -100, type = "constant"),
    "extent|resolution"
  )

  # Invalid type
  expect_error(
    gc_define_surfaces(0, -100, type = "invalid"),
    "not TRUE"
  )

  # Non-numeric constant
  expect_error(
    gc_define_surfaces("not_a_number", -100, type = "constant"),
    "is.numeric"
  )
})

test_that("gc_stratigraphic_transform input validation", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  # Missing required columns
  data_no_z <- data.frame(x = c(100, 200), y = c(100, 200))
  expect_error(
    gc_stratigraphic_transform(data_no_z, surfaces),
    "x.*y.*z"
  )

  # Invalid transform_type
  data <- data.frame(x = c(100, 200), y = c(100, 200), z = c(-50, -50))
  expect_error(
    gc_stratigraphic_transform(data, surfaces, transform_type = "invalid"),
    "transform_type"
  )

  # Rotational without dip/azimuth
  expect_error(
    gc_stratigraphic_transform(data, surfaces, transform_type = "rotational"),
    "dip|azimuth"
  )
})

test_that("gc_stratigraphic_transform stores metadata correctly", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data <- data.frame(
    x = c(100, 200),
    y = c(100, 200),
    z = c(-50, -50)
  )

  data_strat <- gc_stratigraphic_transform(
    data, surfaces,
    transform_type = "proportional",
    h_const = 150
  )

  # Check metadata
  expect_equal(attr(data_strat, "transform_type"), "proportional")
  expect_equal(attr(data_strat, "h_const"), 150)
  expect_equal(class(attr(data_strat, "surfaces")), "list")
  expect_true("top" %in% names(attr(data_strat, "surfaces")))
})

test_that("gc_stratigraphic_transform with rotational stores rotation params", {
  surfaces <- gc_define_surfaces(
    top_surface = 0,
    bottom_surface = -100,
    type = "constant",
    extent = terra::ext(c(0, 1000, 0, 1000)),
    resolution = c(10, 10)
  )

  data <- data.frame(
    x = c(100, 200),
    y = c(100, 200),
    z = c(-50, -50)
  )

  data_strat <- gc_stratigraphic_transform(
    data, surfaces,
    transform_type = "rotational",
    dip = 15,
    azimuth = 45
  )

  expect_equal(attr(data_strat, "dip"), 15)
  expect_equal(attr(data_strat, "azimuth"), 45)
})
