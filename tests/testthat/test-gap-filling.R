# Gap Filling & Stitching Tests
# Tests for edge case handling in zone boundaries and raster merging

test_that("gap filling works with sparse data", {
  library(terra)
  
  # Create a raster with gaps (NAs)
  r <- terra::rast(nrows = 10, ncols = 10, resolution = 1)
  
  # Create data with strategic gaps
  vals <- c(1:50, NA, NA, NA, 51:97)
  r_sand <- r
  r_sand[] <- vals
  names(r_sand) <- "sand.1.sim1"
  
  # Create complementary layers
  r_silt <- r
  silt_vals <- c(rep(50, 50), rep(NA, 3), rep(40, 47))
  r_silt[] <- silt_vals
  names(r_silt) <- "silt.1.sim1"
  
  r_clay <- r
  clay_vals <- c(rep(50, 50), rep(NA, 3), rep(60, 47))
  r_clay[] <- clay_vals
  names(r_clay) <- "clay.1.sim1"
  
  r_stack <- c(r_sand, r_silt, r_clay)
  
  # Apply gap filling
  filled <- geocoda:::.fill_raster_gaps(r_stack, max_iterations = 3, verbose = FALSE)
  
  # Check that gaps are filled or reduced
  expect_true(methods::is(filled, "SpatRaster"))
  filled_vals <- terra::values(filled[[1]], na.rm = TRUE)
  expect_true(length(filled_vals) > length(terra::values(r_sand, na.rm = TRUE)))
})

test_that("rescaling enforces compositional sum constraints", {
  library(terra)
  
  # Create slightly invalid raster (sums ~ 95-105)
  r <- terra::rast(nrows = 5, ncols = 5, resolution = 1)
  
  r_sand <- r
  r_sand[] <- rnorm(terra::ncell(r), 20, 2)
  r_sand[] <- pmax(r_sand[], 0.1) # Ensure positive
  names(r_sand) <- "sand.1.sim1"
  
  r_silt <- r
  r_silt[] <- rnorm(terra::ncell(r), 40, 2)
  r_silt[] <- pmax(r_silt[], 0.1)
  names(r_silt) <- "silt.1.sim1"
  
  r_clay <- r
  r_clay[] <- rnorm(terra::ncell(r), 40, 2)
  r_clay[] <- pmax(r_clay[], 0.1)
  names(r_clay) <- "clay.1.sim1"
  
  r_stack <- c(r_sand, r_silt, r_clay)
  
  # Rescale to 100
  rescaled <- geocoda:::.rescale_compositions(r_stack, comp_layer_pattern = "\\.sim",
                                              target_sum = 100, tolerance = 0.5, verbose = FALSE)
  
  expect_true(methods::is(rescaled, "SpatRaster"))
  
  # Verify sums are now 100
  sands <- terra::values(rescaled[[1]], na.rm = TRUE)
  silts <- terra::values(rescaled[[2]], na.rm = TRUE)
  clays <- terra::values(rescaled[[3]], na.rm = TRUE)
  
  sums <- sands + silts + clays
  expect_true(all(sums >= 99.5 & sums <= 100.5))
})

test_that("zone boundary validation detects overlaps and gaps", {
  library(sf)
  
  # Create three zones: two with clear overlap, one with gap
  z1 <- sf::st_polygon(list(cbind(c(0,2,2,0,0), c(0,0,2,2,0))))
  z2 <- sf::st_polygon(list(cbind(c(1.5,3.5,3.5,1.5,1.5), c(0,0,2,2,0)))) # Overlaps with z1
  z3 <- sf::st_polygon(list(cbind(c(4,5,5,4,4), c(0,0,2,2,0)))) # Gap before this
  
  zones_sf <- sf::st_sf(
    zone_id = c(1, 2, 3),
    geometry = sf::st_sfc(z1, z2, z3),
    crs = NA
  )
  
  # Validate topology
  issues <- geocoda:::.validate_zone_topology(zones_sf, verbose = FALSE)
  
  # Should detect overlap and gap
  expect_true(length(issues) > 0)
  issue_text <- paste(issues, collapse = ",")
  expect_match(issue_text, "overlap|gap|coverage", perl = TRUE, ignore.case = TRUE)
})

test_that("buffer zone creates expanded boundaries correctly", {
  library(sf)
  library(terra)
  
  # Create simple zone
  z1 <- sf::st_polygon(list(cbind(c(0,1,1,0,0), c(0,0,1,1,0))))
  zones_sf <- sf::st_sf(
    zone_id = 1,
    geometry = sf::st_sfc(z1),
    crs = NA
  )
  
  # Buffer with fixed distance
  buffered <- geocoda:::.buffer_zone(zones_sf[1, ], buffer_distance = 0.2, method = "fixed")
  
  # Check that buffered area is larger
  orig_area <- sf::st_area(zones_sf[1, ])
  buff_area <- sf::st_area(buffered)
  
  expect_true(buff_area > orig_area)
})

test_that("zone raster merging handles overlaps correctly", {
  library(terra)
  
  # Create 3 rasters with different values in overlapping regions
  r_base <- terra::rast(nrows = 10, ncols = 10, resolution = 1)
  
  # Zone 1: left half
  r1 <- r_base
  r1[1:100] <- 10 # Will be overwritten by overlap
  r1[1:50] <- NA
  r1[51:100] <- 10 # Keep these
  
  # Zone 2: right half (overlaps with zone 1)
  r2 <- r_base
  r2[] <- NA
  r2[51:75] <- 20 # Overlap zone
  r2[76:100] <- 20 # Non-overlap zone
  
  zone_rasters <- list(zone1 = r1, zone2 = r2)
  
  # Merge with "buffer" strategy
  merged <- geocoda:::.merge_zone_rasters(zone_rasters, r_base, overlap_strategy = "buffer")
  
  # Verify merge completed
  expect_true(methods::is(merged, "SpatRaster"))
  
  # Check that some cells have data
  merged_vals <- terra::values(merged, na.rm = TRUE)
  expect_true(length(merged_vals) > 0)
})
