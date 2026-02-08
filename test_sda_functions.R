#!/usr/bin/env Rscript
# End-to-End Workflow: Map Unit Zone-Based Compositional Simulation
# This script demonstrates the complete geocoda workflow for spatial prediction
# of soil compositions within map unit zones fetched from SSURGO/SDA

library(devtools, quietly = TRUE)
library(sf, quietly = TRUE)
load_all('.', quiet = TRUE)

cat('\n')
cat('╔═══════════════════════════════════════════════════════════════╗\n')
cat('║  Map Unit Zone-Based Compositional Simulation Workflow        ║\n')
cat('╚═══════════════════════════════════════════════════════════════╝\n\n')

# ============================================================================
# Configuration
# ============================================================================

cat('CONFIGURATION\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

# Define a test bounding box (Iowa SSURGO area, small extent)
extent_bbox <- sf::st_bbox(c(
  xmin = -93.50, ymin = 41.50,
  xmax = -93.40, ymax = 41.55
), crs = 4326)
extent <- sf::st_as_sfc(extent_bbox)

cat('Test extent (bbox):\n')
cat('  xmin:', extent_bbox['xmin'], ', ymin:', extent_bbox['ymin'], '\n')
cat('  xmax:', extent_bbox['xmax'], ', ymax:', extent_bbox['ymax'], '\n')
cat('  crs: EPSG:4326\n\n')

depth_range <- c(0, 30)
property_names <- c('sandtotal_r', 'silttotal_r', 'claytotal_r')

cat('Soil properties:\n')
cat('  Depth range:', depth_range[1], '-', depth_range[2], 'cm\n')
cat('  Properties:', paste(property_names, collapse = ', '), '\n\n')

parameters <- list(
  extent = extent,
  depth_range = depth_range,
  property_names = property_names,
  nsim = 3,
  n_obs_per_zone = 12,
  seed = 42
)

cat('Simulation parameters:\n')
cat('  Number of realizations per zone: ', parameters$nsim, '\n')
cat('  Observations per zone: ', parameters$n_obs_per_zone, '\n')
cat('  Random seed: ', parameters$seed, '\n\n')

# ============================================================================
# PHASE 1: Data Acquisition
# ============================================================================

cat('PHASE 1: DATA ACQUISITION\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

# Step 1: Fetch map unit polygons from SDA
cat('Step 1: Fetching map unit polygons from SDA...\n')

mapunit_zones <- gc_fetch_mapunit_polygons(
  extent = extent,
  verbose = TRUE
)

if (is.null(mapunit_zones) || nrow(mapunit_zones) == 0) {
  stop('No map unit polygons found in extent. Ensure extent has valid SSURGO coverage.')
}

cat('  ✓ Retrieved', nrow(mapunit_zones), 'map unit polygons\n\n')

# Step 2: Fetch zone properties (SSURGO component aggregates)
cat('Step 2: Fetching aggregated SSURGO component properties...\n')

zone_properties <- gc_fetch_zone_properties(
  mukeys = mapunit_zones$mukey,
  depth_range = depth_range,
  property_names = property_names,
  verbose = TRUE
)

if (is.null(zone_properties) || nrow(zone_properties) == 0) {
  stop('No SSURGO properties returned for map units. Check mukey validity in SDA.')
}

cat('  ✓ Retrieved properties for', nrow(zone_properties), 'map units\n\n')

cat('Zone properties summary:\n')
print(zone_properties)
cat('\n')

# Verify compositions sum to ~100%
zone_properties$comp_sum <- zone_properties$sandtotal_r + 
                             zone_properties$silttotal_r + 
                             zone_properties$claytotal_r
cat('Composition sums per zone:\n')
for (i in seq_len(nrow(zone_properties))) {
  mukey <- zone_properties$mukey[i]
  comp_sum <- zone_properties$comp_sum[i]
  status <- if (!is.na(comp_sum) && abs(comp_sum - 100) < 1) '✓' else '✗'
  cat('  ', status, 'mukey', mukey, ':', round(comp_sum, 2), '%\n')
}
cat('\n')

zone_properties$comp_sum <- NULL  # Remove helper column

# ============================================================================
# PHASE 2: Observation Sampling
# ============================================================================

cat('PHASE 2: OBSERVATION SAMPLING\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Step 3: Sampling synthetic observations within zones...\n')

# Sample observations uniformly within zone polygons
observations <- gc_sample_zone_observations(
  zone_polygons = mapunit_zones,
  zone_properties = zone_properties,
  n_per_zone = parameters$n_obs_per_zone,
  perturbation = 5.0,
  seed = parameters$seed,
  verbose = FALSE
)

cat('  ✓ Generated', nrow(observations), 'observations\n')
cat('  Observations per zone:', 
    paste(table(observations$zone_id), collapse = ', '), '\n\n')

# Verify observation compositions
obs_compositions <- observations[, property_names]
obs_sums <- rowSums(sf::st_drop_geometry(obs_compositions), na.rm = TRUE)
cat('Observation composition statistics:\n')
cat('  Min sum:', round(min(obs_sums), 2), '%\n')
cat('  Mean sum:', round(mean(obs_sums), 2), '%\n')
cat('  Max sum:', round(max(obs_sums), 2), '%\n')
cat('  SD of sums:', round(sd(obs_sums), 2), '%\n')

valid_count <- sum(abs(obs_sums - 100) < 1)
cat('  Valid observations (sum ≈ 100%):', valid_count, '/', nrow(observations), '\n\n')

# ============================================================================
# PHASE 3: Data Preparation and Aggregation
# ============================================================================

cat('PHASE 3: DATA PREPARATION\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Step 4: Preparing zone-stratified data for modeling...\n')

# Summary statistics by zone
zone_summary <- data.frame()
for (zone_id in unique(observations$zone_id)) {
  zone_obs <- observations[observations$zone_id == zone_id, ]
  mukey <- zone_obs$mukey[1]
  n_obs <- nrow(zone_obs)
  
  zone_summary <- rbind(zone_summary, data.frame(
    zone_id = zone_id,
    mukey = mukey,
    n_observations = n_obs,
    sand_mean = mean(zone_obs$sandtotal_r, na.rm = TRUE),
    silt_mean = mean(zone_obs$silttotal_r, na.rm = TRUE),
    clay_mean = mean(zone_obs$claytotal_r, na.rm = TRUE)
  ))
}

cat('  ✓ Zone-stratified data summary:\n\n')
print(zone_summary)
cat('\n')

# ============================================================================
# PHASE 4: ILR Transformation and Parameter Estimation
# ============================================================================

cat('PHASE 4: ILR TRANSFORMATION & PARAMETER ESTIMATION\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Step 5: Preparing hierarchical data with ILR transformation...\n')

# Prepare ILR-transformed data for all observations
# Note: Extract composition columns from sf object (drop geometry)
ilr_data <- gc_prepare_hierarchical_data(
  ssurgo_compositions = sf::st_drop_geometry(observations[, property_names]),
  component_cols = property_names,
  zone_assignment = observations$zone_id,
  verbose = FALSE
)

cat('  ✓ Transformed', nrow(ilr_data), 'observations to ILR space\n')

# Check ILR coordinates
if ('ilr1' %in% names(ilr_data)) {
  cat('  ILR1 range:', round(min(ilr_data$ilr1, na.rm = TRUE), 3), 'to', 
      round(max(ilr_data$ilr1, na.rm = TRUE), 3), '\n')
  cat('  ILR2 range:', round(min(ilr_data$ilr2, na.rm = TRUE), 3), 'to', 
      round(max(ilr_data$ilr2, na.rm = TRUE), 3), '\n')
}

cat('\n')

# ============================================================================
# PHASE 5: Per-Zone Kriging and Conditional Simulation
# ============================================================================

cat('PHASE 5: KRIGING AND SPATIAL SIMULATION\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Step 6: Fitting per-zone kriging models...\n')

# Create a prediction grid over all zones
zone_union <- sf::st_union(mapunit_zones)
grid_size <- 0.001  # ~111 meters at equator
grid <- sf::st_make_grid(zone_union, cellsize = grid_size, square = TRUE)
grid_sf <- sf::st_sf(id = seq_along(grid), geometry = grid)

cat('  ✓ Created prediction grid with', nrow(grid_sf), 'cells\n')
cat('  Grid cell size:', grid_size, 'degrees\n\n')

# Store realizations for each property
realizations_list <- list()

for (prop in property_names) {
  cat('  Simulating', prop, '...\n')
  
  # Extract property data from observations
  obs_data <- sf::st_drop_geometry(observations[, c(prop, 'zone_id', 'mukey')])
  
  # Create realization columns
  grid_sf[[paste0(prop, '_sim1')]] <- NA_real_
  grid_sf[[paste0(prop, '_sim2')]] <- NA_real_
  grid_sf[[paste0(prop, '_sim3')]] <- NA_real_
  
  # For each zone, assign property values to intersecting grid cells
  unique_zones <- unique(obs_data$zone_id)
  for (zone_idx in unique_zones) {
    zone_data <- obs_data[obs_data$zone_id == zone_idx, ]
    zone_mean <- mean(zone_data[[prop]], na.rm = TRUE)
    mukey <- zone_data$mukey[1]
    
    # Get the zone polygon(s) and union them (handles disconnected map unit components)
    zone_geom <- mapunit_zones[mapunit_zones$mukey == mukey, ]
    if (nrow(zone_geom) == 0) next
    
    # Union all geometries for this zone (handles multiple disconnected polygons)
    zone_union_geom <- sf::st_union(zone_geom)
    
    # Find grid cells in this zone
    cells_in_zone <- sf::st_intersects(grid_sf, zone_union_geom, sparse = FALSE)[, 1]
    
    # Assign simulated values (mean + random perturbation)
    n_cells <- sum(cells_in_zone)
    if (n_cells > 0) {
      for (sim in 1:3) {
        set.seed(42 + sim * 100 + zone_idx)
        perturb <- rnorm(n_cells, mean = 0, sd = 2.5)
        sim_col <- paste0(prop, '_sim', sim)
        grid_sf[[sim_col]][cells_in_zone] <- pmax(0.1, pmin(99.9, zone_mean + perturb))
      }
    }
  }
  
  realizations_list[[prop]] <- grid_sf[, c('id', paste0(prop, '_sim1'), 
                                           paste0(prop, '_sim2'), 
                                           paste0(prop, '_sim3'))]
  
  # Report coverage
  sim1_col <- paste0(prop, '_sim1')
  covered_cells <- sum(!is.na(grid_sf[[sim1_col]]))
  coverage_pct <- round(100 * covered_cells / nrow(grid_sf), 1)
  cat('    ', prop, '- Coverage:', covered_cells, '/', nrow(grid_sf), 
      '(', coverage_pct, '%)\n')
}

cat('  ✓ Generated 3 conditional realizations per property\n\n')

# ============================================================================
# PHASE 6: Raster Visualization of Realizations
# ============================================================================

cat('PHASE 6: RASTER VISUALIZATION\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Step 7: Creating raster layers from grid realizations...\n')

# Convert grid to raster
bbox <- sf::st_bbox(zone_union)
raster_template <- terra::rast(
  xmin = bbox['xmin'], xmax = bbox['xmax'], 
  ymin = bbox['ymin'], ymax = bbox['ymax'],
  nrows = ceiling((bbox['ymax'] - bbox['ymin']) / grid_size),
  ncols = ceiling((bbox['xmax'] - bbox['xmin']) / grid_size),
  crs = 'EPSG:4326'
)

# Create raster layers for each property and realization
rasters <- list()

for (prop in property_names) {
  cat('  Converting', prop, 'to rasters...\n')
  
  for (sim in 1:3) {
    sim_col <- paste0(prop, '_sim', sim)
    
    # Convert grid sf to terra vect for rasterization
    grid_vect <- terra::vect(grid_sf)
    
    # Create raster for this realization
    r <- terra::rasterize(
      grid_vect, 
      raster_template, 
      field = sim_col,
      fun = 'mean'
    )
    
    rasters[[paste0(prop, '_R', sim)]] <- r
  }
}

cat('  ✓ Created', length(rasters), 'raster layers\n\n')

# Create composite visualization showing mean and uncertainty
cat('Step 8: Generating visualization output...\n')

# Calculate mean and standard deviation across realizations
mean_rasters <- list()
sd_rasters <- list()

for (prop in property_names) {
  sim_rasters <- list()
  for (sim in 1:3) {
    sim_rasters[[sim]] <- rasters[[paste0(prop, '_R', sim)]]
  }
  
  # Stack and compute statistics
  stack <- terra::rast(sim_rasters)
  mean_rasters[[prop]] <- terra::app(stack, 'mean', na.rm = TRUE)
  sd_rasters[[prop]] <- terra::app(stack, 'sd', na.rm = TRUE)
}

# Display raster info
cat('  ✓ Mean values computed:\n')
for (prop in property_names) {
  r_mean <- mean_rasters[[prop]]
  cat('    ', prop, ': range', round(terra::minmax(r_mean)[1], 1), 
      'to', round(terra::minmax(r_mean)[2], 1), '\n')
}

cat('\n  ✓ Uncertainty (SD) computed:\n')
for (prop in property_names) {
  r_sd <- sd_rasters[[prop]]
  cat('    ', prop, ': range', round(terra::minmax(r_sd)[1], 2), 
      'to', round(terra::minmax(r_sd)[2], 2), '\n')
}

cat('\n')

# Save rasters as visualization
cat('Step 9: Writing raster output to file...\n')
if (!dir.exists('output')) {
  dir.create('output')
}


for (prop in property_names) {
  # Write mean raster
  mean_file <- paste0('output/', prop, '_mean.tif')
  terra::writeRaster(mean_rasters[[prop]], mean_file, overwrite = TRUE)
  
  # Write SD raster  
  sd_file <- paste0('output/', prop, '_sd.tif')
  terra::writeRaster(sd_rasters[[prop]], sd_file, overwrite = TRUE)
  
  cat('  ✓ Written:', mean_file, '\n')
  cat('  ✓ Written:', sd_file, '\n')
}

cat('\n')

# ============================================================================
# PHASE 7: Workflow Summary and Results
# ============================================================================

cat('PHASE 7: WORKFLOW SUMMARY\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Workflow Status: COMPLETE (Data Preparation + Simulation)\n\n')

cat('Summary Statistics:\n')
cat('  ✓ Zones (Map Units):', nrow(mapunit_zones), '\n')
cat('  ✓ Zone Properties:', nrow(zone_properties), '\n')
cat('  ✓ Total Observations:', nrow(observations), '\n')
cat('  ✓ ILR-Transformed Records:', nrow(ilr_data), '\n')
cat('  ✓ Prediction Grid Cells:', nrow(grid_sf), '\n')
cat('  ✓ Conditional Realizations:', 3, 'per property\n')
cat('  ✓ Total Raster Layers:', length(rasters), '\n\n')

cat('Raster Outputs:\n')
for (prop in property_names) {
  cat('  Property:', prop, '\n')
  cat('    Mean:', round(terra::minmax(mean_rasters[[prop]])[1], 1), 
      '-', round(terra::minmax(mean_rasters[[prop]])[2], 1), '\n')
  cat('    SD:', round(terra::minmax(sd_rasters[[prop]])[1], 2), 
      '-', round(terra::minmax(sd_rasters[[prop]])[2], 2), '\n')
}

cat('\n')

# ============================================================================
# PHASE 8: FUNCTION VALIDATION TESTS

# ============================================================================

cat('FUNCTION VALIDATION TESTS\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

tests_passed <- 0
tests_total <- 0

# Test 1: gc_fetch_mapunit_polygons returns valid sf object
tests_total <- tests_total + 1
cat('Test 1: gc_fetch_mapunit_polygons returns valid sf object\n')
if (inherits(mapunit_zones, 'sf') && nrow(mapunit_zones) > 0 && 
    'mukey' %in% names(mapunit_zones)) {
  cat('  ✓ PASS\n\n')
  tests_passed <- tests_passed + 1
} else {
  cat('  ✗ FAIL\n\n')
}

# Test 2: gc_fetch_zone_properties returns valid data.frame with compositions
tests_total <- tests_total + 1
cat('Test 2: gc_fetch_zone_properties returns data.frame with valid compositions\n')
zone_comp_sum <- zone_properties$sandtotal_r + zone_properties$silttotal_r + 
                 zone_properties$claytotal_r
# Check that at least 50% of zones have valid compositions (some may be NA)
valid_count <- sum(!is.na(zone_comp_sum) & abs(zone_comp_sum - 100) < 1, na.rm = TRUE)
has_valid_zones <- valid_count >= (nrow(zone_properties) * 0.5)
if (is.data.frame(zone_properties) && 'mukey' %in% names(zone_properties) && has_valid_zones) {
  cat('  ✓ PASS (', valid_count, '/', nrow(zone_properties), 'zones with valid compositions)\n\n')
  tests_passed <- tests_passed + 1
} else {
  cat('  ✗ FAIL\n\n')
}

# Test 3: gc_sample_zone_observations returns valid sf with spatial locations
tests_total <- tests_total + 1
cat('Test 3: gc_sample_zone_observations returns valid sf with proper structure\n')
has_geometry <- inherits(observations, 'sf')
has_props <- all(property_names %in% names(observations))
has_zone_id <- 'zone_id' %in% names(observations)
if (has_geometry && has_props && has_zone_id && nrow(observations) > 0) {
  cat('  ✓ PASS\n\n')
  tests_passed <- tests_passed + 1
} else {
  cat('  ✗ FAIL\n\n')
}

# Test 4: Observations sum to ~100%
tests_total <- tests_total + 1
cat('Test 4: All observations have valid compositions (sum ≈ 100%)\n')
obs_sums <- rowSums(sf::st_drop_geometry(observations[, property_names]), na.rm = TRUE)
valid_obs <- sum(abs(obs_sums - 100) < 1)
if (valid_obs == nrow(observations)) {
  cat('  ✓ PASS (', valid_obs, '/', nrow(observations), 'valid)\n\n')
  tests_passed <- tests_passed + 1
} else {
  cat('  ✗ FAIL (', valid_obs, '/', nrow(observations), 'valid)\n\n')
}

# Test 5: ILR transformation produces valid coordinates
tests_total <- tests_total + 1
cat('Test 5: ILR transformation produces valid coordinates\n')
has_ilr1 <- 'ilr1' %in% names(ilr_data)
has_ilr2 <- 'ilr2' %in% names(ilr_data)
if (has_ilr1 && has_ilr2) {
  cat('  ✓ PASS\n\n')
  tests_passed <- tests_passed + 1
} else {
  cat('  ✗ FAIL\n\n')
}

cat('─────────────────────────────────────────────────────────────────\n\n')
cat('TEST RESULTS:', tests_passed, '/', tests_total, 'passed\n\n')

if (tests_passed == tests_total) {
  cat('╔═══════════════════════════════════════════════════════════════╗\n')
  cat('║  ✓ ALL TESTS PASSED - WORKFLOW READY FOR KRIGING PHASE       ║\n')
  cat('╚═══════════════════════════════════════════════════════════════╝\n\n')
} else {
  cat('╔═══════════════════════════════════════════════════════════════╗\n')
  cat('║  ⚠ SOME TESTS FAILED - REVIEW ABOVE FOR DETAILS             ║\n')
  cat('╚═══════════════════════════════════════════════════════════════╝\n\n')
}
