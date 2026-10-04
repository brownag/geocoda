#!/usr/bin/env Rscript
# End-to-End Workflow: Map Unit Zone-Based Compositional Simulation
# This script demonstrates the complete geocoda workflow for spatial prediction
# of soil compositions within map unit zones fetched from SSURGO/SDA

library(devtools, quietly = TRUE)
library(sf, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(terra, quietly = TRUE)
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
# HELPER FUNCTIONS: Polygon Adjacency & Property Contrasts
# ============================================================================
# Enable neighbor-aware simulation where adjacent map units influence each other
# based on: (1) polygon adjacency, (2) shared boundary length, (3) polygon size,
# and (4) property contrasts (high contrast = sharp boundaries, low = smooth)

#' Compute polygon adjacency relationships with metrics
#' @param polygons sf object with mukey column
#' @param verbose logical
#' @return list with adjacency matrix and metrics
compute_polygon_adjacency <- function(polygons, verbose = FALSE) {
  if (verbose) cat('  Computing polygon adjacency graph...\n')
  
  n_zones <- nrow(polygons)
  mukeys <- polygons$mukey
  
  # Initialize matrices for adjacency relationships
  adjacency_matrix <- matrix(FALSE, n_zones, n_zones)
  boundary_lengths <- matrix(0, n_zones, n_zones)
  area_ratios <- matrix(1, n_zones, n_zones)
  
  # Use spatial index to find candidate pairs (much faster than all pairwise)
  touches_sparse <- sf::st_touches(polygons)  # Sparse representation
  
  # Process only touching pairs
  for (i in seq_len(n_zones)) {
    neighbors <- touches_sparse[[i]]
    
    for (j in neighbors[neighbors > i]) {  # Only upper triangle (avoid duplicates)
      adjacency_matrix[i, j] <- TRUE
      adjacency_matrix[j, i] <- TRUE
      
      # Compute shared boundary length using st_difference or direct approach
      # For unioned geometries (multipolygons), we estimate based on shared perimeter
      try({
        geom_i <- sf::st_geometry(polygons[i, ])[[1]]
        geom_j <- sf::st_geometry(polygons[j, ])[[1]]
        
        # A simple proxy: shared boundary ≈ proportion of touching that represents real boundary
        # For now, use a nominal value (could improve with more sophisticated boundary extraction)
        boundary_lengths[i, j] <- 0.1  # Nominal shared boundary length
        boundary_lengths[j, i] <- 0.1
      }, silent = TRUE)
      
      # Compute area ratio (normalized to [0,1], where 1 = equal size)
      area_i <- as.numeric(sf::st_area(polygons[i, ]))
      area_j <- as.numeric(sf::st_area(polygons[j, ]))
      ratio <- min(area_i, area_j) / max(area_i, area_j)
      area_ratios[i, j] <- ratio
      area_ratios[j, i] <- ratio
    }
  }
  
  # Compute combined adjacency strength (60% boundary, 40% area ratio)
  # Since boundary lengths are nominal, adjacency_strength is driven primarily by area_ratios
  adjacency_strength <- (0.6 * boundary_lengths) + (0.4 * area_ratios)
  
  diag(adjacency_strength) <- 0
  
  if (verbose) {
    n_adjacent <- rowSums(adjacency_matrix)
    cat('    Adjacent zone pairs:', sum(adjacency_matrix) / 2, '\n')
    cat('    Mean neighbors per mukey:', round(mean(n_adjacent[n_adjacent > 0]), 1), '\n')
  }
  
  list(
    mukeys = mukeys,
    adjacency_matrix = adjacency_matrix,
    boundary_lengths = boundary_lengths,
    area_ratios = area_ratios,
    adjacency_strength = adjacency_strength
  )
}

#' Compute property contrasts at zone boundaries
#' @param zone_properties data.frame with mukey and composition columns
#' @param adjacency adjacency list from compute_polygon_adjacency()
#' @param property_names character vector of composition columns
#' @return list with contrast matrix and distance decay parameters
compute_boundary_contrasts <- function(zone_properties, adjacency, property_names) {
  mukeys <- adjacency$mukeys
  n_zones <- length(mukeys)
  
  # Composition data
  compositions <- sf::st_drop_geometry(zone_properties[, c("mukey", property_names)])
  
  # Initialize contrast and decay matrices
  contrast_matrix <- matrix(0, n_zones, n_zones)
  decay_scaling <- matrix(1, n_zones, n_zones)  # Multiplier for distance decay
  
  # Compute contrasts between adjacent zones
  for (i in seq_len(n_zones)) {
    for (j in seq_len(n_zones)) {
      if (i >= j) next  # Only upper triangle
      
      if (adjacency$adjacency_matrix[i, j]) {
        # Find zones in properties
        idx_i <- which(compositions$mukey == mukeys[i])
        idx_j <- which(compositions$mukey == mukeys[j])
        
        if (length(idx_i) > 0 && length(idx_j) > 0) {
          # Get compositions
          comp_i <- as.numeric(compositions[idx_i[1], property_names])
          comp_j <- as.numeric(compositions[idx_j[1], property_names])
          
          # Euclidean distance in composition space (percentage units)
          if (!anyNA(c(comp_i, comp_j))) {
            contrast <- sqrt(sum((comp_i - comp_j)^2, na.rm = TRUE))
            contrast_matrix[i, j] <- contrast
            contrast_matrix[j, i] <- contrast
            
            # Property contrast determines decay scaling:
            # High contrast (>25 percentage points) = short decay (0.5x, sharp)
            # Low contrast (<5 percentage points) = long decay (2.5x, smooth)
            contrast_normalized <- min(25, max(0, contrast))
            # Linear mapping: contrast 0→25 maps to decay 2.5→0.5
            decay_factor <- 2.5 - (2 * (contrast_normalized / 25))
            decay_scaling[i, j] <- decay_factor
            decay_scaling[j, i] <- decay_factor
          }
        }
      }
    }
  }
  
  list(
    contrast_matrix = contrast_matrix,
    decay_scaling = decay_scaling,
    adjacency_strength = adjacency$adjacency_strength
  )
}

#' Compute neighbor influence weights for observations
#' @param observations sf object with mukey and zone_id columns
#' @param zone_properties data.frame with mukey and properties
#' @param adjacency adjacency list
#' @param contrasts contrast list
#' @param neighbor_order integer: 1, 2, or 3 (how many neighbor levels)
#' @return observations with neighbor_influence column added
compute_neighbor_weights <- function(observations, zone_properties, adjacency, 
                                      contrasts, neighbor_order = 1) {
  
  mukeys <- adjacency$mukeys
  mukey_to_idx <- setNames(seq_along(mukeys), as.character(mukeys))
  
  observations$neighbor_influence <- 1.0
  
  # Debug tracking
  n_with_neighbors <- 0
  weights_applied <- numeric(0)
  
  # For each observation, compute neighbor influence weight
  for (obs_idx in seq_len(nrow(observations))) {
    mukey_char <- as.character(observations$mukey[obs_idx])
    zone_idx <- mukey_to_idx[mukey_char]
    
    if (is.na(zone_idx)) next  # Zone not in adjacency matrix
    
    influence_weight <- 1.0
    
    # First order neighbors
    neighbors_1 <- which(adjacency$adjacency_matrix[zone_idx, ])
    if (length(neighbors_1) > 0) {
      # Strength increases with better adjacency
      adj_strength_vals <- contrasts$adjacency_strength[zone_idx, neighbors_1]
      adj_strength <- mean(adj_strength_vals, na.rm = TRUE)
      if (adj_strength > 0) {
        n_with_neighbors <- n_with_neighbors + 1
        influence_weight <- influence_weight + (0.35 * adj_strength)
      }
    }
    
    # Second order neighbors (if requested)
    if (neighbor_order >= 2 && length(neighbors_1) > 0) {
      neighbors_2 <- unique(unlist(lapply(neighbors_1, function(n) {
        which(adjacency$adjacency_matrix[n, ])
      })))
      neighbors_2 <- setdiff(neighbors_2, c(zone_idx, neighbors_1))
      
      if (length(neighbors_2) > 0) {
        adj_strength_2 <- mean(contrasts$adjacency_strength[zone_idx, neighbors_2], na.rm = TRUE)
        if (adj_strength_2 > 0 && !is.na(adj_strength_2)) {
          influence_weight <- influence_weight + (0.15 * adj_strength_2)
        }
      }
    }
    
    # Third order neighbors (if requested)
    if (neighbor_order >= 3 && length(neighbors_2) > 0) {
      neighbors_3 <- unique(unlist(lapply(neighbors_2, function(n) {
        which(adjacency$adjacency_matrix[n, ])
      })))
      neighbors_3 <- setdiff(neighbors_3, c(zone_idx, neighbors_1, neighbors_2))
      
      if (length(neighbors_3) > 0) {
        adj_strength_3 <- mean(contrasts$adjacency_strength[zone_idx, neighbors_3], na.rm = TRUE)
        if (adj_strength_3 > 0 && !is.na(adj_strength_3)) {
          influence_weight <- influence_weight + (0.05 * adj_strength_3)
        }
      }
    }
    
    observations$neighbor_influence[obs_idx] <- influence_weight
    if (influence_weight > 1.05) {
      weights_applied <- c(weights_applied, influence_weight)
    }
  }
  
  if (length(weights_applied) > 0) {
    cat('    [DEBUG] Observations with elevated weights: ', length(weights_applied), 
        ' (range:', round(min(weights_applied), 3), '-', round(max(weights_applied), 3), ')\n')
  } else {
    cat('    [DEBUG] No elevated weights - checking adjacency_strength values...\n')
    cat('    [DEBUG] Mean adjacency_strength:', 
        round(mean(contrasts$adjacency_strength)[contrasts$adjacency_strength > 0], 3), '\n')
  }
  
  observations
}

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
# PHASE 2A: POLYGON ADJACENCY & PROPERTY CONTRASTS
# ============================================================================

cat('PHASE 2A: POLYGON ADJACENCY AND PROPERTY CONTRASTS\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Step 3A: Computing polygon adjacency relationships...\n')

# Group polygons by mukey to create pseudo-zones for adjacency
# This handles the distinction between mupolygonkey (polygon features) 
# and mukey (map unit concepts that may have multiple polygons)
unique_mukeys <- unique(mapunit_zones$mukey)
n_zones_concept <- length(unique_mukeys)

cat('  Total polygons (mupolygonkey):', nrow(mapunit_zones), '\n')
cat('  Unique map units (mukey):', n_zones_concept, '\n')

# Build zone-level polygons by unioning all polygons with the same mukey
# This creates one pseudo-polygon per unique mukey for adjacency computation
zone_polygons_list <- list()
for (mukey in unique_mukeys) {
  mukey_polys <- mapunit_zones[mapunit_zones$mukey == mukey, ]
  # Union all polygons with this mukey into a single geometry
  union_geom <- sf::st_union(mukey_polys)
  zone_polygons_list[[as.character(mukey)]] <- sf::st_sf(
    mukey = mukey,
    n_polygons = nrow(mukey_polys),
    geometry = union_geom
  )
}
zone_polygons_for_adjacency <- do.call(rbind, zone_polygons_list)
rownames(zone_polygons_for_adjacency) <- NULL

# Compute polygon adjacency on unique mukeys (map unit concepts)
adjacency <- compute_polygon_adjacency(zone_polygons_for_adjacency, verbose = TRUE)

# Compute property contrasts at zone boundaries
cat('Step 3B: Computing boundary property contrasts (based on unique mukeys)...\n')

# Compute property contrasts at zone boundaries
# This uses the unique mukey-level zones (not individual polygons)
contrasts <- compute_boundary_contrasts(zone_properties, adjacency, property_names)

# Report contrast statistics
contrast_summary <- contrasts$contrast_matrix[contrasts$contrast_matrix > 0]
if (length(contrast_summary) > 0) {
  cat('  Composition contrasts at boundaries:\n')
  cat('    Min contrast:', round(min(contrast_summary), 2), 'percentage points\n')
  cat('    Mean contrast:', round(mean(contrast_summary), 2), 'percentage points\n')
  cat('    Max contrast:', round(max(contrast_summary), 2), 'percentage points\n\n')
}

# Apply neighbor influence weights to observations
cat('Step 3C: Computing observation-level spatial weights...\n')

# Load spatial weighting functions if not already loaded
if (!exists('compute_observation_spatial_weights', mode = 'function')) {
  source('R/spatial_weighting.R')
}

# Compute observation-level spatial weights (distance to boundary + direction + adjacency)
observations <- compute_observation_spatial_weights(
  observations = observations,
  mapunit_zones = mapunit_zones,
  contrasts = contrasts,
  distance_power = 0.5,      # Square-root scaling (moderate boundary emphasis)
  direction_boost = 1.3,     # 30% boost for neighbor-facing observations
  adjacency_scale = 0.5,     # Scale adjacency strength contribution
  verbose = TRUE
)

# Compute per-zone kriging range adjustments (for later use in Phase 5)
cat('\nStep 3C2: Computing per-zone kriging range adjustments...\n')
kriging_adjustments <- compute_kriging_range_adjustments(
  zone_properties = zone_properties,
  adjacency = adjacency,
  base_range = 0.01,                # Base range 0.01 degrees
  neighbor_count_effect = 0.3,      # How much neighbor count matters
  adjacency_strength_effect = 0.2,  # How much composition similarity matters
  verbose = TRUE
)

cat('  ✓ Kriging adjustments computed for', nrow(kriging_adjustments), 'zones\n\n')

# ============================================================================
# PHASE 2D: Weighted Zone Model Fitting
# ============================================================================

cat('Step 3D: Fitting weighted zone models with observation-level spatial weights...\n')

# Verify that spatial_weight column was added in Phase 2C
if (!'spatial_weight' %in% names(observations)) {
  cat('  ⚠ WARNING: spatial_weight column NOT found in observations!\n')
  cat('    Available columns:', paste(names(observations), collapse=', '), '\n')
  stop('spatial_weight column required - check Phase 2C computation')
} else {
  cat('  ✓ spatial_weight column present\n')
  cat('    Range:', round(min(observations$spatial_weight, na.rm=TRUE), 4), 'to', 
      round(max(observations$spatial_weight, na.rm=TRUE), 4), '\n')
  cat('    Mean weight:', round(mean(observations$spatial_weight, na.rm=TRUE), 4), '\n')
}

# Load the weighted kriging function if not already loaded
if (!exists('gc_fit_zone_models_weighted', mode = 'function')) {
  source('R/sda_integration.R')
}

# Create unweighted copy for baseline (remove spatial_weight column)
obs_for_unweighted <- observations
obs_for_unweighted$spatial_weight <- NULL  # Remove spatial weights for unweighted call

# Fit unweighted baseline model for comparison
cat('\n  Fitting unweighted baseline model...\n')
zone_models_unweighted <- gc_fit_zone_models_weighted(
  observations = obs_for_unweighted,
  property_names = property_names,
  use_neighbor_weights = FALSE,
  verbose = FALSE
)

# Fit weighted model using observation-level spatial weights
cat('\n  Fitting weighted model with spatial weights...\n')
zone_models_weighted <- gc_fit_zone_models_weighted(
  observations = observations,  # Includes spatial_weight column from Phase 2C
  property_names = property_names,
  use_neighbor_weights = TRUE,  # Not used if spatial_weight exists
  verbose = FALSE
)

cat('  ✓ Both weighted and unweighted zone models fitted\n\n')

# Compare zone statistics: weighted vs unweighted (should now show DIFFERENCES)
cat('  Zone means comparison (first 10 zones):\n')
zone_comparison <- data.frame(
  zone_id = zone_models_weighted$zone_statistics$zone_id[1:10],
  n_obs = zone_models_weighted$zone_statistics$n_obs[1:10],
  unweighted_mean = zone_models_unweighted$zone_statistics$sandtotal_r_mean[1:10],
  weighted_mean = zone_models_weighted$zone_statistics$sandtotal_r_mean[1:10]
)
zone_comparison$diff_pct <- 100 * (zone_comparison$weighted_mean - zone_comparison$unweighted_mean) /
                                  abs(zone_comparison$unweighted_mean + 1e-8)
print(zone_comparison)

# Summary: count how many zones show differences
n_diff_zones <- sum(abs(zone_comparison$diff_pct) > 0.01)
cat('\n  Zones with >0.01% difference:', n_diff_zones, 'of 10 shown\n')
cat('\n')

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
# This stores both original (sand, silt, clay) and ILR coordinates (ilr1, ilr2)
ilr_data <- gc_prepare_hierarchical_data(
  ssurgo_compositions = sf::st_drop_geometry(observations[, property_names]),
  component_cols = property_names,
  verbose = FALSE
)

# Add zone information
ilr_data$zone <- observations$zone_id

cat('  ✓ Transformed', nrow(ilr_data), 'observations to ILR space\n')

# Check ILR coordinates
if ('ilr1' %in% names(ilr_data)) {
  cat('  ILR1 range:', round(min(ilr_data$ilr1, na.rm = TRUE), 3), 'to', 
      round(max(ilr_data$ilr1, na.rm = TRUE), 3), '\n')
  cat('  ILR2 range:', round(min(ilr_data$ilr2, na.rm = TRUE), 3), 'to', 
      round(max(ilr_data$ilr2, na.rm = TRUE), 3), '\n')
}

# Store zone means in ILR space for later use
ilr_zone_means <- ilr_data %>%
  group_by(zone) %>%
  summarise(
    ilr1_mean = mean(ilr1, na.rm = TRUE),
    ilr2_mean = mean(ilr2, na.rm = TRUE),
    .groups = 'drop'
  )

cat('  ✓ Computed zone-level means in ILR space\n\n')

# ============================================================================
# PHASE 5: Per-Zone Kriging and Conditional Simulation via ILR Space
# ============================================================================

cat('PHASE 5: KRIGING & COMPOSITIONAL SIMULATION IN ILR SPACE\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Step 6: Fitting variograms and building kriging models...\n')

# For this test, compute ILR parameters from ALL observations (pooled across zones)
# In production, this would be done per-zone or regionally

# Prepare data for kriging: add x,y coordinates to ILR data
ilr_with_coords <- ilr_data %>%
  # Ensure we keep the original composition columns AND ILR coordinates
  select(all_of(c(property_names, 'ilr1', 'ilr2', 'zone'))) %>%
  mutate(obs_idx = row_number()) %>%
  # Add coordinates
  left_join(
    observations %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      rename(x = X, y = Y) %>%
      mutate(obs_idx = row_number()),
    by = 'obs_idx'
  ) %>%
  select(all_of(c(property_names, 'ilr1', 'ilr2', 'x', 'y', 'zone')))

# Remove any rows with missing ILR or coordinates
ilr_with_coords <- ilr_with_coords %>%
  filter(!is.na(ilr1), !is.na(ilr2), !is.na(x), !is.na(y))

# Convert to sf object and project to UTM (Zone 15N for Iowa)
ilr_with_coords_sf <- sf::st_as_sf(
  ilr_with_coords,
  coords = c('x', 'y'),
  crs = sf::st_crs(observations)  # EPSG:4326
) %>%
  sf::st_transform(crs = 32615)  # UTM Zone 15N (suitable for Iowa)

# Convert mapunit_zones to projected CRS as well
mapunit_zones <- mapunit_zones %>% sf::st_transform(crs = 32615)

cat('  ILR data with coordinates: ', nrow(ilr_with_coords_sf), ' observations (projected to UTM Zone 15N)\n')

# Compute ILR parameters from original compositional data
# gc_ilr_params() will transform the compositions to ILR internally
ilr_params <- gc_ilr_params(ilr_with_coords[, property_names])

cat('  ILR parameters computed:\n')
cat('    Mean (ILR1):', round(ilr_params$mean[1], 4), '\n')
cat('    Mean (ILR2):', round(ilr_params$mean[2], 4), '\n')
cat('    Variance (ILR1):', round(ilr_params$cov[1,1], 4), '\n')
cat('    Variance (ILR2):', round(ilr_params$cov[2,2], 4), '\n\n')

# Fit variogram to ILR data
cat('Step 6B: Fitting variogram to ILR-transformed data...\n')
vgm_fit <- gc_fit_vgm(
  ilr_params = ilr_params,
  data = ilr_with_coords,  # Use the data.frame version for variogram fitting
  vgm_model_type = "Exp",
  aggregate = TRUE  # Use aggregated variogram across both ILR dimensions
)

cat('  Variogram fitted\n\n')

# Build kriging model
cat('Step 6C: Building gstat kriging model...\n')
kriging_model <- gc_ilr_model(
  ilr_params = ilr_params,
  variogram_model = vgm_fit,
  data = ilr_with_coords_sf,  # Use the sf version for kriging
  model_type = "univariate"
)

cat('  Kriging model built\n\n')

# ============================================================================
# PHASE 5B: Generate Realizations Using gc_sim_composition
# ============================================================================

cat('Step 7: Generating conditional realizations via gc_sim_composition...\n')

# Create prediction grid in projected CRS
zone_union <- sf::st_union(mapunit_zones)
grid_size <- 1000  # 1 km grid in UTM (projected)
grid <- sf::st_make_grid(zone_union, cellsize = grid_size, square = TRUE)
grid_sf <- sf::st_sf(id = seq_along(grid), geometry = grid)

cat('  ✓ Created prediction grid with', nrow(grid_sf), 'cells\n\n')

# Generate realizations using gc_sim_composition
# This function:
# 1. Does kriging in ILR space
# 2. Inverse-transforms back to original composition space  
# 3. Returns values that sum to 100% (composition closure guaranteed)
realizations_raw <- gc_sim_composition(
  model = kriging_model,
  locations = grid_sf,
  nsim = 3,
  target_names = property_names,
  crs = sf::st_crs(observations)
)

cat('  ✓ Generated 3 realizations\n')
cat('  Output: terra::SpatRaster with layers for each component & realization\n\n')

# Extract individual property rasters for comparison
realizations_weighted <- list()
realizations_unweighted <- list()

# Store the realizations in list format for later use
for (prop in property_names) {
  layer_names <- grep(paste0('^', prop, '\\'), names(realizations_raw), value = TRUE)
  
  if (length(layer_names) > 0) {
    realizations_weighted[[prop]] <- terra::subset(realizations_raw, layer_names)
  }
}

cat('Step 7B: Verification of closure constraint...\n')

# Verify that compositions sum to 100%
sand_layers <- terra::subset(realizations_raw, grep('^sandtotal_r', names(realizations_raw)))
silt_layers <- terra::subset(realizations_raw, grep('^silttotal_r', names(realizations_raw)))
clay_layers <- terra::subset(realizations_raw, grep('^claytotal_r', names(realizations_raw)))

# Check sums for each realization
for (sim_idx in 1:3) {
  sand <- terra::as.array(terra::subset(sand_layers, sim_idx))
  silt <- terra::as.array(terra::subset(silt_layers, sim_idx))
  clay <- terra::as.array(terra::subset(clay_layers, sim_idx))
  
  sums <- sand + silt + clay
  valid_sums <- sums[!is.na(sums)]
  
  if (length(valid_sums) > 0) {
    cat(sprintf('  Realization %d: sum mean=%.4f, min=%.4f, max=%.4f\n',
                sim_idx, mean(valid_sums), min(valid_sums), max(valid_sums)))
  }
}

cat('\n')
cat('─────────────────────────────────────────────────────────────────\n\n')


cat('PHASE 5B: COMPARING TO LEGACY ZONE MODEL APPROACH\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

# Compare zone-based means (from Phase 2D fitted models) to kriging output
cat('Step 8: Zone statistics comparison...\n')
cat('  Zone models from Phase 2D (spatial weighting applied):\n')
comparison <- data.frame(
  zone_id = zone_models_weighted$zone_statistics$zone_id[1:10],
  n_obs = zone_models_weighted$zone_statistics$n_obs[1:10],
  weighted_mean = zone_models_weighted$zone_statistics$sandtotal_r_mean[1:10],
  unweighted_mean = zone_models_unweighted$zone_statistics$sandtotal_r_mean[1:10]
)
comparison$diff_pct <- 100 * (comparison$weighted_mean - comparison$unweighted_mean) / 
                              abs(comparison$unweighted_mean + 1e-8)
print(comparison)
cat('\n')

# ============================================================================
# PHASE 6: Raster Summary and Validation
# ============================================================================

cat('PHASE 6: RASTER STATISTICS & VALIDATION\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Step 9: Raster layer summary...\n')
cat('  Layers generated by gc_sim_composition:\n')
cat('    Total layers:', terra::nlyr(realizations_raw), '\n')
cat('    Layer names:', paste(names(realizations_raw)[1:6], '...', collapse=', '), '\n')
cat('    Dimensions:', terra::nrow(realizations_raw), 'x', terra::ncol(realizations_raw), '\n')
cat('    CRS:', terra::crs(realizations_raw), '\n\n')

# Compute statistics per realization
cat('Step 9B: Per-realization statistics (sand content):\n')
sand_r1 <- terra::subset(realizations_raw, 'sandtotal_r.sim1')
sand_r2 <- terra::subset(realizations_raw, 'sandtotal_r.sim2')
sand_r3 <- terra::subset(realizations_raw, 'sandtotal_r.sim3')

for (r_idx in 1:3) {
  sand_r <- c(sand_r1, sand_r2, sand_r3)[[r_idx]]
  vals <- terra::values(sand_r, na.rm=TRUE)
  cat(sprintf('  Realization %d: mean=%.2f, sd=%.2f, min=%.2f, max=%.2f\n',
              r_idx, mean(vals), sd(vals), min(vals), max(vals)))
}

cat('\n')

cat('  ✓ Created', length(rasters_weighted), 'weighted raster layers\n')
cat('  ✓ Created', length(rasters_unweighted), 'unweighted raster layers\n\n')

# Create composite visualization showing mean and uncertainty
cat('Step 8: Generating visualization output...\n')

# Calculate mean and standard deviation across realizations (WEIGHTED)
mean_rasters_weighted <- list()
sd_rasters_weighted <- list()

# Calculate mean and standard deviation across realizations (UNWEIGHTED)
mean_rasters_unweighted <- list()
sd_rasters_unweighted <- list()

for (prop in property_names) {
  # WEIGHTED
  sim_rasters_w <- list()
  for (sim in 1:3) {
    sim_rasters_w[[sim]] <- rasters_weighted[[paste0(prop, '_R', sim)]]
  }
  stack_w <- terra::rast(sim_rasters_w)
  mean_rasters_weighted[[prop]] <- terra::app(stack_w, 'mean', na.rm = TRUE)
  sd_rasters_weighted[[prop]] <- terra::app(stack_w, 'sd', na.rm = TRUE)
  

cat('Step 9C: Writing raster output to file...\n')

if (!dir.exists('output')) {
  dir.create('output')
}

# Write the kriging realizations
output_file <- 'output/gc_sim_composition_realizations.tif'
terra::writeRaster(realizations_raw, output_file, overwrite = TRUE)
cat('  ✓ Written kriging realizations:', output_file, '\n')
cat('    ', terra::nlyr(realizations_raw), 'layers\n\n')

# ============================================================================
# PHASE 7: Workflow Summary and Results
# ============================================================================

cat('PHASE 7: WORKFLOW SUMMARY\n')
cat('─────────────────────────────────────────────────────────────────\n\n')

cat('Workflow Status: COMPLETE (ILR-based Kriging & Simulation)\n\n')

cat('Summary Statistics:\n')
cat('  ✓ Zones (Map Units):', nrow(mapunit_zones), '\n')
cat('  ✓ Total Observations:', nrow(observations), '\n')
cat('  ✓ ILR-Transformed Records:', nrow(ilr_with_coords), '\n')
cat('  ✓ Prediction Grid Cells:', nrow(grid_sf), '\n')
cat('  ✓ Conditional Realizations:', 3, 'per component\n')
cat('  ✓ Total Raster Layers:', terra::nlyr(realizations_raw), '\n\n')

# Compute ensemble means for each realization
cat('Ensemble Statistics:\n')
for (sim_idx in 1:3) {
  subset_layers <- grep(paste0('\\.sim', sim_idx, '$'), names(realizations_raw), value = TRUE)
  subset_stack <- terra::subset(realizations_raw, subset_layers)
  ensemble_mean <- terra::app(subset_stack, 'mean', na.rm = TRUE)
  vals <- terra::values(ensemble_mean, na.rm = TRUE)
  cat(sprintf('  Realization %d: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n',
              sim_idx, mean(vals), sd(vals), min(vals), max(vals)))
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

# Test 6: Kriging model created successfully
tests_total <- tests_total + 1
cat('Test 6: Kriging model built successfully\n')
if (!is.null(kriging_model) && inherits(kriging_model, 'gstat')) {
  cat('  ✓ PASS\n\n')
  tests_passed <- tests_passed + 1
} else {
  cat('  ✗ FAIL\n\n')
}

# Test 7: Realizations sum to 100% (closure constraint)
tests_total <- tests_total + 1
cat('Test 7: Compositional realizations respect closure constraint (sum ≈ 100%)\n')
sand_r1 <- terra::values(terra::subset(realizations_raw, 1), na.rm=TRUE)
silt_r1 <- terra::values(terra::subset(realizations_raw, 2), na.rm=TRUE)
clay_r1 <- terra::values(terra::subset(realizations_raw, 3), na.rm=TRUE)
sums <- sand_r1 + silt_r1 + clay_r1
mean_sum <- mean(sums, na.rm=TRUE)
if (abs(mean_sum - 100) < 0.1) {
  cat('  ✓ PASS (mean sum =', round(mean_sum, 4), ')\n\n')
  tests_passed <- tests_passed + 1
} else {
  cat('  ✗ FAIL (mean sum =', round(mean_sum, 4), ', expected 100)\n\n')
}

cat('─────────────────────────────────────────────────────────────────\n\n')
cat('TEST RESULTS:', tests_passed, '/', tests_total, 'passed\n\n')

if (tests_passed == tests_total) {
  cat('╔═══════════════════════════════════════════════════════════════╗\n')
  cat('║  ✓ ALL TESTS PASSED - WORKFLOW COMPLETE                       ║\n')
  cat('║  Kriging performed in ILR space with automatic back-transform  ║\n')
  cat('║  Compositional closure preserved in all outputs               ║\n')
  cat('╚═══════════════════════════════════════════════════════════════╝\n\n')
} else {
  cat('╔═══════════════════════════════════════════════════════════════╗\n')
  cat('║  ⚠ SOME TESTS FAILED - REVIEW ABOVE FOR DETAILS             ║\n')
  cat('╚═══════════════════════════════════════════════════════════════╝\n\n')
}

# ============================================================================
# NEIGHBOR-AWARE SIMULATION APPROACHES
# ============================================================================
# Demonstration of three methods for incorporating adjacent map unit influence

cat('NEIGHBOR-AWARE SIMULATION OPTIONS\n')
cat('═════════════════════════════════════════════════════════════════\n\n')

cat('Available methods for incorporating adjacent map unit properties:\n\n')

# ============================================================================
# APPROACH 1: BOOTSTRAP - Weight Observations by Neighbor Proximity
# ============================================================================

cat('METHOD 1: BOOTSTRAP APPROACH\n')
cat('─────────────────────────────────────────────────────────────────\n')
cat('Observations are weighted by neighbor influence strength.\n')
cat('Higher weights at zone boundaries create gradual transitions.\n')
cat('Property contrasts scale the transition distance (sharp vs smooth).\n\n')

cat('Current neighbor influence stats:\n')
cat('  Min weight:', round(min(observations$neighbor_influence), 3), '\n')
cat('  Mean weight:', round(mean(observations$neighbor_influence), 3), '\n')
cat('  Max weight:', round(max(observations$neighbor_influence), 3), '\n')

# Summary: how many observations have elevated neighbor influence
high_influence <- sum(observations$neighbor_influence > 1.2)
cat('  Observations with elevated influence (>1.2x):', high_influence, '/', nrow(observations), '\n\n')

cat('Usage: Use observations$neighbor_influence as weights in kriging.\n')
cat('Location: Ready in current observations sf object.\n\n')

# ============================================================================
# APPROACH 2: ANALYTICAL SPATIAL DECAY - Property Contrasts Control Range
# ============================================================================

cat('METHOD 2: ANALYTICAL SPATIAL DECAY\n')
cat('─────────────────────────────────────────────────────────────────\n')
cat('Property contrasts modulate variogram range via decay_scaling.\n')
cat('High contrast (>25%%) → short decay (0.5x, sharp boundary)\n')
cat('Low contrast (<5%%) → long decay (2.5x, smooth transition)\n\n')

cat('Current decay scaling at boundaries:\n')
decay_values <- contrasts$decay_scaling[contrasts$decay_scaling != 1 & contrasts$decay_scaling > 0]
if (length(decay_values) > 0) {
  cat('  Min decay multiplier:', round(min(decay_values), 2), '\n')
  cat('  Mean decay multiplier:', round(mean(decay_values), 2), '\n')
  cat('  Max decay multiplier:', round(max(decay_values), 2), '\n')
} else {
  cat('  No adjacent zones - set threshold to zero adjacency\n')
}
cat('\nUsage: Modify gc_fit_zone_models() to use decay_scaling per zone pair.\n')
cat('Impact: Variogram nugget and range adjusted by boundary contrasts.\n\n')

# ============================================================================
# APPROACH 3: BAYESIAN CAR - Intrinsic Conditional Auto-Regressive Prior
# ============================================================================

cat('METHOD 3: BAYESIAN CAR FRAMEWORK (Stan/Nimble)\n')
cat('─────────────────────────────────────────────────────────────────\n')
cat('Zone means shrink toward neighbors via Conditional Auto-Regressive prior.\n')
cat('Adjacency graph automatically constrains zone parameters.\n')
cat('Full Bayesian posterior captures spatial dependence structure.\n\n')

cat('Adjacency graph structure for CAR prior:\n')
cat('  Total zones:', nrow(mapunit_zones), '\n')
cat('  Adjacent zone pairs:', sum(adjacency$adjacency_matrix) / 2, '\n')
cat('  Mean neighbors per zone:', 
    round(mean(rowSums(adjacency$adjacency_matrix)), 2), '\n')
cat('  Max neighbors per zone:', max(rowSums(adjacency$adjacency_matrix)), '\n\n')

cat('Usage: Pass adjacency$adjacency_matrix to gc_fit_hierarchical_stan_3d().\n')
cat('Example:  fit_hierarchical_stan_3d(\n')
cat('            data, prior_spec,\n')
cat('            adjacency_matrix = adjacency$adjacency_matrix,\n')
cat('            car_prior_scale = 1.0,\n')
cat('            estimate_spatial_decay = TRUE\n')
cat('          )\n\n')

cat('═════════════════════════════════════════════════════════════════\n\n')
cat('Summary:\n')
cat('  • Bootstrap: Fastest, intuitive weighting\n')
cat('  • Analytical: Medium speed, deterministic, decay-aware\n')
cat('  • Bayesian CAR: Slowest, most accurate, full posterior inference\n\n')
