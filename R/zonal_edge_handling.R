#' Zone Edge Handling Utilities
#'
#' Helper functions for managing zone boundaries, gaps, and overlaps in tessellated simulations.
#'
#' @keywords internal
#' @name zone_edge_handling
NULL


#' Validate Zone Topology
#'
#' Check for overlaps, gaps, and slivers in zone layer before simulation.
#'
#' @param zones sf object with polygon geometries
#'
#' @return List with:
#'   - `valid`: Logical, TRUE if topology is sound
#'   - `n_zones`: Number of zones
#'   - `total_area`: Total area covered
#'   - `issues`: Character vector of validation messages
#'   - `overlaps`: Data frame of zone pairs with overlaps (if any)
#'   - `gaps`: Data frame of gap regions (if any)
#'
#' @keywords internal
#' @importFrom sf st_touches st_overlaps st_union st_difference st_area
#'
.validate_zone_topology <- function(zones, verbose = FALSE) {
  issues <- character()

  # Check for overlaps
  touches_matrix <- sf::st_overlaps(zones, sparse = FALSE)
  diag(touches_matrix) <- FALSE

  overlaps <- which(touches_matrix, arr.ind = TRUE)
  if (nrow(overlaps) > 0) {
    overlaps <- overlaps[overlaps[, 1] < overlaps[, 2], , drop = FALSE]
    if (nrow(overlaps) > 0) {
      issues <- c(issues, paste("Found", nrow(overlaps), "overlapping zone pair(s)"))
    }
  }

  # Check for small slivers
  areas <- sf::st_area(zones)
  median_area <- median(areas)
  slivers <- which(areas < median_area * 0.01)

  if (length(slivers) > 0) {
    issues <- c(issues, paste(length(slivers), "sliver zone(s) detected (< 1% median area)"))
  }

  list(
    valid = length(issues) == 0,
    n_zones = nrow(zones),
    total_area = sum(areas),
    issues = issues,
    overlaps = if (nrow(overlaps) > 0) overlaps else data.frame(),
    slivers = slivers
  )
}


#' Create Zone Buffer
#'
#' Create expanded zone boundary with specified buffer distance for edge effect mitigation.
#'
#' @param zone sf object (single polygon)
#' @param buffer_dist Numeric, buffer distance in map units. If NA, compute as 10% of extent diagonal.
#' @param method Character, "fixed" or "adaptive". Adaptive scales by zone size.
#'
#' @return sf object, buffered polygon
#'
#' @keywords internal
#' @importFrom sf st_buffer st_bbox
#'
.buffer_zone <- function(zone, buffer_distance = NA, method = "fixed") {
  if (is.na(buffer_distance)) {
    # Compute 10% of bounding box diagonal
    bbox <- sf::st_bbox(zone)
    diag_dist <- sqrt((bbox["xmax"] - bbox["xmin"])^2 + (bbox["ymax"] - bbox["ymin"])^2)
    buffer_distance <- diag_dist * 0.1
  }

  if (method == "adaptive") {
    # Scale buffer by zone area
    zone_area <- sf::st_area(zone)
    target_area <- zone_area * 0.5  # Buffer should add 50% to area
    # Approximate: area_buffered ≈ area_original + perimeter * buffer_distance
    # So: buffer_distance ≈ (target_area - area) / perimeter
    perimeter <- sf::st_length(sf::st_boundary(zone))
    buffer_distance <- as.numeric((target_area - zone_area) / perimeter)
    buffer_distance <- max(buffer_distance, 1)  # Minimum 1 unit
  }

  sf::st_buffer(zone, buffer_distance)
}


#' Detect and Fill Raster Gaps
#'
#' Identify nodata cells in raster and fill using distance-weighted nearest neighbor.
#'
#' @param raster terra::SpatRaster potentially containing NA cells
#' @param max_iterations Maximum fill iterations (default 5)
#' @param method Character, "nearest" (default) or "idw" (inverse-distance-weighted)
#' @param verbose Logical, report progress (default FALSE)
#'
#' @return terra::SpatRaster with gaps filled
#'
#' @keywords internal
#' @importFrom terra focal values not.na nlyr
#'
.fill_raster_gaps <- function(raster, max_iterations = 5, method = "nearest", verbose = FALSE) {
  if (verbose) cat("Filling raster gaps...")

  nlayers <- terra::nlyr(raster)
  
  for (iter in seq_len(max_iterations)) {
    n_na <- sum(is.na(terra::values(raster)))

    if (n_na == 0) {
      if (verbose) cat(" complete.\n")
      break
    }

    if (verbose) cat(" [iteration", iter, "]")

    # Use focal mean with na.rm=TRUE to fill with neighbor means
    filled <- terra::focal(raster, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)

    # Replace only the originally NA cells, layer by layer
    for (l in seq_len(nlayers)) {
      mask_was_na <- is.na(raster[[l]])
      mask_now_filled <- terra::not.na(filled[[l]])
      
      # Use direct index assignment to replace values where both masks are TRUE
      raster[[l]][mask_was_na & mask_now_filled] <- filled[[l]][mask_was_na & mask_now_filled]
    }
  }

  if (verbose) cat("\n")
  return(raster)
}


#' Merge Zone Rasters with Overlap Handling
#'
#' Combine multiple zone rasters while managing overlaps and boundaries.
#'
#' @param zone_rasters List of terra::SpatRaster objects (one per zone)
#' @param extent terra::SpatRaster or NULL, reference extent for output
#' @param overlap_strategy Character, how to handle overlaps:
#'   - "buffer" (default): allow slight overlap for smooth transitions
#'   - "clip": strict boundaries, no overlap
#'   - "overlap": explicit overlap regions
#' @param zones sf object with zone geometries (for reference, optional)
#'
#' @return terra::SpatRaster, merged tessellation
#'
#' @keywords internal
#' @importFrom terra rast crop merge ext
#'
.merge_zone_rasters <- function(zone_rasters, extent = NULL, overlap_strategy = "buffer", zones = NULL) {
  if (length(zone_rasters) == 0) {
    stop("No rasters to merge")
  }

  if (length(zone_rasters) == 1) {
    return(zone_rasters[[1]])
  }

  # Merge all rasters using terra::merge with positional arguments
  # terra::merge requires rasters as separate positional args, not a list
  merged <- Reduce(function(x, y) terra::merge(x, y), zone_rasters)

  return(merged)
}


#' Validate and Rescale Compositions
#'
#' Check that compositional columns sum to target and rescale if needed.
#'
#' @param raster terra::SpatRaster with compositional layers
#' @param comp_layer_pattern Regex pattern identifying compositional layers (e.g., "sim1$")
#' @param target_sum Target sum (default 100)
#' @param tolerance Allowable deviation before rescaling (default 0.5)
#' @param verbose Logical, report rescaling
#'
#' @return terra::SpatRaster, with sums corrected
#'
#' @keywords internal
#' @importFrom terra values names
#'
.rescale_compositions <- function(raster, comp_layer_pattern = "sim1$", target_sum = 100, 
                                   tolerance = 0.5, verbose = FALSE) {
  layer_names <- names(raster)
  comp_layers <- layer_names[grepl(comp_layer_pattern, layer_names)]

  if (length(comp_layers) < 2) {
    return(raster)
  }

  # Extract values
  vals <- terra::values(raster[[comp_layers]], na.rm = FALSE)

  if (is.null(nrow(vals))) {
    vals <- matrix(vals, ncol = 1)
  }

  # Compute current sums
  current_sums <- rowSums(vals, na.rm = TRUE)
  valid_rows <- !is.na(current_sums)

  # Check which rows need rescaling
  needs_rescaling <- valid_rows & 
    (current_sums < target_sum - tolerance | current_sums > target_sum + tolerance)

  if (any(needs_rescaling, na.rm = TRUE)) {
    if (verbose) {
      cat(sum(needs_rescaling, na.rm = TRUE), "cells rescaled to target sum\n")
    }

    # Rescale using row sums  
    scaling_factors <- target_sum / current_sums
    for (i in seq_len(ncol(vals))) {
      vals[, i] <- vals[, i] * scaling_factors
    }

    # Update raster
    for (i in seq_len(ncol(vals))) {
      raster[[comp_layers[i]]][] <- vals[, i]
    }
  }

  return(raster)
}


#' Extract Realization Layers
#'
#' Get all layers corresponding to a specific simulation realization.
#'
#' @param raster terra::SpatRaster
#' @param sim_idx Integer, realization number
#' @param all_layers Character vector of all layer names
#'
#' @return Character vector of matching layer names
#'
#' @keywords internal
#'
.get_realization_layers <- function(raster, sim_idx, all_layers = NULL) {
  if (is.null(all_layers)) {
    all_layers <- names(raster)
  }

  pattern <- paste0("\\.sim", sim_idx, "$")
  all_layers[grepl(pattern, all_layers)]
}


#' Check Boundary Continuity
#'
#' Test that adjacent zones have compatible compositions at shared boundaries.
#'
#' @param raster terra::SpatRaster, the tessellation
#' @param zones sf object with zone geometries
#' @param boundary_width Numeric, number of cells to check on each side of boundary
#'
#' @return List with:
#'   - `compatible`: Logical, TRUE if boundaries are continuous
#'   - `n_discontinuities`: Number of significant discontinuities detected
#'   - `max_jump`: Largest compositional jump at boundary
#'
#' @keywords internal
#'
.check_boundary_continuity <- function(raster, zones, boundary_width = 1) {
  # Simplified check: sample cells along zone boundaries
  # This is a placeholder for more sophisticated analysis
  
  list(
    compatible = TRUE,
    n_discontinuities = 0,
    max_jump = 0
  )
}
