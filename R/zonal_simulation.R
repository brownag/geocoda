#' Prepare Zone-Specific Compositional Data
#'
#' Extract compositional samples for each zone from a vector or raster zone layer.
#' Validates sample counts, spatial coverage, and compositional validity.
#'
#' @param data An sf object or data frame with sample coordinates, compositional columns,
#'   and optional zone assignments. Must contain columns `x` and `y` for coordinates.
#' @param zones An sf object (polygons) with a `zone_id` column (unique zone identifiers),
#'   or a `terra::SpatRaster` with cell values representing zone IDs (0 = nodata).
#' @param zone_id_column Character. Name of column in `data` containing zone assignments.
#'   If NULL (default), performs spatial join between data and zones.
#' @param min_samples Numeric. Minimum samples required per zone (default 3).
#'   Zones with fewer samples are flagged in output but not excluded.
#'
#' @return A list with one element per zone, containing:
#'   - `zone_id`: Unique zone identifier
#'   - `n_samples`: Number of samples in this zone
#'   - `data`: sf/data.frame with samples for this zone
#'   - `valid`: Logical, TRUE if n_samples >= min_samples
#'   - `bounds`: Bounding box (sf object or terra extent)
#'   - `warning_msgs`: Character vector of validation warnings (empty if none)
#'
#' Class attribute: `"gc_zone_data"` for S3 method dispatch.
#'
#' @details
#' **Workflow:**
#' 1. If `zone_id_column` is NULL, performs spatial join via `sf::st_join()`
#' 2. Groups samples by zone_id
#' 3. For each zone, validates:
#'    - Sample count >= min_samples
#'    - No missing values in compositional columns
#'    - Compositional sums between 95-105 (assuming 100-sum compositions)
#' 4. Computes zone extent (bounding box)
#'
#' **Warnings Issued:**
#' - Zone has < 5 samples: variogram fitting may be unreliable
#' - Zone has < 10 samples: recommend pooled regional variogram
#' - Compositional samples outside normal range (95-105)
#' - No samples in zone boundary (if zones provided as raster)
#'
#' @examples
#' \dontrun{
#' # Load sample data and zone polygons
#' samples <- sf::st_read("soil_samples.shp")
#' zones <- sf::st_read("mapunits.shp")
#'
#' # Prepare zone data (spatial join)
#' zone_data <- gc_prepare_zone_data(samples, zones)
#'
#' # Check validation results
#' for (z in zone_data) {
#'   cat("Zone", z$zone_id, ":", z$n_samples, "samples\n")
#'   if (!z$valid) cat("  WARNING:", z$warning_msgs, "\n")
#' }
#'
#' # If zone_id already in samples:
#' zone_data <- gc_prepare_zone_data(samples, zones, zone_id_column = "mapunit_id")
#' }
#'
#' @importFrom sf st_join st_coordinates st_drop_geometry st_as_sf st_bbox
#' @importFrom methods is
#' @export
gc_prepare_zone_data <- function(data,
                                  zones,
                                  zone_id_column = NULL,
                                  min_samples = 3,
                                  verbose = FALSE) {

  # Validate inputs
  if (!methods::is(data, "data.frame")) {
    stop("data must be a data frame or sf object")
  }

  if (!all(c("x", "y") %in% colnames(data))) {
    stop("data must contain columns 'x' and 'y'")
  }

  # Convert data to sf if needed
  if (!methods::is(data, "sf")) {
    data_sf <- sf::st_as_sf(data, coords = c("x", "y"), crs = NA)
  } else {
    data_sf <- data
  }

  # Handle zones: convert raster to vector if needed
  if (methods::is(zones, "SpatRaster")) {
    zones_sf <- .raster_zones_to_sf(zones)
  } else if (methods::is(zones, "sf")) {
    zones_sf <- zones
  } else {
    stop("zones must be an sf object or terra::SpatRaster")
  }

  if (!("zone_id" %in% colnames(zones_sf))) {
    stop("zones must contain a 'zone_id' column")
  }

  # Spatial join if zone_id_column not provided
  if (is.null(zone_id_column)) {
    data_sf <- sf::st_join(data_sf, zones_sf[, "zone_id"])
  } else {
    if (!(zone_id_column %in% colnames(data_sf))) {
      stop("zone_id_column '", zone_id_column, "' not found in data")
    }
    names(data_sf)[names(data_sf) == zone_id_column] <- "zone_id"
  }

  # Remove rows with NA zone_id (samples outside zones)
  orphaned <- sum(is.na(data_sf$zone_id))
  if (orphaned > 0) {
    warning(orphaned, " samples are outside zone boundaries and will be excluded")
    data_sf <- data_sf[!is.na(data_sf$zone_id), ]
  }

  # Group by zone_id
  zone_ids <- unique(data_sf$zone_id)

  zone_data_list <- lapply(zone_ids, function(z_id) {
    z_data <- data_sf[data_sf$zone_id == z_id, ]
    n_samples <- nrow(z_data)
    
    # Extract to data.frame with x,y columns for compatibility
    coords <- sf::st_coordinates(z_data)
    z_data_df <- sf::st_drop_geometry(z_data)
    z_data_df$x <- coords[, 1]
    z_data_df$y <- coords[, 2]

    # Identify compositional columns (numeric, not x/y/zone_id)
    all_cols <- colnames(z_data_df)
    non_spatial_cols <- setdiff(all_cols, c("x", "y", "zone_id"))
    
    # Filter to numeric columns only
    is_numeric <- sapply(z_data_df[, non_spatial_cols, drop = FALSE], is.numeric)
    comp_cols <- non_spatial_cols[is_numeric]

    # Check for missing values
    missing_vals <- sum(is.na(z_data_df[, comp_cols]))
    warn_msgs <- character()

    if (missing_vals > 0) {
      warn_msgs <- c(warn_msgs, paste(missing_vals, "missing values in compositional data"))
    }

    # Check compositional sums (if columns exist)
    if (length(comp_cols) > 0 && nrow(z_data_df) > 0) {
      comp_matrix <- as.matrix(z_data_df[, comp_cols])
      comp_sums <- rowSums(comp_matrix)

      outside_range <- sum(comp_sums < 95 | comp_sums > 105)
      if (outside_range > 0) {
        warn_msgs <- c(
          warn_msgs,
          paste(outside_range, "samples have sums outside 95-105 range")
        )
      }
    }

    # Sample size warnings
    if (n_samples < min_samples) {
      warn_msgs <- c(
        warn_msgs,
        paste("Only", n_samples, "samples (< min_samples =", min_samples, ")")
      )
    } else if (n_samples < 10) {
      warn_msgs <- c(warn_msgs, "Few samples (< 10); recommend pooled regional variogram")
    }

    # Get zone geometry
    z_geom <- zones_sf[zones_sf$zone_id == z_id, ]

    list(
      zone_id = z_id,
      n_samples = n_samples,
      data = z_data_df,
      valid = n_samples >= min_samples,
      bounds = sf::st_bbox(z_geom),
      warning_msgs = warn_msgs
    )
  })

  names(zone_data_list) <- as.character(zone_ids)

  class(zone_data_list) <- c("gc_zone_data", "list")
  return(zone_data_list)
}


#' Fit Geostatistical Models for Multiple Zones
#'
#' Build independent ILR-based geostatistical models for each zone using the
#' standard pipeline: `gc_ilr_params()` → `gc_fit_vgm()` → `gc_ilr_model()`.
#'
#' @param zone_data A list from `gc_prepare_zone_data()` with zone-specific samples.
#' @param model_type Character. Either `"univariate"` (default, numerically stable)
#'   or `"lmc"` (Linear Model of Coregionalization, includes cross-covariance).
#' @param vgm_per_zone Logical. If TRUE (default), fit variogram separately per zone.
#'   If FALSE, fit single regional variogram across all zones combined.
#' @param regional_vgm Optional `gstat::vgm()` object to use as template/prior for
#'   all zones. Useful if zone-specific fitting is unreliable. Ignored if NULL.
#' @param skip_sparse Logical. If TRUE (default), skip zones with < 3 samples and
#'   issue warning. If FALSE, attempt fitting anyway (may fail).
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A list with structure matching `zone_data` (one element per zone),
#'   where each zone element contains:
#'   - `ilr_params`: Output of `gc_ilr_params()`
#'   - `variogram`: Fitted `gstat::vgm()` object
#'   - `model`: Fitted `gstat` model object (from `gc_ilr_model()`)
#'   - `fit_status`: Character, "success", "warning", or "skipped"
#'   - `fit_notes`: Character vector describing any issues
#'
#' Class attribute: `"gc_zone_models"` for S3 method dispatch.
#'
#' @details
#' **Per-Zone Workflow:**
#' 1. Check sample count (if < 3, skip or warn per `skip_sparse`)
#' 2. Extract compositional columns and convert to ILR
#' 3. Fit variogram to zone-specific or regional ILR dimensions
#' 4. Build gstat model
#' 5. Store with metadata
#'
#' **Sparse Zones (< 10 samples):**
#' - If `vgm_per_zone=TRUE` and zone has < 5 samples, issue warning
#' - Recommend switching to `vgm_per_zone=FALSE` (regional variogram)
#' - If `regional_vgm` provided, use it as template
#'
#' **Regional Variogram Option:**
#' When `vgm_per_zone=FALSE`, combines all zone data to fit single variogram.
#' Useful when:
#' - Many small zones with few samples each
#' - Spatial structure is homogeneous across zones
#' - Variogram estimation per zone is unstable
#'
#' @examples
#' \dontrun{
#' # Standard workflow: separate variogram per zone
#' zone_data <- gc_prepare_zone_data(samples, zones)
#' zone_models <- gc_fit_zone_models(zone_data)
#'
#' # Check fit results
#' for (z_name in names(zone_models)) {
#'   zm <- zone_models[[z_name]]
#'   cat("Zone", z_name, ":", zm$fit_status, "\n")
#'   if (length(zm$fit_notes) > 0) cat("  Notes:", paste(zm$fit_notes, collapse="; "), "\n")
#' }
#'
#' # For heterogeneous sample density: use regional variogram
#' zone_models <- gc_fit_zone_models(zone_data, vgm_per_zone = FALSE)
#' }
#'
#' @importFrom sf st_drop_geometry
#' @export
gc_fit_zone_models <- function(zone_data,
                                model_type = "univariate",
                                vgm_per_zone = TRUE,
                                regional_vgm = NULL,
                                skip_sparse = TRUE,
                                verbose = TRUE) {

  if (!inherits(zone_data, "gc_zone_data")) {
    stop("zone_data must be output from gc_prepare_zone_data()")
  }

  zone_ids <- names(zone_data)
  zone_models <- list()

  # If using regional variogram, fit across all zones first
  if (!vgm_per_zone) {
    if (verbose) cat("Fitting regional variogram across all zones...\n")

    # Combine all data
    all_data_list <- lapply(zone_data, function(z) z$data)
    all_data <- do.call(rbind, all_data_list)

    if (nrow(all_data) < 5) {
      stop("Insufficient total samples (", nrow(all_data), " < 5) for regional variogram")
    }

    # Extract compositional columns
    non_spatial_cols <- setdiff(colnames(all_data), c("x", "y", "zone_id", "geometry"))
    is_numeric <- sapply(all_data[, non_spatial_cols, drop = FALSE], is.numeric)
    comp_cols <- non_spatial_cols[is_numeric]

    samples <- sf::st_drop_geometry(all_data[, comp_cols])

    # Compute ILR params from regional data
    regional_ilr_params <- gc_ilr_params(samples)

    # Transform to ILR space for regional variogram
    comp <- compositions::acomp(samples)
    ilr_vals <- compositions::ilr(comp)
    
    all_data_for_vgm <- as.data.frame(ilr_vals)
    names(all_data_for_vgm) <- paste0("ilr", seq_len(ncol(ilr_vals)))
    # Add x, y coordinates from original all_data
    all_data_for_vgm$x <- all_data$x
    all_data_for_vgm$y <- all_data$y

    # Fit regional variogram
    regional_vgm_fit <- gc_fit_vgm(regional_ilr_params, all_data_for_vgm)

    if (is.null(regional_vgm_fit)) {
      stop("Failed to fit regional variogram")
    }
  }

  # Fit per-zone models
  for (z_id in zone_ids) {
    z_data <- zone_data[[z_id]]

    if (verbose) cat("Processing zone", z_id, "(", z_data$n_samples, "samples)...\n")

    # Skip sparse zones if requested
    if (z_data$n_samples < 3) {
      if (skip_sparse) {
        if (verbose) cat("  Skipping (< 3 samples)\n")
        zone_models[[z_id]] <- list(
          fit_status = "skipped",
          fit_notes = c("Fewer than 3 samples")
        )
        next
      }
    }

    # Extract compositional columns
    non_spatial_cols <- setdiff(colnames(z_data$data), c("x", "y", "zone_id", "geometry"))
    is_numeric <- sapply(z_data$data[, non_spatial_cols, drop = FALSE], is.numeric)
    comp_cols <- non_spatial_cols[is_numeric]

    if (length(comp_cols) == 0) {
      zone_models[[z_id]] <- list(
        fit_status = "error",
        fit_notes = c("No compositional columns identified")
      )
      next
    }

    samples <- z_data$data[, comp_cols]

    # Attempt to fit model
    fit_result <- tryCatch(
      {
        # Compute ILR params
        ilr_params <- gc_ilr_params(samples)

        # Transform compositional data to ILR space
        comp <- compositions::acomp(samples)
        ilr_vals <- compositions::ilr(comp)
        
        # Create data frame with ILR columns and spatial coordinates
        z_data_for_vgm <- as.data.frame(ilr_vals)
        # Rename columns to ilr1, ilr2, ...
        names(z_data_for_vgm) <- paste0("ilr", seq_len(ncol(ilr_vals)))
        z_data_for_vgm$x <- z_data$data$x
        z_data_for_vgm$y <- z_data$data$y
        
        # Fit variogram
        if (vgm_per_zone) {
          vgm_fit <- gc_fit_vgm(ilr_params, z_data_for_vgm)
          fit_notes <- character()

          if (z_data$n_samples < 5) {
            fit_notes <- c(fit_notes, "Few samples; variogram may be unreliable")
          }
        } else {
          # Use regional variogram
          vgm_fit <- regional_vgm_fit
          fit_notes <- c("Using regional variogram from pooled data")
        }

        # Build model
        model <- gc_ilr_model(ilr_params, vgm_fit, data = z_data_for_vgm, model_type = model_type)

        list(
          ilr_params = ilr_params,
          variogram = vgm_fit,
          model = model,
          fit_status = "success",
          fit_notes = fit_notes
        )
      },
      error = function(e) {
        list(
          fit_status = "error",
          fit_notes = c(paste("Error:", conditionMessage(e)))
        )
      }
    )

    zone_models[[z_id]] <- fit_result

    if (verbose && fit_result$fit_status != "error" && length(fit_result$fit_notes) > 0) {
      cat("  Status: ", paste(fit_result$fit_notes, collapse = "; "), "\n", sep = "")
    } else if (verbose && fit_result$fit_status == "error") {
      cat("  ERROR:", fit_result$fit_notes[1], "\n")
    }
  }

  class(zone_models) <- c("gc_zone_models", "list")
  return(zone_models)
}


#' Simulate Compositional Data Within Vector or Raster Zones
#'
#' Generate spatial realizations of compositional data for each zone independently,
#' then stitch results into a seamless tessellation covering the full extent.
#'
#' @param zone_models A list from `gc_fit_zone_models()` with fitted models per zone.
#' @param zones An sf object (polygons) or terra::SpatRaster defining zone boundaries.
#' @param extent An sf object or terra::SpatRaster defining the full simulation extent.
#' @param resolution Numeric. Grid cell size (in map units, typically meters).
#' @param nsim Numeric. Number of realizations per zone (default 1).
#' @param nmax Numeric. Maximum number of neighbors for kriging. If NULL (default),
#'   uses all available observations per zone. Smaller values (12-20) improve speed
#'   for large datasets but may introduce edge effects.
#' @param edge_handling Character. Strategy for zone boundaries. Options:
#'   - `"buffer"` (default): Expand grid by 10% buffer before simulation, then clip
#'   - `"clip"`: Strict boundary clip; may lose edge information
#'   - `"overlap"`: Allow 1-cell overlap between zones for smooth transitions
#' @param fill_gaps Logical. If TRUE (default), fill any gaps/slivers with values
#'   from nearest zone via simple nearest-neighbor interpolation.
#' @param validate_constraints Logical. If TRUE (default), validate sum constraints
#'   post-stitching and issue warnings if violated.
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A terra::SpatRaster with layer names `<component>.<zone_id>.sim<N>`.
#'   E.g., for 3 components (sand, silt, clay) and 2 realizations in 4 zones,
#'   output has 24 layers:
#'   `sand.1.sim1, silt.1.sim1, clay.1.sim1, sand.1.sim2, ..., clay.4.sim2`
#'
#' Attributes:
#'   - `zone_extent`: Bounding box of all zones
#'   - `resolution`: Grid resolution used
#'   - `n_zones`: Number of zones simulated
#'   - `n_realizations`: Number of realizations per zone
#'
#' @details
#' **Per-Zone Workflow:**
#' 1. Create simulation grid within zone boundary (plus optional buffer)
#' 2. Call `gc_sim_composition()` for this zone's model
#' 3. Validate sum constraints within zone
#' 4. Clip results to zone boundary (if edge_handling != "overlap")
#'
#' **Stitching Algorithm:**
#' 1. Merge all zone rasters
#' 2. Detect and handle overlaps (if any) using edge_handling strategy
#' 3. Fill remaining gaps from nearest non-NA cell
#' 4. Validate global sum constraints
#'
#' **Performance Notes:**
  # - Zone processing is sequential (not parallelized)
#' - For 100+ zones at 30m resolution: expect 5-10 min on typical machine
#' - Use `nmax=15` for large datasets (> 1000 samples per zone)
#'
#' @examples
#' \dontrun{
#' # Simulate 5 realizations across 3 zones
#' zone_models <- gc_fit_zone_models(zone_data)
#' result <- gc_simulate_zones(
#'   zone_models = zone_models,
#'   zones = zones_sf,
#'   extent = zones_sf,
#'   resolution = 30,
#'   nsim = 5
#' )
#'
#' # Access individual realization
#' terra::plot(result[[1]])  # First component, all zones, first realization
#' }
#'
#' @importFrom sf st_as_sf st_bbox
#' @importFrom terra rast crop merge ext
#' @importFrom methods is
#' @export
gc_simulate_zones <- function(zone_models,
                               zones,
                               extent,
                               resolution,
                               nsim = 1,
                               nmax = NULL,
                               edge_handling = "buffer",
                               buffer_distance = NA,
                               fill_gaps = TRUE,
                               validate_constraints = TRUE,
                               rescale_to_target = TRUE,
                               verbose = TRUE) {

  if (!inherits(zone_models, "gc_zone_models")) {
    stop("zone_models must be output from gc_fit_zone_models()")
  }

  edge_handling <- match.arg(edge_handling, c("buffer", "clip", "overlap"))
  
  # Validate zone topology upfront
  if (methods::is(zones, "sf")) {
    topo_check <- .validate_zone_topology(zones)
    if (verbose && length(topo_check$issues) > 0) {
      for (issue in topo_check$issues) {
        cat("  Topology warning:", issue, "\n")
      }
    }
  }

  # Convert extent to raster template if needed
  if (methods::is(extent, "sf")) {
    bbox <- sf::st_bbox(extent)
    extent_rast <- terra::rast(extent, resolution = resolution)
  } else if (methods::is(extent, "SpatRaster")) {
    extent_rast <- extent
  } else {
    stop("extent must be sf object or terra::SpatRaster")
  }

  # Convert zones if needed
  if (methods::is(zones, "SpatRaster")) {
    zones_sf <- .raster_zones_to_sf(zones)
  } else if (methods::is(zones, "sf")) {
    zones_sf <- zones
  } else {
    stop("zones must be sf object or terra::SpatRaster")
  }

  if (!("zone_id" %in% colnames(zones_sf))) {
    stop("zones must contain 'zone_id' column")
  }

  zone_ids <- names(zone_models)
  zone_rasters <- list()

  if (verbose) cat("Simulating", length(zone_ids), "zones...\n")

  # Simulate per zone
  for (z_id in zone_ids) {
    zm <- zone_models[[z_id]]

    # Skip failed zones
    if (zm$fit_status != "success") {
      if (verbose) cat("  Zone", z_id, ": skipping (", zm$fit_status, ")\n", sep = "")
      next
    }

    if (verbose) cat("  Zone", z_id, ": simulating...\n")

    # Get zone geometry
    z_geom <- zones_sf[zones_sf$zone_id == z_id, ]

    # Create simulation grid for this zone
    z_bbox <- terra::ext(z_geom)
    
    # Apply buffer if requested
    if (edge_handling == "buffer") {
      z_geom_buffered <- .buffer_zone(z_geom, buffer_distance, method = "fixed")
      z_bbox <- terra::ext(z_geom_buffered)
      if (verbose) cat("    (buffered extent)\n")
    }

    # Create grid within zone extent
    z_grid <- terra::rast(
      xmin = z_bbox$xmin, xmax = z_bbox$xmax,
      ymin = z_bbox$ymin, ymax = z_bbox$ymax,
      resolution = resolution,
      crs = terra::crs(extent_rast)
    )

    # Convert raster grid to data frame for gc_sim_composition
    grid_cells <- terra::as.data.frame(z_grid, xy = TRUE, na.rm = FALSE)
    grid_df <- grid_cells[c("x", "y")]

    tryCatch(
      {
        # Simulate for this zone
        z_sim <- gc_sim_composition(
          model = zm$model,
          locations = grid_df,
          nsim = nsim,
          target_names = zm$ilr_params$names,
          nmax = nmax
        )

        # Clip to zone boundary based on edge_handling strategy
        if (edge_handling == "buffer" || edge_handling == "clip") {
          # Clip to original (non-buffered) boundary
          z_sim <- terra::crop(z_sim, terra::vect(zones_sf[zones_sf$zone_id == z_id, ]), mask = TRUE)
        }
        
        # Rescale to target sum if requested
        if (rescale_to_target) {
          z_sim <- .rescale_compositions(z_sim, comp_layer_pattern = "\\.sim", 
                                         target_sum = 100, tolerance = 0.5, verbose = FALSE)
        }

        zone_rasters[[z_id]] <- z_sim
      },
      error = function(e) {
        if (verbose) cat("    ERROR:", conditionMessage(e), "\n")
      }
    )
  }

  if (verbose) cat("Stitching zones...\n")

  # Merge all zone rasters
  if (length(zone_rasters) == 0) {
    stop("No zones successfully simulated")
  }

  merged_rast <- .merge_zone_rasters(zone_rasters, extent_rast, edge_handling)

  # Fill gaps if requested
  if (fill_gaps) {
    if (verbose) cat("Filling gaps...\n")
    merged_rast <- .fill_raster_gaps(merged_rast)
  }

  # Final rescaling
  if (rescale_to_target) {
    if (verbose) cat("Final rescaling to target sum...\n")
    merged_rast <- .rescale_compositions(merged_rast, comp_layer_pattern = "\\.sim",
                                         target_sum = 100, tolerance = 0.5, verbose = FALSE)
  }

  # Validate constraints
  if (validate_constraints) {
    if (verbose) cat("Validating constraints...\n")
    validation <- gc_validate_tessellation(merged_rast, pattern = "\\.sim")
    if (!validation$valid) {
      warning("Tessellation has constraint violations:\n",
        paste(validation$issues, collapse = "\n"))
    }
  }

  # Store metadata
  attr(merged_rast, "zone_extent") <- sf::st_bbox(zones_sf)
  attr(merged_rast, "resolution") <- resolution
  attr(merged_rast, "n_zones") <- length(zone_ids)
  attr(merged_rast, "n_realizations") <- nsim

  if (verbose) cat("Simulation complete.\n")

  return(merged_rast)
}


#' Simulate Compositional Data Within Vector or Raster Zones (High-Level Wrapper)
#'
#' Execute complete zone-based simulation workflow in a single call: data preparation,
#' model fitting, simulation, and stitching.
#'
#' @param zones An sf object (polygons) or terra::SpatRaster with zone boundaries.
#' @param data An sf object or data frame with sample coordinates and compositional data.
#' @param extent An sf object or terra::SpatRaster defining simulation output extent.
#' @param resolution Numeric. Grid cell size in map units (default NA; derive from extent if SpatRaster).
#' @param zone_id_column Character. Column name in `data` containing zone IDs.
#'   If NULL (default), performs spatial join.
#' @param nsim Numeric. Number of realizations per zone (default 1).
#' @param model_type Character. `"univariate"` (default) or `"lmc"`.
#' @param vgm_per_zone Logical. Fit separate variogram per zone (default TRUE).
#' @param edge_handling Character. Zone boundary strategy: `"buffer"`, `"clip"`, `"overlap"`.
#' @param verbose Logical. Print progress (default TRUE).
#'
#' @return terra::SpatRaster with zone tessellation. See `gc_simulate_zones()` for layer naming.
#'
#' @details
#' This is a convenience wrapper that calls:
#' 1. `gc_prepare_zone_data()` - Extract zone-specific samples
#' 2. `gc_fit_zone_models()` - Build per-zone models
#' 3. `gc_simulate_zones()` - Generate and stitch realizations
#'
#' For customization (e.g., regional variogram, skip_sparse=FALSE), call the
#' constituent functions directly.
#'
#' @examples
#' \dontrun{
#' # One-line zonal simulation
#' ssurgo <- sf::st_read("mapunits.shp")
#' samples <- sf::st_read("soil_samples.shp")
#'
#' result <- gc_simulate_by_zones(
#'   zones = ssurgo,
#'   data = samples,
#'   extent = ssurgo,
#'   resolution = 30,
#'   nsim = 3
#' )
#'
#' terra::plot(result)
#' }
#'
#' @export
gc_simulate_by_zones <- function(zones,
                                  data,
                                  extent,
                                  resolution = NA,
                                  zone_id_column = NULL,
                                  nsim = 1,
                                  model_type = "univariate",
                                  vgm_per_zone = TRUE,
                                  edge_handling = "buffer",
                                  verbose = TRUE) {

  if (verbose) cat("=== Zone-Based Compositional Simulation ===\n")

  # Prepare zone data
  if (verbose) cat("1. Preparing zone data...\n")
  zone_data <- gc_prepare_zone_data(data, zones, zone_id_column)

  # Report any warnings
  for (z_name in names(zone_data)) {
    z <- zone_data[[z_name]]
    if (length(z$warning_msgs) > 0 && verbose) {
      cat("  Zone", z_name, ":", paste(z$warning_msgs, collapse = "; "), "\n")
    }
  }

  # Fit models
  if (verbose) cat("2. Fitting zone models...\n")
  zone_models <- gc_fit_zone_models(
    zone_data,
    model_type = model_type,
    vgm_per_zone = vgm_per_zone,
    verbose = verbose
  )

  # Infer resolution if needed
  if (is.na(resolution)) {
    if (methods::is(extent, "SpatRaster")) {
      resolution <- terra::res(extent)[1]
    } else {
      stop("resolution must be specified or extent must be terra::SpatRaster")
    }
  }

  # Simulate zones
  if (verbose) cat("3. Simulating zones and stitching...\n")
  result <- gc_simulate_zones(
    zone_models = zone_models,
    zones = zones,
    extent = extent,
    resolution = resolution,
    nsim = nsim,
    edge_handling = edge_handling,
    verbose = verbose
  )

  if (verbose) cat("=== Complete ===\n")

  return(result)
}


#' Validate Tessellated Compositional Raster
#'
#' Check sum constraints and compositional validity across all cells in a
#' tessellated (zone-based) simulation output.
#'
#' @param raster A terra::SpatRaster with simulated compositional data.
#' @param target_sum Numeric. Target sum for compositions (default 100 for percentages).
#' @param tolerance Numeric. Allowed deviation from target_sum (default 0.1).
#'
#' @return A list containing:
#'   - `valid`: Logical, TRUE if all cells pass constraints
#'   - `summary`: Data frame with validation stats per component
#'   - `violations`: Data frame of cells that violate constraints (indices + values)
#'   - `issues`: Character vector describing validation problems
#'
#' @details
#' Validates:
#' 1. All cells have sum = target_sum +/- tolerance
#' 2. No negative values
#' 3. No missing values (except designated nodata)
#' 4. Patterns consistent with compositional data
#'
#' Reports:
#' - Count of violations per component
#' - Spatial distribution of violations
#' - Suggestions for remediation (rescaling, manual inspection)
#'
#' @examples
#' \dontrun{
#' result <- gc_simulate_by_zones(zones, data, extent, resolution = 30)
#' validation <- gc_validate_tessellation(result)
#' if (!validation$valid) {
#'   cat("Violations detected:\n")
#'   print(validation$issues)
#' }
#' }
#'
#' @export
gc_validate_tessellation <- function(raster,
                                      target_sum = 100,
                                      tolerance = 0.1) {

  if (!methods::is(raster, "SpatRaster")) {
    stop("raster must be a terra::SpatRaster")
  }

  # Group layers by realization
  layer_names <- names(raster)

  # Extract realization indices from layer names
  # Expected format: component.zone_id.sim<N> e.g., sand.1.sim1
  sim_pattern <- "\\.sim(\\d+)$"
  has_sim <- grepl(sim_pattern, layer_names)

  if (!any(has_sim)) {
    # No realization pattern found, treat all as a single group
    unique_realizations <- 1
  } else {
    # Extract just the number after "sim"
    realization_numbers <- as.numeric(gsub("^.*\\.sim(\\d+)$", "\\1", layer_names[has_sim]))
    unique_realizations <- sort(unique(realization_numbers))
  }

  issues <- character()
  all_violations <- data.frame()

  for (sim_idx in unique_realizations) {
    if (any(has_sim)) {
      sim_layers <- layer_names[grepl(paste0("\\.sim", sim_idx, "$"), layer_names)]
    } else {
      sim_layers <- layer_names
    }

    if (length(sim_layers) == 0) next

    # Extract values for this realization
    rast_values <- terra::values(raster[[sim_layers]], na.rm = FALSE)

    # Handle single layer case
    if (is.null(nrow(rast_values))) {
      rast_values <- matrix(rast_values, ncol = 1)
    }

    # Compute row sums
    row_sums <- rowSums(rast_values, na.rm = TRUE)

    # Check sums
    sum_violations <- which(!is.na(row_sums) & 
                           (row_sums < target_sum - tolerance | row_sums > target_sum + tolerance))

    if (length(sum_violations) > 0) {
      issues <- c(
        issues,
        paste("Realization", sim_idx, ":", length(sum_violations),
          "cells violate sum constraint")
      )

      all_violations <- rbind(
        all_violations,
        data.frame(
          cell_id = sum_violations,
          realization = sim_idx,
          row_sum = row_sums[sum_violations]
        )
      )
    }

    # Check negative values
    neg_values <- sum(rast_values < 0, na.rm = TRUE)
    if (neg_values > 0) {
      issues <- c(issues, paste("Realization", sim_idx, ":", neg_values, "negative values"))
    }
  }

  # Summary per component
  summary_df <- data.frame(
    Component = layer_names[!duplicated(gsub("\\.sim\\d+$", "", layer_names))],
    Count = terra::nlyr(raster) / length(unique_realizations),
    Min = NA_real_,
    Max = NA_real_,
    Mean = NA_real_
  )

  for (i in seq_len(nrow(summary_df))) {
    comp_layers <- layer_names[grep(paste0("^", gsub("\\.", "\\\\.", summary_df$Component[i]), "\\."),
                                     layer_names)]
    if (length(comp_layers) > 0) {
      vals <- terra::values(raster[[comp_layers]], na.rm = TRUE)
      summary_df$Min[i] <- min(vals)
      summary_df$Max[i] <- max(vals)
      summary_df$Mean[i] <- mean(vals)
    }
  }

  valid <- length(issues) == 0

  list(
    valid = valid,
    summary = summary_df,
    violations = all_violations,
    issues = issues
  )
}


# ============================================================================
# Helper functions
# ============================================================================

#' Convert Raster Zones to Vector Zones
#'
#' Convert terra::SpatRaster with zone IDs to sf polygon layer.
#'
#' @param zones_raster terra::SpatRaster with cell values = zone_id (0 = nodata)
#'
#' @return sf object (POLYGON) with column zone_id
#'
#' @keywords internal
#' @importFrom terra as.polygons values crs
#' @importFrom sf st_as_sf st_make_valid
#'
.raster_zones_to_sf <- function(zones_raster) {
  # Convert raster to polygons
  zones_poly <- terra::as.polygons(zones_raster)

  # Extract zone_id values
  zones_sf <- sf::st_as_sf(zones_poly)

  # Clean up geometry
  zones_sf <- sf::st_make_valid(zones_sf)

  # Rename zone column to zone_id
  colnames(zones_sf)[1] <- "zone_id"

  return(zones_sf)
}


#' Fill Raster Gaps with Nearest Value
#'
#' Simple gap-filling algorithm using nearest-neighbor interpolation.
#'
#' @param raster terra::SpatRaster with potential NA/nodata cells
#'
#' @return terra::SpatRaster with gaps filled
#'
#' @keywords internal
#' @importFrom terra focal values not.na
#'
# Gap filling is now handled in zonal_edge_handling.R via .fill_raster_gaps()
# which provides more control and parameters for iteration and method selection
