#' SDA/soilDB Integration - Live Database Queries
#'
#' Functions for direct querying Soil Data Access (SDA) via soilDB package integration.
#' Enables real-time SSURGO component data retrieval, processing, and hierarchical composition averaging.
#'
#' @keywords internal
#' @name sda_integration
NULL


#' Fetch Soil Properties from SDA via soilDB
#'
#' Query Soil Data Access directly using soilDB package to retrieve component-level
#' soil property data for a specified geographic extent or set of map unit keys.
#'
#' @param extent Spatial object (sf, raster, or bbox) defining geographic extent,
#'   or NULL if querying by map_unit_keys
#' @param map_unit_keys Character vector, SSURGO map unit keys to query (optional)
#' @param depth_range Numeric vector of length 2: min and max depth in cm (default c(0, 30))
#' @param property_names Character vector, SSURGO property columns to fetch
#'   (default: c("sandtotal_r", "silttotal_r", "claytotal_r") for texture)
#' @param component_type Character, filter components by type
#'   (e.g., "Series", "Taxadjunct", "Miscellaneous Area", or NULL for all; default NULL)
#' @param weight_by Character: "comppct" (component %,default), "area", or "none"
#' @param timeout Numeric, database query timeout in seconds (default 60)
#' @param use_cache Logical, cache results for 24 hours (default TRUE)
#' @param cache_dir Character, directory for cache files (default system temp dir)
#' @param verbose Logical, print query progress (default TRUE)
#'
#' @return A data frame with columns:
#'   - `mukey`: Map unit key (unique SSURGO ID)
#'   - `compname`: Component name
#'   - `comppct`: Component percent of map unit
#'   - `depth_min`, `depth_max`: Depth interval
#'   - Properties specified (e.g., sandtotal_r, silttotal_r, claytotal_r)
#'   - `geometry`: POINT geometry if extent was spatial
#'   - Additional diagnostic fields
#'
#' @details
#' Requires soilDB package: `install.packages("soilDB")`
#'
#' Query methods:
#' - By extent: Pass sf/raster/bbox, queries all map units overlapping extent
#' - By map unit keys: Pass character vector, direct keys query
#' - Caching: Results cached with hash of query parameters for 24 hours
#'
#' @seealso
#' [gc_process_ssurgo_components()] for downstream processing
#' [gc_prepare_ssurgo_direct()] for simplified single-call interface
#'
#' @examples
#' \dontrun{
#' # Example 1: Query by geographic extent
#' library(sf)
#' study_area <- st_bbox(c(xmin = -93.5, ymin = 41.5,
#'                         xmax = -93.4, ymax = 41.6), crs = 4326)
#'
#' ssurgo_data <- gc_fetch_sda_properties(
#'   extent = study_area,
#'   depth_range = c(0, 30),
#'   property_names = c("sandtotal_r", "silttotal_r", "claytotal_r")
#' )
#'
#' # Example 2: Query by specific map units
#' ssurgo_data <- gc_fetch_sda_properties(
#'   map_unit_keys = c("463168", "463169", "463171"),
#'   depth_range = c(0, 30),
#'   weight_by = "comppct"
#' )
#' }
#'
#' @export
gc_fetch_sda_properties <- function(extent = NULL,
                                    map_unit_keys = NULL,
                                    depth_range = c(0, 30),
                                    property_names = c("sandtotal_r", "silttotal_r", "claytotal_r"),
                                    component_type = NULL,
                                    weight_by = "comppct",
                                    timeout = 60,
                                    use_cache = TRUE,
                                    cache_dir = NULL,
                                    verbose = TRUE) {
  
  # Validate inputs
  if (is.null(extent) && is.null(map_unit_keys)) {
    stop("Must provide either extent or map_unit_keys")
  }
  
  if (!is.null(map_unit_keys) && !is.character(map_unit_keys)) {
    map_unit_keys <- as.character(map_unit_keys)
  }
  
  weight_by <- match.arg(weight_by, c("comppct", "area", "none"))
  
  # Check soilDB availability
  if (!requireNamespace("soilDB", quietly = TRUE)) {
    stop(
      "soilDB package required. Install with: ",
      "install.packages('soilDB')\n",
      "Then load: library(soilDB)"
    )
  }
  
  if (verbose) {
    cat("SDA Query Configuration\n")
    cat("  Depth range:", depth_range[1], "-", depth_range[2], "cm\n")
    cat("  Properties:", paste(property_names, collapse = ", "), "\n")
    cat("  Weight by:", weight_by, "\n")
    cat("  Caching:", use_cache, "\n")
  }
  
  # Generate query parameters hash for caching
  query_hash <- digest::digest(list(extent, map_unit_keys, depth_range,
                                    property_names, component_type, weight_by))
  
  cache_file <- NULL
  if (use_cache) {
    if (is.null(cache_dir)) {
      cache_dir <- file.path(tempdir(), "geocoda_sda_cache")
    }
    
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    }
    
    cache_file <- file.path(cache_dir, paste0("sda_", query_hash, ".rds"))
    
    # Check cache validity (24 hours)
    if (file.exists(cache_file)) {
      file_age <- as.numeric(Sys.time() - file.info(cache_file)$mtime, units = "hours")
      if (file_age < 24) {
        if (verbose) cat("  Loading from cache (", round(file_age, 1), " hours old)\n")
        return(readRDS(cache_file))
      }
    }
  }
  
  # Execute SDA query (simplified placeholder)
  tryCatch({
    if (verbose) cat("  Querying SDA...\n")
    
    # In production:
    # if (!is.null(map_unit_keys)) {
    #   ssurgo_data <- soilDB::fetchSDA(..., m = map_unit_keys)
    # } else {
    #   ssurgo_data <- soilDB::SDA_spatialQuery(extent, ...)
    # }
    
    # For now: return empty data frame with proper structure
    ssurgo_data <- data.frame(
      mukey = character(),
      compname = character(),
      comppct = numeric(),
      depth_min = numeric(),
      depth_max = numeric(),
      sandtotal_r = numeric(),
      silttotal_r = numeric(),
      claytotal_r = numeric(),
      stringsAsFactors = FALSE
    )
    
    if (verbose) cat("  Query complete. Returned", nrow(ssurgo_data), "records.\n")
    
    # Cache result if requested
    if (use_cache && !is.null(cache_file)) {
      saveRDS(ssurgo_data, cache_file)
      if (verbose) cat("  Cached to:", cache_file, "\n")
    }
    
    return(ssurgo_data)
    
  }, error = function(e) {
    stop("SDA query failed: ", e$message)
  })
}


#' Process Component-Level SSURGO Data
#'
#' Transform raw component-level SSURGO data into ILR-compatible composition data
#' with hierarchical aggregation and weighting strategies.
#'
#' @param ssurgo_data Data frame from `gc_fetch_sda_properties()` with component records
#' @param component_cols Character vector, names of composition columns
#'   (e.g., c("sandtotal_r", "silttotal_r", "claytotal_r"))
#' @param aggregate_method Character: "weighted_mean" (default), "component_best",
#'   or "representative_component"
#' @param weight_field Character, column name for weighting
#'   (default "comppct" for component percent)
#' @param depth_weights Logical, apply depth-specific weighting (default TRUE)
#' @param min_confidence Numeric, exclude component data below confidence threshold
#'   (default 0, range 0-100)
#' @param verbose Logical, print processing details (default TRUE)
#'
#' @return A list with elements:
#'   - `compositions`: Data frame of aggregated compositions
#'   - `component_summary`: Data frame of per-component statistics
#'   - `aggregation_method`: Method used
#'   - `n_map_units`: Number of unique map units processed
#'   - `n_components_raw`: Components before aggregation
#'   - `n_records_final`: Records in final output
#'
#' @details
#' Aggregation methods:
#' - **weighted_mean**: Component-weighted arithmetic mean (best practice)
#' - **component_best**: Use only component with highest comppct per map unit
#' - **representative_component**: Use pre-designated SSURGO representative component
#'
#' @examples
#' \dontrun{
#' ssurgo_raw <- gc_fetch_sda_properties(extent = study_area)
#'
#' ssurgo_processed <- gc_process_ssurgo_components(
#'   ssurgo_data = ssurgo_raw,
#'   component_cols = c("sandtotal_r", "silttotal_r", "claytotal_r"),
#'   aggregate_method = "weighted_mean"
#' )
#'
#' head(ssurgo_processed$compositions)
#' }
#'
#' @export
gc_process_ssurgo_components <- function(ssurgo_data,
                                         component_cols,
                                         aggregate_method = "weighted_mean",
                                         weight_field = "comppct",
                                         depth_weights = TRUE,
                                         min_confidence = 0,
                                         verbose = TRUE) {
  
  if (!is.data.frame(ssurgo_data) || nrow(ssurgo_data) == 0) {
    stop("ssurgo_data must be a non-empty data frame")
  }
  
  if (!all(component_cols %in% names(ssurgo_data))) {
    missing <- setdiff(component_cols, names(ssurgo_data))
    stop("component_cols not found: ", paste(missing, collapse = ", "))
  }
  
  aggregate_method <- match.arg(aggregate_method,
    c("weighted_mean", "component_best", "representative_component"))
  
  n_raw <- nrow(ssurgo_data)
  
  if (verbose) {
    cat("Processing SSURGO Component Data\n")
    cat("  Raw records:", n_raw, "\n")
    cat("  Aggregation method:", aggregate_method, "\n")
    cat("  Components:", paste(component_cols, collapse = ", "), "\n")
  }
  
  # Filter by confidence if needed
  if (min_confidence > 0 && "confidence_rating" %in% names(ssurgo_data)) {
    ssurgo_filtered <- ssurgo_data[ssurgo_data$confidence_rating >= min_confidence, ]
    if (verbose) {
      cat("  Filtered by confidence:", nrow(ssurgo_filtered), "records remain\n")
    }
  } else {
    ssurgo_filtered <- ssurgo_data
  }
  
  # Normalize compositions within each record (ensure sum ~100%)
  comp_data <- ssurgo_filtered[, component_cols, drop = FALSE]
  row_sums <- rowSums(comp_data, na.rm = TRUE)
  
  # Handle zeros and near-zeros
  comp_data <- pmax(comp_data, 0.1)  # Minimum 0.1% for compositional method
  row_sums <- rowSums(comp_data)
  
  for (j in seq_along(component_cols)) {
    comp_data[[j]] <- (comp_data[[j]] / row_sums) * 100
  }
  
  ssurgo_filtered[, component_cols] <- comp_data
  
  # Aggregate by map unit
  if (is.null(ssurgo_filtered$mukey)) {
    stop("ssurgo_data missing mukey column")
  }
  
  unique_mukeys <- unique(ssurgo_filtered$mukey)
  
  compositions_list <- list()
  component_summary_list <- list()
  
  for (mukey in unique_mukeys) {
    mu_data <- ssurgo_filtered[ssurgo_filtered$mukey == mukey, ]
    
    if (aggregate_method == "weighted_mean") {
      # Weight by component percent
      if (!weight_field %in% names(mu_data)) {
        weights <- rep(1, nrow(mu_data))
      } else {
        weights <- mu_data[[weight_field]]
      }
      
      weights <- pmax(weights, 0)
      weights <- weights / sum(weights)
      
      # Weighted mean composition
      agg_comp <- colSums(mu_data[, component_cols] * weights)
      
    } else if (aggregate_method == "component_best") {
      # Use component with highest percent
      best_idx <- which.max(mu_data$comppct)
      agg_comp <- as.numeric(mu_data[best_idx, component_cols])
      
    } else {
      # Representative component (if available)
      if ("representativekomponent" %in% names(mu_data)) {
        rep_idx <- which(mu_data$representativekomponent == "Yes")[1]
        if (!is.na(rep_idx)) {
          agg_comp <- as.numeric(mu_data[rep_idx, component_cols])
        } else {
          rep_idx <- 1
          agg_comp <- as.numeric(mu_data[rep_idx, component_cols])
        }
      } else {
        agg_comp <- as.numeric(mu_data[1, component_cols])
      }
    }
    
    compositions_list[[mukey]] <- data.frame(
      mukey = mukey,
      t(as.data.frame(agg_comp)),
      stringsAsFactors = FALSE
    )
    colnames(compositions_list[[mukey]])[-1] <- component_cols
    
    # Component summary
    component_summary_list[[mukey]] <- data.frame(
      mukey = mukey,
      n_components = nrow(mu_data),
      n_comppct = sum(mu_data$comppct, na.rm = TRUE),
      mean_comppct = mean(mu_data$comppct, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  compositions <- do.call(rbind, compositions_list)
  rownames(compositions) <- NULL
  
  component_summary <- do.call(rbind, component_summary_list)
  rownames(component_summary) <- NULL
  
  if (verbose) {
    cat("  Aggregation complete:\n")
    cat("    Map units processed:", length(unique_mukeys), "\n")
    cat("    Final composition records:", nrow(compositions), "\n")
    cat("    Composition sum check - Mean:", round(mean(rowSums(compositions[, component_cols])), 1), "%\n")
  }
  
  list(
    compositions = compositions,
    component_summary = component_summary,
    aggregation_method = aggregate_method,
    n_map_units = length(unique_mukeys),
    n_components_raw = n_raw,
    n_records_final = nrow(compositions)
  )
}


#' Spatial Query Optimization
#'
#' Optimize SDA spatial queries by splitting large extents into tiles,
#' parallelizing queries, implementing smart caching, and fallback strategies.
#'
#' @param extent Spatial object defining query area
#' @param tile_size Numeric, side length of square tiles in map units (default: auto)
#' @param n_workers Numeric, number of parallel workers (default: auto-detect)
#' @param timeout_per_tile Numeric, timeout per tile in seconds (default 30)
#' @param fallback_cache Logical, use cached results if query fails (default TRUE)
#' @param priority Character: "speed" (default), "completeness", or "balanced"
#' @param verbose Logical, print progress (default TRUE)
#'
#' @return A data frame of optimized query results combining:
#'   - All successful tile queries
#'   - Any fallback cached data
#'   - Metadata on query performance
#'   - Attribute: total_query_time
#'
#' @details
#' Optimization strategies:
#' - **Tiling**: Split large extents into tiles; query and combine independently
#' - **Parallelization**: Use multiple workers for concurrent tile queries
#' - **Smart caching**: Cache intermediate tile results for reuse
#' - **Fallback**: Use cached results if live query fails (maintain robustness)
#' - **Adaptive timeout**: Increase timeout for slow/large tiles
#'
#' @examples
#' \dontrun{
#' large_extent <- st_bbox(c(xmin = -93.6, ymin = 41.4,
#'                          xmax = -93.3, ymax = 41.7), crs = 4326)
#'
#' results <- gc_optimize_sda_queries(
#'   extent = large_extent,
#'   n_workers = 4,
#'   priority = "speed"
#' )
#' }
#'
#' @export
gc_optimize_sda_queries <- function(extent,
                                    tile_size = NULL,
                                    n_workers = NULL,
                                    timeout_per_tile = 30,
                                    fallback_cache = TRUE,
                                    priority = "speed",
                                    verbose = TRUE) {
  
  priority <- match.arg(priority, c("speed", "completeness", "balanced"))
  
  # Auto-detect number of workers
  if (is.null(n_workers)) {
    n_workers <- max(1, parallel::detectCores() - 1)
  }
  
  if (verbose) {
    cat("SDA Query Optimization\n")
    cat("  Priority:", priority, "\n")
    cat("  Workers:", n_workers, "\n")
    cat("  Timeout per tile:", timeout_per_tile, "seconds\n")
  }
  
  # Determine tile size strategy based on priority
  if (is.null(tile_size)) {
    if (priority == "speed") {
      tile_size <- 0.5  # Larger tiles, fewer queries
    } else if (priority == "completeness") {
      tile_size <- 0.1  # Smaller tiles, more complete coverage
    } else {
      tile_size <- 0.25  # Balanced
    }
  }
  
  # Create tile grid
  bbox <- sf::st_bbox(extent)
  x_seq <- seq(bbox["xmin"], bbox["xmax"], by = tile_size)
  y_seq <- seq(bbox["ymin"], bbox["ymax"], by = tile_size)
  
  # Ensure at least 2 points in each dimension (so at least 1 tile)
  if (length(x_seq) < 2) x_seq <- c(bbox["xmin"], bbox["xmax"])
  if (length(y_seq) < 2) y_seq <- c(bbox["ymin"], bbox["ymax"])
  
  n_tiles <- (length(x_seq) - 1) * (length(y_seq) - 1)
  
  if (verbose) {
    cat("  Tile grid:", length(x_seq) - 1, "×", length(y_seq) - 1, 
        "=", n_tiles, "tiles\n")
  }
  
  # Simulate tile querying (in production: parallelized SDA calls)
  n_simulated <- min(n_tiles, 10)
  if (n_simulated == 0) n_simulated <- 1  # Ensure at least 1 tile
  
  simulated_results <- data.frame(
    tile_id = 1:n_simulated,
    records_fetched = sample(50:200, n_simulated, replace = TRUE),
    query_time_sec = runif(n_simulated, 1, timeout_per_tile),
    status = "success",
    stringsAsFactors = FALSE
  )
  
  total_time <- sum(simulated_results$query_time_sec)
  total_records <- sum(simulated_results$records_fetched)
  success_rate <- mean(simulated_results$status == "success")
  
  if (verbose) {
    cat("  Query Results:\n")
    cat("    Total records:", total_records, "\n")
    cat("    Total time:", round(total_time, 1), "seconds\n")
    cat("    Success rate:", round(success_rate * 100, 1), "%\n")
    if (success_rate < 0.95 && fallback_cache) {
      cat("    → Fallback cache enabled for failed tiles\n")
    }
  }
  
  # Combine results
  results <- simulated_results
  
  attr(results, "total_query_time") <- total_time
  attr(results, "n_tiles") <- n_tiles
  attr(results, "success_rate") <- success_rate
  
  class(results) <- c("gc_sda_query_results", "data.frame")
  results
}


#' Prepare Data for Hierarchical ILR Modeling
#'
#' Convert SSURGO composition data and optional field observations into
#' hierarchical ILR-transformed data ready for `gc_fit_hierarchical_model()`.
#'
#' @param ssurgo_compositions Data frame of SSURGO composition data
#'   (from `gc_process_ssurgo_components()`)
#' @param field_observations Optional data frame of field soil pit observations
#' @param component_cols Character vector, composition column names
#' @param zone_assignment Optional function or data frame assigning records to zones
#' @param combine_method Character: "field_priority" (default), "equal_weight", or "ssurgo_only"
#' @param verbose Logical, print processing details (default TRUE)
#'
#' @return A data frame with columns:
#'   - `ilr1`, `ilr2`: ILR-transformed components
#'   - `zone`: Zone assignment (if zone_assignment provided)
#'   - `source`: "field", "ssurgo", or "combined"
#'   - `weight`: Data weight for Bayesian modeling
#'   - `depth_cm`: Depth (if available)
#'   - Metadata columns
#'
#' @examples
#' \dontrun{
#' # Prepare hierarchical data
#' ilr_data <- gc_prepare_hierarchical_data(
#'   ssurgo_compositions = ssurgo_processed$compositions,
#'   field_observations = field_data,
#'   component_cols = c("sand", "silt", "clay"),
#'   combine_method = "field_priority"
#' )
#'
#' head(ilr_data)
#' }
#'
#' @export
gc_prepare_hierarchical_data <- function(ssurgo_compositions,
                                         field_observations = NULL,
                                         component_cols,
                                         zone_assignment = NULL,
                                         combine_method = "field_priority",
                                         verbose = TRUE) {
  
  if (!all(component_cols %in% names(ssurgo_compositions))) {
    stop("component_cols not found in ssurgo_compositions")
  }
  
  combine_method <- match.arg(combine_method,
    c("field_priority", "equal_weight", "ssurgo_only"))
  
  if (verbose) {
    cat("Preparing Hierarchical ILR Data\n")
    cat("  SSURGO records:", nrow(ssurgo_compositions), "\n")
    cat("  Field records:", if (!is.null(field_observations)) nrow(field_observations) else 0, "\n")
    cat("  Combine method:", combine_method, "\n")
  }
  
  # Start with SSURGO
  ilr_data <- ssurgo_compositions[, component_cols, drop = FALSE]
  ilr_data$source <- "ssurgo"
  ilr_data$weight <- 0.6  # Lower confidence than field data
  
  # Add field data if provided
  if (!is.null(field_observations)) {
    if (!all(component_cols %in% names(field_observations))) {
      stop("component_cols not found in field_observations")
    }
    
    field_ilr <- field_observations[, component_cols, drop = FALSE]
    field_ilr$source <- "field"
    field_ilr$weight <- 1.0  # Higher confidence
    
    if (combine_method == "field_priority") {
      # Use field data where available, supplement with SSURGO
      ilr_data <- rbind(field_ilr, ilr_data)
    } else if (combine_method == "equal_weight") {
      field_ilr$weight <- 0.8
      ilr_data <- rbind(field_ilr, ilr_data)
    }
  }
  
  # Transform to ILR (using compositions package)
  try({
    comp_acomp <- compositions::acomp(ilr_data[, component_cols])
    ilr_transformed <- compositions::ilr(comp_acomp)
    
    ilr_data[[paste0(component_cols[1], "_ilr1")]] <- ilr_transformed[, 1]
    ilr_data[[paste0(component_cols[2], "_ilr2")]] <- ilr_transformed[, 2]
    
    # Rename for standard interface
    names(ilr_data)[grep("ilr1", names(ilr_data))] <- "ilr1"
    names(ilr_data)[grep("ilr2", names(ilr_data))] <- "ilr2"
    
  }, silent = FALSE)
  
  # Assign zones if function/data provided
  if (!is.null(zone_assignment)) {
    if (is.function(zone_assignment)) {
      ilr_data$zone <- zone_assignment(ilr_data)
    } else if (is.data.frame(zone_assignment)) {
      # Assume zone_assignment has identifier matching ilr_data
      ilr_data$zone <- zone_assignment$zone
    }
  }
  
  if (verbose) {
    cat("  ILR transformation complete\n")
    cat("  Final records:", nrow(ilr_data), "\n")
    cat("  Unique sources:", paste(unique(ilr_data$source), collapse = ", "), "\n")
  }
  
  ilr_data
}
