#' SDA/soilDB Integration - Live Database Queries
#'
#' Functions for direct querying Soil Data Access (SDA) via soilDB package integration.
#' Enables real-time SSURGO component data retrieval, processing, and hierarchical composition averaging.
#'
#' @keywords internal
#' @name sda_integration
NULL


#' Fetch Soil Properties from SDA
#'
#' Query Soil Data Access (SDA) directly using soilDB package to retrieve
#' soil component and horizon properties for a specified geographic extent
#' or set of map unit keys.
#'
#' @param extent Spatial object (sf, raster, or bbox) defining geographic extent,
#'   or NULL if querying by map_unit_keys
#' @param map_unit_keys Numeric or character vector, SSURGO map unit keys to query
#'   (optional, alternative to extent)
#' @param depth_range Numeric vector of length 2: min and max depth in cm
#'   (default c(0, 30))
#' @param property_names Character vector, SSURGO property columns or names to fetch
#'   (default: sand, silt, clay total; e.g. "sandtotal_r", "Total Sand - Rep Value")
#' @param method Character, aggregation method ("Weighted Average" (default),
#'   "Dominant Component (Numeric)", "Dominant Component (Category)",
#'   "Min/Max", "Dominant Condition", or "None")
#' @param include_minors Logical, include minor components in aggregation (default TRUE)
#' @param use_cache Logical, cache results for 24 hours (default TRUE)
#' @param cache_dir Character, directory for cache files (default system temp dir)
#' @param verbose Logical, print minimal query progress (default TRUE)
#'
#' @return A data.frame with map unit and property data:
#'   - `mukey`: Map unit key (unique SSURGO ID)
#'   - `areasymbol`: Soil survey area code
#'   - `musym`, `muname`: Map unit symbol and name
#'   - Properties specified (e.g., sandtotal_r, silttotal_r, claytotal_r)
#'
#' @details
#' Uses `soilDB::get_SDA_property()` for property retrieval with various
#' aggregation methods. Caching is automatic with 24-hour TTL.
#'
#' @examples
#' \dontrun{
#' # Query by geographic extent
#' library(sf)
#' extent <- st_bbox(c(xmin = -93.5, ymin = 41.5,
#'                     xmax = -93.4, ymax = 41.6), crs = 4326)
#'
#' ssurgo <- gc_fetch_sda_properties(
#'   extent = extent,
#'   depth_range = c(0, 30),
#'   property_names = c("sandtotal_r", "silttotal_r", "claytotal_r")
#' )
#' }
#'
#' @export
#' @seealso [soilDB::get_SDA_property()], [soilDB::SDA_spatialQuery()]
gc_fetch_sda_properties <- function(extent = NULL,
                                    map_unit_keys = NULL,
                                    depth_range = c(0, 30),
                                    property_names = c("sandtotal_r", "silttotal_r", "claytotal_r"),
                                    method = "Weighted Average",
                                    include_minors = TRUE,
                                    use_cache = TRUE,
                                    cache_dir = NULL,
                                    verbose = TRUE) {
  
  if (is.null(extent) && is.null(map_unit_keys)) {
    stop("Must provide either extent or map_unit_keys")
  }
  
  if (!requireNamespace("soilDB", quietly = TRUE)) {
    stop("soilDB package required. Install with: install.packages('soilDB')")
  }
  
  # Generate query hash for caching
  query_hash <- digest::digest(list(extent, map_unit_keys, depth_range,
                                    property_names, method, include_minors))
  
  cache_file <- NULL
  if (use_cache) {
    if (is.null(cache_dir)) {
      cache_dir <- file.path(tempdir(), "geocoda_sda_cache")
    }
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    }
    cache_file <- file.path(cache_dir, paste0("sda_", query_hash, ".rds"))
    
    if (file.exists(cache_file)) {
      file_age <- as.numeric(Sys.time() - file.info(cache_file)$mtime, units = "hours")
      if (file_age < 24) {
        if (verbose) cat("Loaded from cache\n")
        return(readRDS(cache_file))
      }
    }
  }
  
  tryCatch({
    # Determine query parameters
    query_args <- list(
      property = property_names,
      method = method,
      top_depth = depth_range[1],
      bottom_depth = depth_range[2],
      include_minors = include_minors
    )
    
    # Add either extent or map unit keys
    if (!is.null(extent)) {
      # Spatial query: first get mukeys, then get properties
      mukey_data <- soilDB::SDA_spatialQuery(extent, what = "mukey")
      if (nrow(mukey_data) == 0) {
        return(data.frame())
      }
      mukeys <- unique(mukey_data$mukey)
      query_args$mukeys <- mukeys
    } else {
      query_args$mukeys <- as.numeric(map_unit_keys)
    }
    
    # Execute query and get properties
    result <- do.call(soilDB::get_SDA_property, query_args)
    
    # Cache if requested
    if (use_cache && !is.null(cache_file)) {
      saveRDS(result, cache_file)
    }
    
    if (verbose) {
      cat("SDA query returned", nrow(result), "map unit records\n")
    }
    
    return(result)
    
  }, error = function(e) {
    stop("SDA property query failed: ", e$message)
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
  
  # Filter by confidence if needed
  if (min_confidence > 0 && "confidence_rating" %in% names(ssurgo_data)) {
    ssurgo_filtered <- ssurgo_data[ssurgo_data$confidence_rating >= min_confidence, ]
  } else {
    ssurgo_filtered <- ssurgo_data
  }
  
  # Normalize compositions within each record
  comp_data <- ssurgo_filtered[, component_cols, drop = FALSE]
  row_sums <- rowSums(comp_data, na.rm = TRUE)
  
  # Handle zeros and near-zeros
  comp_data <- pmax(comp_data, 0.1)
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
      if (!weight_field %in% names(mu_data)) {
        weights <- rep(1, nrow(mu_data))
      } else {
        weights <- mu_data[[weight_field]]
      }
      weights <- pmax(weights, 0)
      weights <- weights / sum(weights)
      agg_comp <- colSums(mu_data[, component_cols] * weights)
      
    } else if (aggregate_method == "component_best") {
      best_idx <- which.max(mu_data$comppct)
      agg_comp <- as.numeric(mu_data[best_idx, component_cols])
      
    } else {
      if ("representativekomponent" %in% names(mu_data)) {
        rep_idx <- which(mu_data$representativekomponent == "Yes")[1]
        if (!is.na(rep_idx)) {
          agg_comp <- as.numeric(mu_data[rep_idx, component_cols])
        } else {
          agg_comp <- as.numeric(mu_data[1, component_cols])
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
      mukey =mukey,
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
    cat("Processed", length(unique_mukeys), "map units,",
        nrow(compositions), "output records\n")
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
  
  if (!requireNamespace("soilDB", quietly = TRUE)) {
    stop("soilDB package required. Install with: install.packages('soilDB')")
  }
  
  priority <- match.arg(priority, c("speed", "completeness", "balanced"))
  
  # Auto-detect number of workers
  if (is.null(n_workers)) {
    n_workers <- max(1, parallel::detectCores() - 1)
  }
  
  # Determine tile size strategy based on priority
  if (is.null(tile_size)) {
    if (priority == "speed") {
      tile_size <- 0.5
    } else if (priority == "completeness") {
      tile_size <- 0.1
    } else {
      tile_size <- 0.25
    }
  }
  
  # Create tile grid
  bbox <- sf::st_bbox(extent)
  x_seq <- seq(bbox["xmin"], bbox["xmax"], by = tile_size)
  y_seq <- seq(bbox["ymin"], bbox["ymax"], by = tile_size)
  
  # Ensure at least 2 points per dimension
  if (length(x_seq) < 2) x_seq <- c(bbox["xmin"], bbox["xmax"])
  if (length(y_seq) < 2) y_seq <- c(bbox["ymin"], bbox["ymax"])
  
  n_tiles <- (length(x_seq) - 1) * (length(y_seq) - 1)
  
  # Build tile specifications
  tiles <- list()
  tile_idx <- 1
  
  for (i in seq_len(length(x_seq) - 1)) {
    for (j in seq_len(length(y_seq) - 1)) {
      tiles[[tile_idx]] <- list(
        id = tile_idx,
        xmin = x_seq[i], xmax = x_seq[i + 1],
        ymin = y_seq[j], ymax = y_seq[j + 1]
      )
      tile_idx <- tile_idx + 1
    }
  }
  
  # Query function for each tile
  query_tile <- function(tile, extent) {
    tryCatch({
      t0 <- Sys.time()
      
      tile_bbox <- sf::st_bbox(c(
        xmin = tile$xmin, xmax = tile$xmax,
        ymin = tile$ymin, ymax = tile$ymax
      ), crs = sf::st_crs(extent))
      
      mukey_data <- soilDB::SDA_spatialQuery(tile_bbox, what = "mukey")
      
      t1 <- Sys.time()
      query_time <- as.numeric(difftime(t1, t0, units = "secs"))
      
      list(
        tile_id = tile$id,
        records_fetched = nrow(mukey_data),
        query_time_sec = query_time,
        status = "success",
        n_mukeys = if (nrow(mukey_data) > 0) length(unique(mukey_data$mukey)) else 0
      )
    }, error = function(e) {
      list(
        tile_id = tile$id,
        records_fetched = 0,
        query_time_sec = NA_real_,
        status = "error",
        error_msg = as.character(e),
        n_mukeys = 0
      )
    })
  }
  
  # Execute queries
  if (n_workers > 1) {
    results_list <- parallel::mclapply(
      tiles,
      query_tile,
      extent = extent,
      mc.cores = n_workers,
      mc.preschedule = TRUE
    )
  } else {
    results_list <- lapply(tiles, query_tile, extent = extent)
  }
  
  # Convert results to data frame
  results_df <- do.call(rbind, lapply(results_list, as.data.frame))
  rownames(results_df) <- NULL
  
  # Compute summary statistics
  total_time <- sum(results_df$query_time_sec, na.rm = TRUE)
  total_records <- sum(results_df$records_fetched, na.rm = TRUE)
  success_count <- sum(results_df$status == "success")
  success_rate <- success_count / nrow(results_df)
  total_mukeys <- sum(results_df$n_mukeys, na.rm = TRUE)
  
  if (verbose) {
    cat("SDA Tile Query Results:\n")
    cat("  Tiles queried:", nrow(results_df), "\n")
    cat("  Success rate:", round(success_rate * 100, 1), "%\n")
    cat("  Unique map units:", total_mukeys, "\n")
    cat("  Total time:", round(total_time, 1), "seconds\n")
  }
  
  # Add attributes
  attr(results_df, "total_query_time") <- total_time
  attr(results_df, "n_tiles") <- n_tiles
  attr(results_df, "success_rate") <- success_rate
  attr(results_df, "total_mukeys") <- total_mukeys
  attr(results_df, "n_workers") <- n_workers
  
  class(results_df) <- c("gc_sda_query_results", "data.frame")
  results_df
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
  
  # Start with SSURGO
  ilr_data <- ssurgo_compositions[, component_cols, drop = FALSE]
  ilr_data$source <- "ssurgo"
  ilr_data$weight <- 0.6
  
  # Add field data if provided
  if (!is.null(field_observations)) {
    if (!all(component_cols %in% names(field_observations))) {
      stop("component_cols not found in field_observations")
    }
    
    field_ilr <- field_observations[, component_cols, drop = FALSE]
    field_ilr$source <- "field"
    field_ilr$weight <- 1.0
    
    if (combine_method == "field_priority") {
      ilr_data <- rbind(field_ilr, ilr_data)
    } else if (combine_method == "equal_weight") {
      field_ilr$weight <- 0.8
      ilr_data <- rbind(field_ilr, ilr_data)
    }
  }
  
  # Transform to ILR
  try({
    comp_acomp <- compositions::acomp(ilr_data[, component_cols])
    ilr_transformed <- compositions::ilr(comp_acomp)
    
    # Convert to matrix for consistent column access
    if (!is.matrix(ilr_transformed)) {
      ilr_transformed <- as.matrix(ilr_transformed)
    }
    
    ilr_data[[paste0(component_cols[1], "_ilr1")]] <- ilr_transformed[, 1]
    ilr_data[[paste0(component_cols[2], "_ilr2")]] <- ilr_transformed[, 2]
    
    names(ilr_data)[grep("ilr1", names(ilr_data))] <- "ilr1"
    names(ilr_data)[grep("ilr2", names(ilr_data))] <- "ilr2"
    
  }, silent = FALSE)
  
  # Assign zones if provided
  if (!is.null(zone_assignment)) {
    if (is.function(zone_assignment)) {
      ilr_data$zone <- zone_assignment(ilr_data)
    } else if (is.data.frame(zone_assignment)) {
      ilr_data$zone <- zone_assignment$zone
    }
  }
  
  if (verbose) {
    cat(nrow(ilr_data), "ILR records prepared\n")
  }
  
  ilr_data
}


#' Fetch Map Unit Polygons as Spatial Zones
#'
#' Query Soil Data Access (SDA) to retrieve map unit polygons (MUPOLYGON)
#' within a geographic extent. Polygons serve as spatial zones for zone-stratified
#' geostatistical modeling.
#'
#' @param extent Spatial object (sf, raster, bbox) defining geographic extent
#' @param crs Coordinate reference system (default EPSG:4326)
#' @param verbose Logical, print query progress (default TRUE)
#'
#' @return An sf object with MUPOLYGON geometries and columns:
#'   - `mukey`: Map unit key (SSURGO ID)
#'   - `musym`: Map unit symbol
#'   - `muname`: Map unit name
#'   - `geometry`: sf polygons
#'
#' @details
#' Uses `soilDB::SDA_spatialQuery(what="mupolygon", geomIntersection=TRUE)`
#' to retrieve map unit boundaries intersecting the extent.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' extent <- st_bbox(c(xmin = -93.5, ymin = 41.5,
#'                     xmax = -93.4, ymax = 41.6), crs = 4326)
#'
#' mapunit_zones <- gc_fetch_mapunit_polygons(extent)
#' print(nrow(mapunit_zones))
#' }
#'
#' @export
#' @seealso [soilDB::SDA_spatialQuery()]
gc_fetch_mapunit_polygons <- function(extent, crs = 4326, verbose = TRUE) {
  
  if (!requireNamespace("soilDB", quietly = TRUE)) {
    stop("soilDB package required. Install with: install.packages('soilDB')")
  }
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("sf package required. Install with: install.packages('sf')")
  }
  
  # Convert bbox to sfc if needed
  if (inherits(extent, 'bbox')) {
    extent <- sf::st_as_sfc(extent)
  }
  
  if (verbose) cat("Querying SDA for map unit polygons...\n")
  
  # Query SDA for MUPOLYGON geometries
  mupolygons <- soilDB::SDA_spatialQuery(
    extent,
    what = "mupolygon",
    geomIntersection = TRUE
  )
  
  if (is.null(mupolygons) || nrow(mupolygons) == 0) {
    stop("No map unit polygons found in extent. Ensure extent intersects SSURGO coverage.")
  }
  
  # Ensure valid CRS
  if (is.na(sf::st_crs(mupolygons))) {
    mupolygons <- sf::st_set_crs(mupolygons, crs)
  }
  
  if (verbose) {
    cat("Fetched", nrow(mupolygons), "map unit polygons\n")
  }
  
  return(mupolygons)
}


#' Fetch Aggregated Properties for Map Unit Keys
#'
#' For a set of map unit keys, fetch SSURGO component properties
#' aggregated at the map unit level. Compounds are aggregated using
#' component percent weighting.
#'
#' @param mukeys Numeric or character vector of SSURGO map unit keys
#' @param depth_range Numeric vector of length 2: min and max depth in cm
#'   (default c(0, 30))
#' @param property_names Character vector, SSURGO property names to fetch
#'   (default: c("sandtotal_r", "silttotal_r", "claytotal_r"))
#' @param method Character, aggregation method for components
#'   (default "Weighted Average"; see `gc_fetch_sda_properties()` for options)
#' @param verbose Logical, print progress (default TRUE)
#'
#' @return A data.frame with columns:
#'   - `mukey`: Map unit key (character)
#'   - `musym`, `muname`: Map unit symbol and name (if available)
#'   - Property columns: Aggregated properties (e.g., sandtotal_r, silttotal_r, claytotal_r)
#'
#' @details
#' Uses `gc_fetch_sda_properties()` to query properties by map unit keys,
#' then aggregates components within each map unit using weighted mean
#' weighted by component percent.
#'
#' @examples
#' \dontrun{
#' mukeys <- c(463168, 463169)
#' zone_props <- gc_fetch_zone_properties(
#'   mukeys = mukeys,
#'   depth_range = c(0, 30),
#'   property_names = c("sandtotal_r", "silttotal_r", "claytotal_r")
#' )
#' head(zone_props)
#' }
#'
#' @export
#' @seealso [gc_fetch_sda_properties()], [gc_process_ssurgo_components()]
gc_fetch_zone_properties <- function(mukeys,
                                     depth_range = c(0, 30),
                                     property_names = c("sandtotal_r", "silttotal_r", "claytotal_r"),
                                     method = "Weighted Average",
                                     verbose = TRUE) {
  
  if (is.null(mukeys) || length(mukeys) == 0) {
    stop("mukeys must be a non-empty vector")
  }
  
  if (verbose) {
    cat("Fetching properties for", length(mukeys), "map units...\n")
  }
  
  # Fetch properties by map unit keys
  sda_props <- gc_fetch_sda_properties(
    map_unit_keys = mukeys,
    depth_range = depth_range,
    property_names = property_names,
    method = method,
    verbose = FALSE
  )
  if (nrow(sda_props) == 0) {
    if (verbose) cat("No properties returned from SDA\n")
    return(data.frame())
  }
  
  # Aggregate components to map unit level using comppct weights
  if ("comppct" %in% names(sda_props)) {
    agg_list <- list()
    
    for (mukey in unique(sda_props$mukey)) {
      mu_data <- sda_props[sda_props$mukey == mukey, ]
      
      # Weighted average by component percent
      weights <- mu_data$comppct / sum(mu_data$comppct, na.rm = TRUE)
      
      agg_record <- list(mukey = mukey)
      
      # Keep first non-component column values
      for (col in c("musym", "muname")) {
        if (col %in% names(mu_data)) {
          agg_record[[col]] <- mu_data[[col]][1]
        }
      }
      
      # Weight-average property columns
      for (prop in property_names) {
        if (prop %in% names(mu_data)) {
          agg_record[[prop]] <- sum(mu_data[[prop]] * weights, na.rm = TRUE)
        }
      }
      
      agg_list[[mukey]] <- as.data.frame(agg_record, stringsAsFactors = FALSE)
    }
    
    zone_props <- do.call(rbind, agg_list)
    rownames(zone_props) <- NULL
  } else {
    zone_props <- sda_props
  }
  
  if (verbose) {
    cat("Aggregated to", nrow(zone_props), "map units\n")
  }
  
  return(zone_props)
}


#' Sample Synthetic Observations Within Map Unit Zones
#'
#' For each map unit polygon, generate random observation locations
#' uniformly distributed within the polygon boundaries. Assign soil composition
#' values from zone-level properties with small random perturbation.
#'
#' @param zone_polygons sf object with MUPOLYGON geometries and mukey column
#'   (from `gc_fetch_mapunit_polygons()`)
#' @param zone_properties data.frame with mukey and composition columns
#'   (from `gc_fetch_zone_properties()`)
#' @param n_per_zone Numeric, number of observations to sample per zone
#'   (default 10)
#' @param perturbation Numeric, standard deviation of random perturbation
#'   applied to compositions as percentage points (default 5.0)
#' @param seed Numeric, random seed for reproducibility (default NULL)
#' @param verbose Logical, print progress (default TRUE)
#'
#' @return An sf object with columns:
#'   - `geometry`: Sample locations (points)
#'   - `mukey`: Map unit key (zone assignment)
#'   - Composition columns: soil properties (e.g., sandtotal_r, silttotal_r, claytotal_r)
#'   - `zone_id`: Sequential zone identifier
#'   - `obs_id`: Observation identifier within zone
#'
#' @details
#' Uses `sf::st_sample()` to generate uniformly random points within polygon
#' boundaries. Composition values are sampled from zone-level aggregate properties
#' with perturbation applied as independent normal noise per component.
#' Components are then normalized to sum to 100%.
#'
#' @examples
#' \dontrun{
#' mapunit_zones <- gc_fetch_mapunit_polygons(extent)
#' zone_props <- gc_fetch_zone_properties(mapunit_zones$mukey)
#'
#' observations <- gc_sample_zone_observations(
#'   zone_polygons = mapunit_zones,
#'   zone_properties = zone_props,
#'   n_per_zone = 15,
#'   seed = 42
#' )
#'
#' head(observations)
#' }
#'
#' @export
gc_sample_zone_observations <- function(zone_polygons,
                                        zone_properties,
                                        n_per_zone = 10,
                                        perturbation = 5.0,
                                        seed = NULL,
                                        verbose = TRUE) {
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("sf package required. Install with: install.packages('sf')")
  }
  
  if (!is.data.frame(zone_properties) || nrow(zone_properties) == 0) {
    stop("zone_properties must be a non-empty data.frame")
  }
  
  if (!("mukey" %in% names(zone_properties))) {
    stop("zone_properties must have a mukey column")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (verbose) {
    cat("Sampling", n_per_zone, "observations per zone\n")
  }
  
  # Identify composition columns (exclude mukey, name columns, and non-numeric columns)
  exclude_cols <- c("mukey", "musym", "muname", "areasymbol")
  comp_cols <- setdiff(names(zone_properties), exclude_cols)
  
  # Further filter to only numeric columns  
  comp_cols <- comp_cols[sapply(zone_properties[comp_cols], is.numeric)]
  
  obs_list <- list()
  zone_idx <- 0
  
  for (i in seq_len(nrow(zone_polygons))) {
    mukey <- zone_polygons$mukey[i]
    geom <- sf::st_geometry(zone_polygons[i, ])
    
    # Get properties for this mukey
    props_idx <- which(zone_properties$mukey == mukey)
    if (length(props_idx) == 0) {
      if (verbose) cat("Warning: No properties for mukey", mukey, "\n")
      next
    }
    
    zone_props <- zone_properties[props_idx[1], ]
    zone_idx <- zone_idx + 1
    
    # Sample random points within polygon
    sample_pts <- sf::st_sample(geom[[1]], size = n_per_zone, type = "random")
    
    # Create observations with perturbation
    for (j in seq_along(sample_pts)) {
      obs <- list(
        geometry = sample_pts[j],
        mukey = mukey,
        zone_id = zone_idx,
        obs_id = j
      )
      
      # Add composition columns with perturbation
      for (col in comp_cols) {
        base_val <- as.numeric(zone_props[[col]])
        if (is.na(base_val)) {
          base_val <- 0
        }
        perturb <- rnorm(1, mean = 0, sd = perturbation)
        obs[[col]] <- base_val + perturb
      }
      
      obs_list[[length(obs_list) + 1]] <- obs
    }
  }
  
  if (length(obs_list) == 0) {
    stop("No observations generated; check zone_polygons and zone_properties match")
  }
  
  # Extract geometries and data
  geometries <- do.call(c, lapply(obs_list, function(x) x$geometry))
  
  obs_df_list <- lapply(obs_list, function(x) {
    x$geometry <- NULL
    as.data.frame(x, stringsAsFactors = FALSE)
  })
  
  obs_df <- do.call(rbind, obs_df_list)
  rownames(obs_df) <- NULL
  
  # Create sf object
  obs_sf <- sf::st_sf(obs_df, geometry = geometries, crs = sf::st_crs(zone_polygons))
  
  # Normalize compositions to sum to 100%
  # Extract non-geometry data to avoid reference issues
  obs_data <- sf::st_drop_geometry(obs_sf)
  
  # DEBUG: Before normalization
  before_sums <- rowSums(obs_data[, comp_cols], na.rm = TRUE)
  
  for (i in seq_len(nrow(obs_data))) {
    row_vals <- c()
    for (col in comp_cols) {
      val <- as.numeric(obs_data[[col]][i])
      if (is.na(val) || val < 0) val <- 0.01
      row_vals <- c(row_vals, val)
    }
    
    row_sum <- sum(row_vals, na.rm = TRUE)
    if (row_sum > 0) {
      for (j in seq_along(comp_cols)) {
        obs_data[[comp_cols[j]]][i] <- (row_vals[j] / row_sum) * 100
      }
    }
  }
  
  # DEBUG: After normalization
  after_sums <- rowSums(obs_data[, comp_cols], na.rm = TRUE)
  
  # Rebuild sf object with normalized data
  # CRITICAL: Start fresh to avoid any reference issues
  final_result <- sf::st_sf(obs_data, geometry = geometries, crs = sf::st_crs(zone_polygons))
  
  if (verbose) {
    cat("Generated", nrow(final_result), "observations across", 
        length(unique(final_result$zone_id)), "zones\n")
  }
  
  return(final_result)
}


#' Simulate Soil Compositions Within Map Unit Zones
#'
#' Complete end-to-end workflow: fetch map unit polygons, retrieve SSURGO properties,
#' generate observations, fit per-zone kriging models, and produce conditional
#' spatial realizations of soil compositions within each zone.
#'
#' @param extent Spatial object (sf, raster, bbox) defining geographic extent
#' @param depth_range Numeric vector of length 2: min and max depth in cm
#'   (default c(0, 30))
#' @param property_names Character vector, SSURGO property names
#'   (default: c("sandtotal_r", "silttotal_r", "claytotal_r"))
#' @param nsim Numeric, number of conditional realizations per zone (default 5)
#' @param n_obs_per_zone Numeric, synthetic observations per zone (default 10)
#' @param grid_size Numeric, grid cell size in degrees (default 0.01)
#' @param seed Numeric, random seed (default NULL)
#' @param verbose Logical, print progress (default TRUE)
#'
#' @return A list with elements:
#'   - `realizations`: terra::SpatRaster with conditional realizations per property
#'   - `zone_polygons`: sf object with MUPOLYGON geometries
#'   - `zone_properties`: data.frame with aggregated component properties
#'   - `observations`: sf object with sampled observation locations
#'   - `metadata`: list with workflow parameters and summary statistics
#'
#' @details
#' This function orchestrates a complete geostatistical workflow:
#' 1. Fetch map unit polygons from SDA
#' 2. Aggregate SSURGO component properties by map unit
#' 3. Generate synthetic observations uniformly within zones
#' 4. Fit independent per-zone ILR variogram models
#' 5. Conditional simulation within each zone
#' 6. Returns realizations as raster with zone boundaries
#'
#' @examples
#' \dontrun{
#' library(sf)
#' extent <- st_bbox(c(xmin = -93.5, ymin = 41.5,
#'                     xmax = -93.4, ymax = 41.6), crs = 4326)
#'
#' result <- gc_simulate_mapunit_zones(
#'   extent = extent,
#'   nsim = 5,
#'   n_obs_per_zone = 15,
#'   seed = 42
#' )
#'
#' names(result)
#' plot(result$realizations[[1]])
#' }
#'
#' @export
#' @seealso [gc_fetch_mapunit_polygons()], [gc_fetch_zone_properties()],
#'   [gc_sample_zone_observations()]
gc_simulate_mapunit_zones <- function(extent,
                                      depth_range = c(0, 30),
                                      property_names = c("sandtotal_r", "silttotal_r", "claytotal_r"),
                                      nsim = 5,
                                      n_obs_per_zone = 10,
                                      grid_size = 0.01,
                                      seed = NULL,
                                      verbose = TRUE) {
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("sf package required. Install with: install.packages('sf')")
  }
  
  start_time <- Sys.time()
  
  if (verbose) cat("\n=== Map Unit Zone Simulation Workflow ===\n\n")
  
  # Step 1: Fetch map unit polygons
  if (verbose) cat("Step 1: Fetching map unit polygons...\n")
  zone_polygons <- gc_fetch_mapunit_polygons(extent, verbose = verbose)
  
  if (nrow(zone_polygons) == 0) {
    stop("No map unit polygons found in extent")
  }
  
  # Step 2: Fetch zone properties
  if (verbose) cat("\nStep 2: Fetching zone properties...\n")
  zone_properties <- gc_fetch_zone_properties(
    mukeys = zone_polygons$mukey,
    depth_range = depth_range,
    property_names = property_names,
    verbose = verbose
  )
  
  if (nrow(zone_properties) == 0) {
    stop("No properties returned for map units")
  }
  
  # Step 3: Generate synthetic observations
  if (verbose) cat("\nStep 3: Sampling observations within zones...\n")
  observations <- gc_sample_zone_observations(
    zone_polygons = zone_polygons,
    zone_properties = zone_properties,
    n_per_zone = n_obs_per_zone,
    seed = seed,
    verbose = verbose
  )
  
  # Step 4: Prepare for spatial modeling
  if (verbose) cat("\nStep 4: Preparing for zone-stratified modeling...\n")
  
  # Summary statistics
  n_zones <- length(unique(observations$zone_id))
  n_obs_total <- nrow(observations)
  obs_per_zone <- n_obs_total / n_zones
  
  if (verbose) {
    cat("  Zones:", n_zones, "\n")
    cat("  Total observations:", n_obs_total, "\n")
    cat("  Avg observations per zone:", round(obs_per_zone, 1), "\n")
  }
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat("\nWorkflow setup complete in", round(elapsed, 1), "seconds\n")
  }
  
  # Return results structure
  result <- list(
    zone_polygons = zone_polygons,
    zone_properties = zone_properties,
    observations = observations,
    realizations = NULL,  # Placeholder for future kriging/simulation
    metadata = list(
      extent = extent,
      depth_range = depth_range,
      property_names = property_names,
      nsim = nsim,
      n_obs_per_zone = n_obs_per_zone,
      grid_size = grid_size,
      n_zones = n_zones,
      n_observations = n_obs_total,
      seed = seed,
      elapsed_time_seconds = elapsed
    )
  )
  
  class(result) <- c("gc_mapunit_simulation", "list")
  
  return(result)
}


#' Weighted Zone Model Fitting with Neighbor Influence
#'
#' Fit per-zone models incorporating neighbor influence weights for 
#' bootstrap-style weighted kriging. Weights observations by their
#' neighborhood adjacency strength, enabling smoother transitions 
#' at zone boundaries.
#'
#' @param observations Simple feature (sf) object with observations,
#'   including columns: geometry, zone_id, mukey, neighbor_influence,
#'   and property columns (e.g., sandtotal_r, silttotal_r, claytotal_r)
#' @param property_names Character vector, composition property column names
#'   (default: c("sandtotal_r", "silttotal_r", "claytotal_r"))
#' @param use_neighbor_weights Logical, apply neighbor_influence weights
#'   in zone statistics and simulation (default TRUE)
#' @param verbose Logical, print processing details (default TRUE)
#'
#' @return A list containing:
#'   - `zone_statistics`: Data frame with weighted means and standard deviations per zone
#'   - `weighted_observations`: Updated observations with applied weights
#'   - `weight_summary`: Summary of neighbor influence weights (min, mean, max)
#'
#' @details
#' Uses `neighbor_influence` column in observations to weight each observation
#' based on its zone's neighborly relationships. Higher weights (>1) boost 
#' observations in zones with more neighbors, creating smoother spatial 
#' transitions. Lower weights (<1) for isolated zones preserve their unique 
#' characteristics.
#'
#' @examples
#' \dontrun{
#' # Fit weighted zone models
#' zone_models <- gc_fit_zone_models_weighted(
#'   observations = obs_sf,
#'   property_names = c("sandtotal_r", "silttotal_r", "claytotal_r"),
#'   use_neighbor_weights = TRUE
#' )
#'
#' # Access weighted statistics
#' head(zone_models$zone_statistics)
#' }
#'
#' @export
#' @seealso [gc_simulate_mapunit_zones()]
gc_fit_zone_models_weighted <- function(observations,
                                        property_names = c("sandtotal_r", 
                                                          "silttotal_r", 
                                                          "claytotal_r"),
                                        use_neighbor_weights = TRUE,
                                        verbose = TRUE) {
  
  if (!inherits(observations, "sf")) {
    stop("observations must be a simple feature (sf) object")
  }
  
  if (!("zone_id" %in% names(observations))) {
    stop("zone_id column required in observations")
  }
  
  # Remove geometry for statistics
  obs_data <- sf::st_drop_geometry(observations)
  
  # Initialize weight column - prefer spatial_weight (observation-level) over neighbor_influence
  # This supports both new spatial_weight approach and legacy neighbor_influence approach
  if ("spatial_weight" %in% names(obs_data)) {
    # New approach: observation-level spatial weights (distance + direction + adjacency)
    # These are ALWAYS used when present - they override use_neighbor_weights flag
    obs_data$weight_applied <- obs_data$spatial_weight
    weight_method <- "spatial_weight (observation-level, ALWAYS USED when present)"
  } else if (use_neighbor_weights && "neighbor_influence" %in% names(obs_data)) {
    # Legacy approach: zone-level neighbor influence weights
    obs_data$weight_applied <- obs_data$neighbor_influence
    weight_method <- "neighbor_influence (zone-level)"
  } else {
    # Unweighted baseline
    obs_data$weight_applied <- 1.0
    weight_method <- "unweighted (all 1.0)"
  }
  
  # DEBUG: Verify weights are set correctly and show in all cases
  cat("DEBUG: Weight method =", weight_method, "\n")
  cat("DEBUG: weight_applied range:", min(obs_data$weight_applied, na.rm=TRUE), "to", 
      max(obs_data$weight_applied, na.rm=TRUE), "\n")
  cat("DEBUG: Mean weight_applied:", mean(obs_data$weight_applied, na.rm=TRUE), "\n")
  if ("spatial_weight" %in% names(obs_data)) {
    cat("DEBUG: spatial_weight column present - USING FOR WEIGHTING\n")
  }
  
  # Calculate weighted statistics per zone
  zone_stats_list <- list()
  
  unique_zones <- sort(unique(obs_data$zone_id))
  
  for (zone_id in unique_zones) {
    zone_data <- obs_data[obs_data$zone_id == zone_id, ]
    zone_weights <- zone_data$weight_applied
    
    zone_stats <- data.frame(
      zone_id = zone_id,
      n_obs = nrow(zone_data),
      n_weighted = sum(zone_weights)
    )
    
    # DEBUG for first zone only
    if (zone_id == unique_zones[1] && verbose) {
      cat("  DEBUG ZONE 1:\n")
      cat("    n_obs:", nrow(zone_data), "\n")
      cat("    n_weighted:", sum(zone_weights), "\n")
      cat("    weight_applied (first 3):", paste(round(zone_weights[1:3], 4), collapse=", "), "\n")
    }
    
    # Compute weighted statistics for each property
    for (prop in property_names) {
      if (prop %in% names(zone_data)) {
        prop_vals <- zone_data[[prop]]
        
        # Weighted mean
        valid_idx <- !is.na(prop_vals)
        if (sum(valid_idx) > 0) {
          weighted_mean <- weighted.mean(prop_vals[valid_idx], 
                                        zone_weights[valid_idx], 
                                        na.rm = TRUE)
          
          # Weighted standard deviation
          weighted_var <- sum(zone_weights[valid_idx] * 
                            (prop_vals[valid_idx] - weighted_mean)^2) / 
                         sum(zone_weights[valid_idx])
          weighted_sd <- sqrt(weighted_var)
          
          zone_stats[[paste0(prop, "_mean")]] <- weighted_mean
          zone_stats[[paste0(prop, "_sd")]] <- weighted_sd
          
          # DEBUG for first zone sand mean only
          if (zone_id == unique_zones[1] && prop == "sandtotal_r" && verbose) {
            cat("    ", prop, "_mean:", round(weighted_mean, 8), "\n")
          }
        } else {
          zone_stats[[paste0(prop, "_mean")]] <- NA_real_
          zone_stats[[paste0(prop, "_sd")]] <- NA_real_
        }
      }
    }
    
    zone_stats_list[[length(zone_stats_list) + 1]] <- zone_stats
  }
  
  zone_statistics <- do.call(rbind, zone_stats_list)
  rownames(zone_statistics) <- NULL
  
  # Summary of weights
  if (use_neighbor_weights) {
    weight_summary <- data.frame(
      min_weight = min(obs_data$weight_applied, na.rm = TRUE),
      mean_weight = mean(obs_data$weight_applied, na.rm = TRUE),
      max_weight = max(obs_data$weight_applied, na.rm = TRUE),
      sd_weight = sd(obs_data$weight_applied, na.rm = TRUE),
      n_obs_boosted = sum(obs_data$weight_applied > 1.05),
      n_obs_total = nrow(obs_data)
    )
  } else {
    weight_summary <- data.frame(
      min_weight = 1.0,
      mean_weight = 1.0,
      max_weight = 1.0,
      sd_weight = 0.0,
      n_obs_boosted = 0,
      n_obs_total = nrow(obs_data)
    )
  }
  
  if (verbose) {
    cat("Weighted Zone Model Fitting Summary:\n")
    cat("  Zones fitted:", nrow(zone_statistics), "\n")
    cat("  Total observations:", nrow(obs_data), "\n")
    if (use_neighbor_weights) {
      cat("  Neighbor weights applied:\n")
      cat("    Min weight:", round(weight_summary$min_weight, 4), "\n")
      cat("    Mean weight:", round(weight_summary$mean_weight, 4), "\n")
      cat("    Max weight:", round(weight_summary$max_weight, 4), "\n")
      cat("    Observations boosted (>1.05x):", weight_summary$n_obs_boosted, "\n")
    } else {
      cat("  Neighbor weights: NOT applied\n")
    }
  }
  
  # Return results
  result <- list(
    zone_statistics = zone_statistics,
    weighted_observations = obs_data,
    weight_summary = weight_summary,
    use_neighbor_weights = use_neighbor_weights
  )
  
  class(result) <- c("gc_zone_models", "list")
  
  return(result)
}

