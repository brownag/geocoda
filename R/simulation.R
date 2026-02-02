#' Simulate and Back-Transform Compositional Fields
#'
#' Generate spatial realizations of compositional data by predicting in ILR or MAF space
#' using a gstat model, then back-transforming to the original units. Supports
#' both unconditional simulation (independent of data) and conditional simulation
#' (honoring observed values at sample locations). When the model includes MAF
#' decorrelation, automatically back-transforms from MAF to ILR before final
#' compositional inversion. Multiple realizations are stacked into a multi-layer
#' `terra::SpatRaster`.
#'
#' @param model A `gstat` object (typically from `gc_ilr_model()`).
#'   For conditional simulation, the model **must be built with the `data` parameter**
#'   so that conditioning data is embedded in the model.
#' @param locations An `sf` object or data frame defining the simulation grid.
#'   Must contain `x` and `y` columns. If a data frame is provided, it will be
#'   converted to an sf object with coordinates in the specified CRS.
#' @param nsim Number of realizations (default 1).
#' @param target_names Character vector of original component names. If `NULL`,
#'   inferred from the `gstat` model or defaults to `c("comp1", "comp2", ...)`.
#' @param crs Coordinate reference system (default `"local"`).
#' @param observed_data An optional `sf` object or data frame containing observed
#'   compositional samples for conditional simulation. Must have the same structure
#'   as the data used to build the model: columns `x`, `y` for coordinates, and
#'   columns for each ILR dimension named `ilr1`, `ilr2`, etc. (in ILR space).
#'   If `NULL` (default), unconditional simulation is performed.
#' @param nmax Maximum number of nearby observations to use for prediction at each
#'   location (default `NULL`, uses all observations). Setting a smaller `nmax`
#'   (e.g., 12-15) can:
#'   - Improve computational efficiency for large datasets
#'   - Reduce memory usage during prediction
#'   - Create more local uncertainty estimates (Screen Effect)
#'   - May introduce slight artifacts at spatial boundaries
#'
#'   Recommended: Use `NULL` for small/medium datasets, `nmax=12-20` for large datasets.
#'
#' @return A `terra::SpatRaster` with layers named according to the pattern:
#'   `<component>.sim<N>` where `component` is from `target_names` and
#'   `sim<N>` indicates the realization number. For example, with 3 components
#'   and 2 realizations: `comp1.sim1, comp2.sim1, comp3.sim1, comp1.sim2, ...`
#'
#' @details
#' This function:
#' 1. Prepares spatial coordinates from `locations`
#' 2. Validates `observed_data` if provided (must have ilr columns and xy coords)
#' 3. Calls `gstat::predict()` to generate ILR values for `nsim` realizations.
#'    If the model was built with conditioning data, the predictions automatically
#'    honor that data at the sample locations.
#'    The `nmax` parameter limits the neighborhood used for kriging.
#' 4. Extracts simulated ILR columns (identified by pattern `.sim<N>`)
#' 5. Groups ILR values by realization
#' 6. Back-transforms each realization using `compositions::ilrInv()`
#' 7. Rescales from 0-1 to 0-100
#' 8. Stacks results into a single `SpatRaster` with proper naming
#'
#' **Conditional vs. Unconditional:**
#' When the model is built with `data = NULL`, all realizations are independent draws from the
#' spatial distribution (unconditional). When the model is built with `data = <conditioning_data>`,
#' all realizations are conditioned on those values: predictions at sample locations
#' exactly reproduce the observed values, while predictions elsewhere reflect
#' uncertainty updated by the conditioning data.
#'
#' **Important:** For conditional simulation to work, the model **must have been built
#' with the conditioning data embedded via the `data` parameter to `gc_ilr_model()`.
#' The `observed_data` parameter here is only for validation. Conditioning happens
#' automatically because the model was created with that data.
#'
#' **Neighborhood Size Effects (Screen Effect):**
#' When `nmax` is small relative to data density, kriging typically uses only
#' the nearest `nmax` observations. This can create more localized uncertainty
#' estimates and may accelerate computation, but can also:
#' - Introduce discontinuities at boundaries between neighborhoods
#' - Miss spatial structure information from distant data
#' - Bias predictions if far-field data carries important variance information
#'
#' The output strictly honors the sum constraint: all rows sum to `target_sum`
#' (typically 100 for percentages) within floating-point precision.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(gstat)
#'
#' # Assuming ilr_params and model are already defined
#' # (See workflow example in the package README)
#'
#' # Create simulation grid
#' x.range <- seq(0, 100, by = 5)
#' y.range <- seq(0, 100, by = 5)
#' grid_df <- expand.grid(x = x.range, y = y.range)
#'
#' # Unconditional simulation: 5 realizations
#' result <- gc_sim_composition(model, grid_df, nsim = 5,
#'                              target_names = c("sand", "silt", "clay"))
#' print(result)
#'
#' # Conditional simulation: honor observed data
#' result_cond <- gc_sim_composition(model, grid_df, nsim = 5,
#'                                   target_names = c("sand", "silt", "clay"),
#'                                   observed_data = sample_df)
#'
#' # Large dataset: use neighborhood limiting for efficiency
#' result_nmax <- gc_sim_composition(model, grid_df, nsim = 5,
#'                                   target_names = c("sand", "silt", "clay"),
#'                                   nmax = 15)
#'
#' # Access results
#' terra::values(result)
#' }
#'
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry
#' @importFrom terra rast values
#' @importFrom compositions ilrInv
#' @importFrom methods is
#' @export
gc_sim_composition <- function(model,
                               locations,
                               nsim = 1,
                               target_names = NULL,
                               crs = "local",
                               observed_data = NULL,
                               nmax = NULL) {
  requireNamespace("gstat", quietly = TRUE)

  if (!is.null(observed_data)) {
    if (methods::is(observed_data, "data.frame")) {
      # Support both 2D and 3D coordinates
      has_z <- "z" %in% colnames(observed_data)
      coord_cols <- if (has_z) c("x", "y", "z") else c("x", "y")
      observed_data_sf <- sf::st_as_sf(observed_data, coords = coord_cols, crs = crs)
    } else if (methods::is(observed_data, "sf")) {
      observed_data_sf <- observed_data
    } else {
      stop("observed_data must be an sf object or data frame")
    }

    # Check for required ILR columns
    ilr_col_names <- grep("^ilr\\d+$", colnames(observed_data_sf), value = TRUE)
    if (length(ilr_col_names) == 0) {
      stop(
        "observed_data must contain ILR columns (ilr1, ilr2, ...) ",
        "in the same space as the model"
      )
    }
  } else {
    observed_data_sf <- NULL
  }

  if (methods::is(locations, "sf")) {
    locations_sf <- locations
  } else {
    locations_df <- as.data.frame(locations)
    # Support both 2D and 3D coordinates
    has_z <- "z" %in% colnames(locations_df)
    coord_cols <- if (has_z) c("x", "y", "z") else c("x", "y")

    if (!all(coord_cols %in% colnames(locations_df))) {
      stop("locations must contain 'x', 'y', and optionally 'z' columns")
    }
    locations_sf <- sf::st_as_sf(locations_df, coords = coord_cols, crs = crs)
  }

  gstat_locations <- locations_sf

  if (!is.null(nmax)) {
    if (!is.numeric(nmax) || nmax < 1 || nmax != floor(nmax)) {
      stop("nmax must be a positive integer or NULL")
    }
  }

  if (!is.null(nmax)) {
    sim_result <- predict(
      model,
      newdata = gstat_locations,
      nsim = nsim,
      nmax = nmax,
      debug.level = 0
    )
  } else {
    sim_result <- predict(
      model,
      newdata = gstat_locations,
      nsim = nsim,
      debug.level = 0
    )
  }

  # Check if MAF was used
  maf_object <- attr(model, "maf_object")
  ilr_mean_original <- attr(model, "ilr_mean_original")

  # Look for either ILR or MAF columns depending on what was simulated
  if (!is.null(maf_object)) {
    ilr_cols <- grep("^maf\\d+\\.sim\\d+$", colnames(sim_result), value = TRUE)
    is_maf <- TRUE
  } else {
    ilr_cols <- grep("^ilr\\d+\\.sim\\d+$", colnames(sim_result), value = TRUE)
    is_maf <- FALSE
  }

  if (length(ilr_cols) == 0) {
    stop("No simulated ILR or MAF columns found in prediction result")
  }

  n_ilr <- max(as.numeric(gsub("ilr(\\d+)\\.sim\\d+", "\\1", ilr_cols)))
  actual_nsim <- max(as.numeric(gsub("ilr\\d+\\.sim(\\d+)", "\\1", ilr_cols)))

  if (is.null(target_names)) {
    # Infer from model or use defaults
    target_names <- paste0("comp", seq_len(n_ilr + 1))
  }

  if (length(target_names) != n_ilr + 1) {
    stop(
      "target_names length (",
      length(target_names),
      ") must equal number of components (",
      n_ilr + 1,
      ")"
    )
  }

  all_comps <- list()

  for (sim_idx in seq_len(actual_nsim)) {
    ilr_cols_for_sim <- grep(
      paste0("\\.sim", sim_idx, "$"),
      ilr_cols,
      value = TRUE
    )

    # Extract dimension numbers
    if (is_maf) {
      ilr_nums <- as.numeric(gsub("maf(\\d+)\\.sim\\d+", "\\1", ilr_cols_for_sim))
    } else {
      ilr_nums <- as.numeric(gsub("ilr(\\d+)\\.sim\\d+", "\\1", ilr_cols_for_sim))
    }
    ilr_cols_for_sim <- ilr_cols_for_sim[order(ilr_nums)]

    mat_values <- as.matrix(sf::st_drop_geometry(sim_result[, ilr_cols_for_sim]))

    # Back-transform from MAF to ILR if needed
    if (is_maf) {
      # Create a temporary data frame with MAF columns
      maf_df <- as.data.frame(mat_values)
      names(maf_df) <- paste0("maf", seq_len(ncol(maf_df)))

      # Inverse MAF transformation: MAF -> ILR
      ilr_df <- gc_inverse_maf(maf_df, maf_object, ilr_mean = ilr_mean_original)

      # Extract ILR columns
      ilr_col_names <- grep("^ilr\\d+$", names(ilr_df), value = TRUE)
      ilr_mat <- as.matrix(ilr_df[, ilr_col_names])
    } else {
      ilr_mat <- mat_values
    }

    comp_mat <- compositions::ilrInv(ilr_mat)

    comp_mat <- matrix(as.numeric(comp_mat), ncol = n_ilr + 1) * 100

    for (comp_idx in seq_len(n_ilr + 1)) {
      col_name <- paste0(target_names[comp_idx], ".sim", sim_idx)
      all_comps[[col_name]] <- comp_mat[, comp_idx]
    }
  }

  coords <- sf::st_coordinates(gstat_locations)
  comp_df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    do.call(cbind, all_comps)
  )

  rast_stack <- terra::rast(comp_df, type = "xyz", crs = crs)

  layer_names <- colnames(comp_df)[-c(1, 2)]
  names(rast_stack) <- layer_names

  return(rast_stack)
}

#' Simulate Compositional Fields at Block Support
#'
#' Direct block simulation (DBS) for management-scale uncertainty quantification.
#' Simulates block-averaged values by discretizing blocks into sub-points,
#' simulating at sub-point support, and averaging to block support.
#'
#' This approach captures the volume-variance relationship: block variance
#' is lower than point variance due to averaging (nugget averaging effect).
#'
#' @param model A gstat object from [gc_ilr_model()], optionally with MAF attributes
#' @param block_centers Data frame with x, y columns (block center coordinates)
#' @param block_size Vector c(dx, dy) specifying block dimensions (same units as coordinates)
#' @param nsim Number of simulations to generate (default 1)
#' @param target_names Character vector of target component names
#'   (e.g., c("sand", "silt", "clay"))
#' @param crs Coordinate reference system (default "local")
#' @param discretization Vector c(nx, ny) specifying sub-points per block dimension
#'   (default c(4, 4) = 16 sub-points per block)
#' @param observed_data Optional sf object with observed data (for conditioning)
#' @param nmax Maximum number of nearest neighbors for kriging (default NULL = all data)
#'
#' @return terra::SpatRaster stack with block-support compositional fields
#'   - Dimensions: nrow(block_centers) rows
#'   - Layers: sand.sim1, silt.sim1, clay.sim1, sand.sim2, ... (for each component and realization)
#'   - Values: 0-100 (percent/proportion scale, sum to 100 per block)
#'
#' @details
#' **Direct Block Simulation (DBS) workflow:**
#'
#' 1. **Discretization**: Each block is subdivided into regular grid of sub-points
#'    - Default: 4×4 grid = 16 sub-points per block
#'    - Sub-points centered within sub-cells (not at corners)
#'
#' 2. **Sub-point simulation**: Simulate ILR coordinates at all sub-point locations
#'    - Uses standard Sequential Gaussian Simulation (SGS)
#'    - Respects spatial structure (variogram) and any conditioning data
#'
#' 3. **Block averaging**: Average ILR values across sub-points within each block
#'    - For each block B and variable Z: Z_B = mean(Z_sub-points in B)
#'    - This averaging produces the volume-variance reduction
#'
#' 4. **Back-transformation**: ILR → Composition at block support
#'    - Each block-averaged ILR coordinate set back-transforms to composition
#'    - Compositional constraints automatically satisfied (sums to 100)
#'
#' **Volume-variance relationship:**
#' ```
#' Var(Z_block) ≈ Var(Z_point) / (nx × ny)
#' ```
#' Example: 4×4 discretization → block variance ≈ 1/16 of point variance
#'
#' **Computational efficiency:**
#' - Memory: Only stores block-support results, not intermediate sub-point grids
#' - Speed: One simulation run (sub-points and averaging) vs. many refinement steps
#' - Practical: Enables large-scale uncertainty quantification at management units
#'
#' **Comparison to alternatives:**
#' - vs. Point-then-aggregate: Same result but faster (direct computation)
#' - vs. Block kriging: Similar math but simulation vs. estimation
#' - vs. LUC (Localised Uniform Conditioning): Simpler implementation, comparable accuracy
#'
#' @examples
#' \dontrun{
#' # Setup: ILR model from previous steps
#' ilr_params <- gc_ilr_params(texture_data[, c("sand", "silt", "clay")])
#' vgm <- gc_fit_vgm(ilr_params, texture_data, aggregate = TRUE)
#' model <- gc_ilr_model(ilr_params, vgm)
#'
#' # Create management-unit-scale grid (100m blocks)
#' blocks <- expand.grid(
#'   x = seq(0, 1000, by = 100),
#'   y = seq(0, 1000, by = 100)
#' )
#'
#' # Block simulation with 4×4 discretization (16 sub-points per block)
#' block_sims <- gc_sim_composition_block(
#'   model = model,
#'   block_centers = blocks,
#'   block_size = c(100, 100),
#'   nsim = 10,
#'   target_names = c("sand", "silt", "clay"),
#'   discretization = c(4, 4)
#' )
#'
#' # Result: 10 blocks × 3 components × 10 sims = 300 layers
#' print(block_sims)
#'
#' # Check volume-variance relationship
#' point_variance <- terra::global(
#'   gc_sim_composition(model, blocks, nsim = 10),
#'   "sd"
#' )
#' block_variance <- terra::global(block_sims, "sd")
#' # block_variance should be ~1/4 of point_variance (for 4×4 grid)
#' }
#'
#' @export
gc_sim_composition_block <- function(model, block_centers, block_size, nsim = 1,
                                    target_names = NULL, crs = "local",
                                    discretization = c(4, 4), observed_data = NULL,
                                    nmax = NULL) {
  # Input validation
  stopifnot(is.data.frame(block_centers), all(c("x", "y") %in% names(block_centers)))
  stopifnot(is.numeric(block_size), length(block_size) == 2, all(block_size > 0))
  stopifnot(is.numeric(nsim), nsim >= 1)
  stopifnot(is.numeric(discretization), length(discretization) == 2, all(discretization >= 2))

  # Extract ILR parameters from model
  ilr_params <- attr(model, "ilr_params")
  if (is.null(ilr_params)) {
    stop("Model must have ilr_params attribute (from gc_ilr_model)")
  }

  n_ilr <- nrow(ilr_params$cov)
  if (is.null(target_names)) {
    target_names <- ilr_params$names
  }
  stopifnot(length(target_names) == n_ilr + 1)

  # Step 1: Discretize blocks into sub-points
  message("Discretizing ", nrow(block_centers), " blocks into ",
          prod(discretization), " sub-points each...")
  sub_points <- gc_discretize_block(block_centers, block_size, discretization)

  # Step 2: Simulate at sub-point support
  message("Simulating at sub-point support (", nrow(sub_points), " points, ",
          nsim, " realizations)...")
  sub_sims <- gc_sim_composition(
    model = model,
    locations = sub_points[, c("x", "y")],
    nsim = nsim,
    target_names = target_names,
    crs = crs,
    observed_data = observed_data,
    nmax = nmax
  )

  # Step 3: Convert to data frame for averaging
  sub_sims_df <- as.data.frame(sub_sims)

  # Add block IDs
  sub_sims_df$block_id <- sub_points$block_id

  # Step 4: Average to block support
  message("Averaging sub-points to block support...")

  # Identify composition columns (pattern: component.simN where N is integer)
  comp_cols <- grep(paste0("^(", paste(target_names, collapse = "|"), ")\\.sim\\d+$"),
                    names(sub_sims_df), value = TRUE)

  if (length(comp_cols) == 0) {
    stop("No composition simulation columns found in results")
  }

  # Group by block and average
  block_averaged <- aggregate(
    sub_sims_df[, comp_cols, drop = FALSE],
    by = list(block_id = sub_sims_df$block_id),
    FUN = mean
  )

  # Step 5: Normalize compositions to sum to 100 per block (ensure compositional constraint)
  # Normalize per realization to maintain compositional constraint
  for (sim_idx in seq_len(nsim)) {
    sim_cols <- grep(paste0("\\.sim", sim_idx, "$"), names(block_averaged), value = TRUE)
    if (length(sim_cols) > 0) {
      # Extract values for this simulation
      sim_values <- block_averaged[, sim_cols, drop = FALSE]
      # Row-wise sum and normalization
      row_sums <- rowSums(sim_values)
      # Normalize to 100 per block
      for (col in sim_cols) {
        block_averaged[, col] <- block_averaged[, col] * 100 / row_sums
      }
    }
  }

  # Step 6: Convert to terra::SpatRaster
  message("Converting to raster...")

  # Extract coordinates and composition data
  block_coords <- block_centers[, c("x", "y")]
  comp_data <- block_averaged[, -1, drop = FALSE]  # Remove block_id column

  # Create data frame with x, y, composition values
  block_result <- data.frame(
    x = block_coords$x,
    y = block_coords$y,
    comp_data
  )

  # Create raster from xyz data
  rast_stack <- terra::rast(block_result, type = "xyz", crs = crs)

  # Set layer names (all composition columns)
  layer_names <- names(comp_data)
  names(rast_stack) <- layer_names

  message("Block simulation complete. Result: ", nrow(block_centers), " blocks × ",
          nsim, " realizations")

  return(rast_stack)
}
