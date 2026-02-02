#' Define Soil Horizon Surfaces for Stratigraphic Transformation
#'
#' Creates surface definitions for top and bottom of stratigraphic unit.
#' Surfaces can be constant (flat), DEM-based, or interpolated from point data.
#'
#' @param top_surface Numeric (constant elevation), terra::SpatRaster (DEM),
#'   or data.frame with x, y, z columns (points to interpolate)
#' @param bottom_surface Same format as top_surface
#' @param type Character: "constant", "raster", or "interpolate"
#' @param interpolation_method If type = "interpolate": "idw", "kriging", or "spline"
#' @param extent terra::ext object defining spatial extent (required if type = "interpolate")
#' @param resolution Vector c(dx, dy) for interpolated rasters (required if type = "interpolate")
#' @param crs Coordinate reference system (default "local")
#'
#' @return List with:
#'   - top: terra::SpatRaster of top surface elevations
#'   - bottom: terra::SpatRaster of bottom surface elevations
#'   - thickness: terra::SpatRaster of (top - bottom)
#'   - type: Transformation type used
#'   - crs: Coordinate reference system
#'
#' @details
#' **Surface Types**:
#' - constant: Flat surfaces (e.g., top = 0, bottom = -100 for 0-100 cm depth)
#' - raster: Existing DEM or interpolated horizon surfaces
#' - interpolate: Point observations of horizon depths, interpolate to raster
#'
#' **Validation Checks**:
#' - top > bottom everywhere (no inversions)
#' - No NA values in overlap region
#' - Thickness > 0 everywhere
#' - Consistent CRS if rasters provided
#'
#' @examples
#' \dontrun{
#' # Constant surfaces (0-100 cm depth)
#' surfaces <- gc_define_surfaces(
#'   top_surface = 0,
#'   bottom_surface = -100,
#'   type = "constant",
#'   extent = terra::ext(c(0, 1000, 0, 1000)),
#'   resolution = c(10, 10)
#' )
#'
#' # DEM-based (topsoil 0-30 cm)
#' dem <- terra::rast("elevation.tif")
#' surfaces <- gc_define_surfaces(
#'   top_surface = dem,
#'   bottom_surface = dem - 30,  # 30 cm below surface
#'   type = "raster"
#' )
#' }
#'
#' @export
gc_define_surfaces <- function(top_surface, bottom_surface, type = "constant",
                              interpolation_method = "idw",
                              extent = NULL, resolution = NULL, crs = "local") {
  # Input validation
  stopifnot(type %in% c("constant", "raster", "interpolate"))

  # Build top raster
  if (type == "constant") {
    stopifnot(is.numeric(top_surface), length(top_surface) == 1)
    stopifnot(!is.null(extent), !is.null(resolution))
    top_rast <- terra::rast(extent, res = resolution, crs = crs)
    terra::values(top_rast) <- top_surface
  } else if (type == "raster") {
    stopifnot(inherits(top_surface, "SpatRaster"))
    top_rast <- top_surface
  } else {  # interpolate
    stopifnot(is.data.frame(top_surface))
    stopifnot(all(c("x", "y", "z") %in% names(top_surface)))
    stopifnot(!is.null(extent), !is.null(resolution))
    top_rast <- gc_interpolate_surface(top_surface, interpolation_method, extent, resolution, crs)
  }

  # Build bottom raster (similar logic)
  if (type == "constant") {
    bottom_rast <- terra::rast(extent, res = resolution, crs = crs)
    terra::values(bottom_rast) <- bottom_surface
  } else if (type == "raster") {
    stopifnot(inherits(bottom_surface, "SpatRaster"))
    bottom_rast <- bottom_surface
  } else {
    stopifnot(is.data.frame(bottom_surface))
    stopifnot(all(c("x", "y", "z") %in% names(bottom_surface)))
    bottom_rast <- gc_interpolate_surface(bottom_surface, interpolation_method, extent, resolution, crs)
  }

  # Validation: top > bottom everywhere
  thickness_rast <- top_rast - bottom_rast
  invalid_mask <- terra::values(thickness_rast) <= 0
  if (any(invalid_mask, na.rm = TRUE)) {
    n_invalid <- sum(invalid_mask, na.rm = TRUE)
    stop("Invalid surfaces: top must be > bottom everywhere. Found ",
         n_invalid, " cells with inversions or zero thickness.")
  }

  if (any(is.na(terra::values(thickness_rast)))) {
    warning("NA values found in thickness raster. Check surface alignment.")
  }

  return(list(
    top = top_rast,
    bottom = bottom_rast,
    thickness = thickness_rast,
    type = type,
    crs = crs
  ))
}


#' Internal: Interpolate Surface from Points
#'
#' @param points Data frame with x, y, z columns
#' @param method Interpolation method: "idw", "kriging", or "spline"
#' @param extent terra::ext object
#' @param resolution Vector c(dx, dy)
#' @param crs Coordinate reference system
#'
#' @return terra::SpatRaster of interpolated surface
#'
#' @keywords internal
gc_interpolate_surface <- function(points, method, extent, resolution, crs) {
  # Create target raster
  target <- terra::rast(extent, res = resolution, crs = crs)

  # Convert points to sf
  points_sf <- sf::st_as_sf(points, coords = c("x", "y"), crs = crs)

  if (method == "idw") {
    # Inverse distance weighting
    interpolated <- gstat::idw(z ~ 1, locations = points_sf, newdata = target)
    result <- terra::rast(interpolated)
  } else if (method == "kriging") {
    # Simple kriging (fit variogram and krige)
    vgm <- gstat::vgm(psill = stats::var(points$z), model = "Exp", range = 100, nugget = 0)
    interpolated <- gstat::krige(z ~ 1, locations = points_sf, newdata = target, model = vgm)
    result <- terra::rast(interpolated)
  } else if (method == "spline") {
    # Thin plate spline - requires fields package
    if (!requireNamespace("fields", quietly = TRUE)) {
      stop("fields package required for spline interpolation. Install with install.packages('fields')")
    }
    fit <- fields::Tps(points[, c("x", "y")], points$z)
    grid_xy <- terra::xyFromCell(target, seq_len(terra::ncell(target)))
    z_pred <- predict(fit, grid_xy)
    terra::values(target) <- z_pred
    result <- target
  } else {
    stop("Unknown interpolation method: ", method)
  }

  return(result)
}


#' Transform Coordinates to Stratigraphic Space
#'
#' Convert (x, y, z) Cartesian coordinates to (x, y, z_strat) stratigraphic coordinates
#' using one of four transformation methods.
#'
#' @param data Data frame with x, y, z columns (z = depth in cm, negative downward)
#' @param surfaces Output from gc_define_surfaces()
#' @param transform_type Character: "proportional", "top_onlap", "bottom_erosional", "rotational"
#' @param h_const Numeric: Standard thickness for proportional transform (default 100 cm)
#' @param dip Numeric: Dip angle in degrees (0-90) for rotational transform
#' @param azimuth Numeric: Azimuth angle in degrees (0-360) for rotational transform
#' @param validate Logical: Check for points outside stratigraphic unit (default TRUE)
#'
#' @return Data frame with original columns plus z_strat column
#'
#' @details
#' **Transformation Formulas**:
#' 1. Proportional: z_strat = (z - z_base) / (z_top - z_base) × h_const
#' 2. Top-Onlap: z_strat = z_top - z
#' 3. Bottom-Erosional: z_strat = z - z_base
#' 4. Rotational: z_strat = R * [x, y, z]^T (rotation matrix from dip/azimuth)
#'
#' **Validation**:
#' - Points must be within stratigraphic unit: z_base <= z <= z_top
#' - Warns if points outside unit (can truncate or error)
#'
#' @examples
#' \dontrun{
#' # Proportional transform
#' surfaces <- gc_define_surfaces(0, -100, type = "constant", ...)
#' data_strat <- gc_stratigraphic_transform(
#'   data, surfaces,
#'   transform_type = "proportional",
#'   h_const = 100
#' )
#'
#' # Rotational for steep slope
#' data_strat <- gc_stratigraphic_transform(
#'   data, surfaces,
#'   transform_type = "rotational",
#'   dip = 15,
#'   azimuth = 180
#' )
#' }
#'
#' @export
gc_stratigraphic_transform <- function(data, surfaces,
                                      transform_type = "proportional",
                                      h_const = 100, dip = NULL, azimuth = NULL,
                                      validate = TRUE) {
  # Input validation
  stopifnot(is.data.frame(data))
  stopifnot(all(c("x", "y", "z") %in% names(data)))
  stopifnot(transform_type %in% c("proportional", "top_onlap", "bottom_erosional", "rotational"))

  if (transform_type == "rotational") {
    stopifnot(!is.null(dip), !is.null(azimuth))
    stopifnot(dip >= 0, dip <= 90)
    stopifnot(azimuth >= 0, azimuth <= 360)
  }

  # Extract surface values at data locations
  data_sf <- sf::st_as_sf(data, coords = c("x", "y"), crs = surfaces$crs, remove = FALSE)
  z_top <- terra::extract(surfaces$top, terra::vect(data_sf))[, 2]
  z_bottom <- terra::extract(surfaces$bottom, terra::vect(data_sf))[, 2]

  # Validation: check if points within unit
  if (validate) {
    outside <- data$z < z_bottom | data$z > z_top
    if (any(outside, na.rm = TRUE)) {
      n_outside <- sum(outside, na.rm = TRUE)
      warning(paste(n_outside, "points outside stratigraphic unit. Consider truncating or checking surfaces."))
    }
  }

  # Apply transformation
  if (transform_type == "proportional") {
    thickness <- z_top - z_bottom
    data$z_strat <- ((z_top - data$z) / thickness) * h_const

  } else if (transform_type == "top_onlap") {
    data$z_strat <- z_top - data$z

  } else if (transform_type == "bottom_erosional") {
    data$z_strat <- data$z - z_bottom

  } else {  # rotational
    # Build rotation matrix
    dip_rad <- dip * pi / 180
    azimuth_rad <- azimuth * pi / 180

    # Rotation around dip direction
    R <- matrix(c(
      cos(azimuth_rad), -sin(azimuth_rad), 0,
      sin(azimuth_rad) * cos(dip_rad), cos(azimuth_rad) * cos(dip_rad), sin(dip_rad),
      -sin(azimuth_rad) * sin(dip_rad), -cos(azimuth_rad) * sin(dip_rad), cos(dip_rad)
    ), nrow = 3, byrow = TRUE)

    # Apply rotation
    coords <- as.matrix(data[, c("x", "y", "z")])
    coords_rot <- t(R %*% t(coords))

    data$x_strat <- coords_rot[, 1]
    data$y_strat <- coords_rot[, 2]
    data$z_strat <- coords_rot[, 3]
  }

  # Store transformation metadata
  attr(data, "transform_type") <- transform_type
  attr(data, "surfaces") <- surfaces
  attr(data, "h_const") <- h_const
  if (transform_type == "rotational") {
    attr(data, "dip") <- dip
    attr(data, "azimuth") <- azimuth
  }

  return(data)
}


#' Transform from Stratigraphic to Cartesian Coordinates
#'
#' Inverse of gc_stratigraphic_transform() - converts (x, y, z_strat) back to (x, y, z).
#'
#' @param data_strat Data frame with stratigraphic coordinates (output from gc_stratigraphic_transform())
#'   Must have transformation metadata stored as attributes.
#'
#' @return Data frame with x, y, z columns restored to Cartesian coordinates
#'
#' @details
#' Reads transformation type and parameters from attributes, applies inverse transformation.
#'
#' @export
gc_inverse_stratigraphic_transform <- function(data_strat) {
  # Extract metadata
  transform_type <- attr(data_strat, "transform_type")
  surfaces <- attr(data_strat, "surfaces")
  h_const <- attr(data_strat, "h_const")

  if (is.null(transform_type) || is.null(surfaces)) {
    stop("data_strat must have transformation metadata. Use output from gc_stratigraphic_transform().")
  }

  # Extract surface values
  data_sf <- sf::st_as_sf(data_strat, coords = c("x", "y"),
                         crs = surfaces$crs, remove = FALSE)
  z_top <- terra::extract(surfaces$top, terra::vect(data_sf))[, 2]
  z_bottom <- terra::extract(surfaces$bottom, terra::vect(data_sf))[, 2]

  # Inverse transformation
  if (transform_type == "proportional") {
    thickness <- z_top - z_bottom
    data_strat$z <- z_top - (data_strat$z_strat / h_const) * thickness

  } else if (transform_type == "top_onlap") {
    data_strat$z <- z_top - data_strat$z_strat

  } else if (transform_type == "bottom_erosional") {
    data_strat$z <- data_strat$z_strat + z_bottom

  } else {  # rotational
    dip <- attr(data_strat, "dip")
    azimuth <- attr(data_strat, "azimuth")

    # Inverse rotation matrix (transpose of original)
    dip_rad <- dip * pi / 180
    azimuth_rad <- azimuth * pi / 180

    R <- matrix(c(
      cos(azimuth_rad), -sin(azimuth_rad), 0,
      sin(azimuth_rad) * cos(dip_rad), cos(azimuth_rad) * cos(dip_rad), sin(dip_rad),
      -sin(azimuth_rad) * sin(dip_rad), -cos(azimuth_rad) * sin(dip_rad), cos(dip_rad)
    ), nrow = 3, byrow = TRUE)

    R_inv <- t(R)

    coords_strat <- as.matrix(data_strat[, c("x_strat", "y_strat", "z_strat")])
    coords_cart <- t(R_inv %*% t(coords_strat))

    data_strat$x <- coords_cart[, 1]
    data_strat$y <- coords_cart[, 2]
    data_strat$z <- coords_cart[, 3]
  }

  return(data_strat)
}
