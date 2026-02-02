#' Discretize Block into Sub-Points
#'
#' Subdivide a block (grid cell) into internal sub-points for covariance calculation.
#' Uses regular grid points centered within the block for uniform coverage.
#'
#' @param block_centers Data frame with x, y columns (block center coordinates)
#' @param block_size Vector c(dx, dy) specifying block dimensions (same units as coordinates)
#' @param discretization Vector c(nx, ny) specifying number of sub-points per dimension
#'   (default c(4, 4) gives 16 sub-points per block)
#'
#' @return Data frame with columns:
#'   - x, y: Sub-point coordinates
#'   - block_id: Index linking sub-point to parent block
#'
#' @details
#' **Discretization strategy:**
#' - Blocks are subdivided into regular grid of sub-points
#' - Sub-points are centered within each sub-cell (not at corners)
#' - Example: 4×4 discretization creates 16 evenly-spaced sub-points per block
#' - Total sub-points = nrow(block_centers) × (nx × ny)
#'
#' **Block-support covariance computation:**
#' Block covariances are computed via double integration using sub-points:
#' ```
#' C_BB(h) = (1/V²) ∑∑ C(u_i - v_j) (approximate via sub-points)
#' ```
#'
#' @examples
#' \dontrun{
#' # Create a 3×3 grid of blocks centered at regular intervals
#' blocks <- expand.grid(x = c(0, 100, 200), y = c(0, 100, 200))
#'
#' # Discretize each block into 4×4 sub-points
#' sub_points <- gc_discretize_block(blocks, block_size = c(100, 100),
#'                                    discretization = c(4, 4))
#'
#' # Total: 9 blocks × 16 sub-points = 144 sub-points
#' nrow(sub_points)  # 144
#'
#' # Sub-points within first block
#' sub_points[sub_points$block_id == 1, ]
#' }
#'
#' @export
gc_discretize_block <- function(block_centers, block_size, discretization = c(4, 4)) {
  stopifnot(is.data.frame(block_centers), all(c("x", "y") %in% names(block_centers)))
  stopifnot(is.numeric(block_size), length(block_size) == 2, all(block_size > 0))
  stopifnot(is.numeric(discretization), length(discretization) == 2, all(discretization >= 2))

  nx <- discretization[1]
  ny <- discretization[2]
  dx_block <- block_size[1]
  dy_block <- block_size[2]

  # Sub-cell dimensions
  dx_sub <- dx_block / nx
  dy_sub <- dy_block / ny

  # Sub-point offsets within block (centered in each sub-cell)
  # Range: [-block_size/2 + sub_size/2, ..., block_size/2 - sub_size/2]
  x_offsets <- seq(-dx_block / 2 + dx_sub / 2, dx_block / 2 - dx_sub / 2, length.out = nx)
  y_offsets <- seq(-dy_block / 2 + dy_sub / 2, dy_block / 2 - dy_sub / 2, length.out = ny)

  # Generate sub-points for all blocks
  sub_points_list <- lapply(seq_len(nrow(block_centers)), function(i) {
    grid <- expand.grid(
      x_offset = x_offsets,
      y_offset = y_offsets
    )
    data.frame(
      x = block_centers$x[i] + grid$x_offset,
      y = block_centers$y[i] + grid$y_offset,
      block_id = i
    )
  })

  sub_points <- do.call(rbind, sub_points_list)
  rownames(sub_points) <- NULL

  return(sub_points)
}


#' Average Sub-Point ILR Values to Block Support
#'
#' Compute block-averaged ILR values from simulated sub-point results.
#' Averages across all sub-points within each block for each simulation realization.
#'
#' @param sim_result Simulation result data frame with sub-point locations and ILR values.
#'   Columns should include ILR simulation results (e.g., ilr1.sim1, ilr2.sim1, etc.)
#' @param sub_points Data frame from [gc_discretize_block()] with block_id column
#'   linking each row to its parent block
#'
#' @return Data frame with averaged ILR values per block per realization
#'   - Rows: One per block
#'   - Columns: ilr1.sim1, ilr1.sim2, ..., ilr2.sim1, ilr2.sim2, ...
#'   (or equivalent MAF columns if MAF transformation was used)
#'
#' @details
#' **Averaging operation:**
#' For each block B and simulation S:
#' ```
#' ILR_B^S = mean(ILR_ij^S : (i,j) ∈ sub-points of B)
#' ```
#'
#' This averaging is what creates the volume-variance relationship:
#' - Block variance = Point variance / (number of sub-points)
#' - Blocks are smoother than points (lower variance)
#'
#' @examples
#' \dontrun{
#' # After discretizing and simulating at sub-point support
#' sub_points <- gc_discretize_block(blocks, block_size = c(100, 100))
#' sim_points <- gc_sim_composition(model, sub_points, nsim = 10)
#'
#' # Average to block support
#' block_sims <- gc_average_block_ilr(sim_points, sub_points)
#' }
#'
#' @export
gc_average_block_ilr <- function(sim_result, sub_points) {
  stopifnot(is.data.frame(sim_result), is.data.frame(sub_points))
  stopifnot("block_id" %in% names(sub_points))
  stopifnot(nrow(sim_result) == nrow(sub_points))

  # Add block_id to simulation results
  sim_result$block_id <- sub_points$block_id

  # Identify ILR or MAF columns (columns matching ilr/maf pattern with .simN suffix)
  ilr_cols <- grep("^(ilr|maf)\\d+\\.sim\\d+$", names(sim_result), value = TRUE)

  if (length(ilr_cols) == 0) {
    stop("No ILR/MAF simulation columns found (expected pattern: ilr1.sim1, maf1.sim1, etc.)")
  }

  # Group by block_id and average
  block_averaged <- aggregate(
    sim_result[, ilr_cols, drop = FALSE],
    by = list(block_id = sim_result$block_id),
    FUN = mean
  )

  return(block_averaged)
}
