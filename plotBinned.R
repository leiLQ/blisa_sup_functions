# Plot gene expression or cell-type composition from a BLISA binned object.
# Requires sf and ggplot2; sparse expression matrices also require Matrix.

plotBinned <- function(x, gene = NULL, cell_type = NULL,
                       cell_stat = c("proportion", "number"),
                       log1p = FALSE, as_points = FALSE, size = 1.5,
                       title = NULL,
                       colours = NULL) {
  if (is.null(gene) == is.null(cell_type)) {
    stop("Supply exactly one of 'gene' or 'cell_type'.")
  }
  feature <- if (is.null(gene)) cell_type else gene
  if (!is.character(feature) || length(feature) != 1L || is.na(feature)) {
    stop("Select one gene or cell-type label.")
  }
  cell_stat <- match.arg(cell_stat)
  if (is.null(colours)) {
    colours <- if (!is.null(gene)) {
      # Match scider::plotGrid's continuous default: rev(scider:::col.spec).
      # Keep the exact colour stops here to avoid depending on a private object.
      c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
        "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F")
    } else {
      c("#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")
    }
  }
  bins <- x$bins
  if (!inherits(bins, "sf") || nrow(bins) == 0L ||
      is.null(bins$bin_id) || anyNA(bins$bin_id) || anyDuplicated(bins$bin_id)) {
    stop("x$bins must contain sf polygons with unique, non-missing bin_id values.")
  }
  bin_ids <- paste0("bin_factor", bins$bin_id)
  limits <- NULL

  if (!is.null(gene)) {
    counts <- x$counts_matrix
    if (sum(rownames(counts) == gene) != 1L) {
      stop("Gene must occur exactly once in x$counts_matrix: ", gene)
    }
    if (is.null(colnames(counts)) || anyDuplicated(colnames(counts))) {
      stop("Expression bin names must be present and unique.")
    }
    index <- match(bin_ids, colnames(counts))
    if (anyNA(index)) stop("Some polygons have no matching expression bin.")
    values <- as.numeric(counts[gene, index, drop = FALSE])
    legend_title <- "Transcript count"
    if (log1p) {
      if (any(!is.finite(values)) || any(values < 0)) {
        stop("log1p requires finite, non-negative gene counts.")
      }
      values <- base::log1p(values)
      legend_title <- "log1p(count)"
    }
  } else {
    if (log1p) stop("log1p is only supported for gene expression.")
    metadata <- x$bins_metadata
    prefix <- if (cell_stat == "proportion") "cell_proportion::" else "cell_count::"
    column <- paste0(prefix, cell_type)
    if (is.null(metadata) || !column %in% colnames(metadata)) {
      stop("Cell-type metadata column not found: ", column)
    }
    if (is.null(metadata$sample_id) || anyNA(metadata$sample_id) ||
        anyDuplicated(metadata$sample_id)) {
      stop("bins_metadata must have unique, non-missing sample_id values.")
    }
    index <- match(bin_ids, metadata$sample_id)
    if (anyNA(index)) stop("Some polygons have no matching metadata row.")
    values <- metadata[[column]][index]
    legend_title <- if (cell_stat == "proportion") "Cell proportion" else "Cell count"
    if (cell_stat == "proportion") limits <- c(0, 1)
  }
  if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0)) {
    stop("Plot values must be finite, non-negative numbers.")
  }
  if (!is.null(limits) && any(values > 1)) {
    stop("Cell-type proportions must lie between zero and one.")
  }
  bins$.bin_value <- values

  if (as_points) {
    centres <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(bins)))
    points <- data.frame(x = centres[, 1], y = centres[, 2], .bin_value = values)
    p <- ggplot2::ggplot(points, ggplot2::aes(x = x, y = y, fill = .bin_value)) +
      ggplot2::geom_point(shape = 21, colour = NA, size = size) +
      ggplot2::coord_equal()
  } else {
    p <- ggplot2::ggplot(bins) +
      ggplot2::geom_sf(ggplot2::aes(fill = .bin_value), colour = NA)
  }
  p + ggplot2::scale_fill_gradientn(colours = colours, limits = limits,
                                    name = legend_title) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top")) +
    ggplot2::labs(title = if (is.null(title)) feature else title) +
    ggplot2::theme_void() + ggplot2::theme(legend.position = "right")
}
