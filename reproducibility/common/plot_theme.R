# =============================================================================
# plot_theme.R
#
# The style of every figure:
#   * large fonts (base size 14) and capitalized axis titles.
#   * a fixed color, point shape, line type, line width and point size for
#     each method, so that a method looks the same in every figure.
#     LEWIS (complexity) is a heavy black dashed line and LEWIS (uniform) a
#     heavy red dot-dashed line. The other methods are light solid lines in
#     colors that are safe for color-blind readers (Paul Tol's palettes).
#   * save_figure() writes a vector PDF and a 300-dpi PNG. It removes the
#     figure title, subtitle and caption, because the paper gives the
#     caption in LaTeX.
#
# Usage. Map all five method aesthetics in aes(), give the geoms no fixed
# size or linewidth, and add scale_method_all() and theme_lewis():
#
#     p <- ggplot(d, aes(x, y, color = Method, shape = Method,
#                        linetype = Method, linewidth = Method,
#                        size = Method, group = Method)) +
#            geom_line() + geom_point() +
#            scale_method_all() + theme_lewis()
#     save_figure(p, file.path(lewis_fig_dir("comparison"), "name"),
#                 width = 11, height = 5.5)
#
# The calling script must define COMMON_DIR before it sources this file.
# Dependency: ggplot2 (and patchwork for composite figures).
# =============================================================================

stopifnot(exists("COMMON_DIR"))
source(file.path(COMMON_DIR, "settings.R"))
source(file.path(COMMON_DIR, "paths.R"))

# Load ggplot2.
suppressPackageStartupMessages(library(ggplot2))

# Colors. Black and red are used only by the two versions of LEWIS.
METHOD_COLORS <- c(
  "LEWIS (complexity)" = "#000000",
  "LEWIS (uniform)"    = "#CC3311",
  "PC"                 = "#DDCC77",
  "GES"                = "#44AA99",
  "ICA-LiNGAM"         = "#88CCEE",
  "DirectLiNGAM"       = "#4477AA",
  "RCD"                = "#117733",
  "TL"                 = "#999933",
  "NOTEARS"            = "#AA4499",
  "Pairwise-LiNGAM"    = "#661100"
)

# Point shapes, a second channel that does not depend on color.
METHOD_SHAPES <- c(
  "LEWIS (complexity)" = 16,
  "LEWIS (uniform)"    = 1,
  "PC"                 = 17,
  "GES"                = 15,
  "ICA-LiNGAM"         = 3,
  "DirectLiNGAM"       = 4,
  "RCD"                = 18,
  "TL"                 = 8,
  "NOTEARS"            = 6,
  "Pairwise-LiNGAM"    = 7
)

# Line types. LEWIS (complexity) uses the short dash "22" so that the
# legend key shows a dash on each side of the marker.
METHOD_LINETYPES <- c(
  "LEWIS (complexity)" = "22",
  "LEWIS (uniform)"    = "dotdash",
  "PC"                 = "solid",
  "GES"                = "solid",
  "ICA-LiNGAM"         = "solid",
  "DirectLiNGAM"       = "solid",
  "RCD"                = "solid",
  "TL"                 = "solid",
  "NOTEARS"            = "solid",
  "Pairwise-LiNGAM"    = "solid"
)

# Colors for groups that are not methods, for example the values of p in
# the calibration figures (Paul Tol's bright palette).
GROUP_COLORS <- c("#4477AA", "#EE6677", "#228833", "#CCBB44",
                   "#66CCEE", "#AA3377")

# Line widths and point sizes. The two versions of LEWIS are drawn heavier.
.IS_LEWIS <- names(METHOD_COLORS) %in% c("LEWIS (complexity)", "LEWIS (uniform)")
METHOD_LINEWIDTHS <- setNames(ifelse(.IS_LEWIS, 1.7, 0.65), names(METHOD_COLORS))
METHOD_POINTSIZES <- setNames(ifelse(.IS_LEWIS, 2.8, 1.7), names(METHOD_COLORS))

# The five manual scales of a method plot, merged into one legend with two
# rows. With fill = TRUE the fill scale is added as well (used by the SID
# ribbons).
scale_method_all <- function(fill = FALSE) {
  leg <- guide_legend(nrow = 2)
  out <- list(scale_color_manual(values = METHOD_COLORS),
              scale_shape_manual(values = METHOD_SHAPES),
              scale_linetype_manual(values = METHOD_LINETYPES),
              scale_linewidth_manual(values = METHOD_LINEWIDTHS),
              scale_size_manual(values = METHOD_POINTSIZES),
              guides(color = leg, shape = leg, linetype = leg,
                     linewidth = leg, size = leg))
  if (isTRUE(fill)) {
    out <- c(out, list(scale_fill_manual(values = METHOD_COLORS)))
  }
  out
}

# The theme: large fonts, bold titles, dark facet strips, legend at the
# bottom.
theme_lewis <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = rel(1.15)),
      plot.subtitle    = element_text(size = rel(0.95)),
      strip.background = element_rect(fill = "#222222"),
      strip.text       = element_text(color = "white", face = "bold",
                                      size = rel(0.95)),
      axis.title       = element_text(size = rel(1.05)),
      axis.text        = element_text(size = rel(0.95)),
      legend.position  = "bottom",
      legend.text      = element_text(size = rel(0.95)),
      legend.title     = element_blank(),
      legend.key.width = grid::unit(2.6, "lines"),
      legend.spacing.x = grid::unit(0.6, "lines"),
      panel.grid.minor = element_blank()
    )
}

# Remove the figure title, subtitle and caption. A patchwork composite keeps
# them in plot_annotation(). The subtitles of its panels are kept, because
# they act as row labels.
strip_figure_titles <- function(plot) {
  if (inherits(plot, "patchwork")) {
    plot <- plot + patchwork::plot_annotation(title = NULL, subtitle = NULL,
                                              caption = NULL)
  } else {
    plot$labels$title <- NULL
    plot$labels$subtitle <- NULL
    plot$labels$caption <- NULL
  }
  plot
}

# Write <stem>.pdf (vector) and <stem>.png (300 dpi).
save_figure <- function(plot, stem, width, height, strip_titles = TRUE) {
  if (isTRUE(strip_titles)) plot <- strip_figure_titles(plot)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         device = cairo_pdf)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = 300)
  cat(sprintf("Wrote %s.pdf and %s.png\n", stem, stem))
}
