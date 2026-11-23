# =============================================================================
# replot_uq_figures.R
#
# Draws the two posterior displays of the uncertainty quantification study
# from the saved tables. No chain is run. It takes seconds on a laptop.
#
#   posterior_marginal_annotated   Figure 17 (Appendix H): the posterior
#       mass of the 15 most probable RPDAGs, colored by class, each with its
#       edits relative to D_hat
#   hasse_partition_compact        Figure 11(b) (Section 7.3): one bubble per
#       graph at (SHD from D_hat, edges minus edges of D_hat), with the
#       radius of the credible ball
#
# Input (in results/tables/ of this folder):
#   posterior_uncertainty_table.csv   columns rank, mass, category,
#       shd_to_Dhat, n_edges_added, n_edges_removed, n_edges_reversed,
#       added_edges_str, removed_edges_str, reversed_edges_str,
#       nG_added_str, nG_removed_str, and n_edges (optional, used for a
#       check)
#   decision_rule_summary.csv          columns eps_alpha and alpha_ball, and
#       D_hat_edges (optional, used for a check)
#
# Format of the edit lists. Entries are separated by "|". A directed edge is
# "a->b" and an undirected edge is "a - b". A reversed edge is written as its
# orientation in D_hat, then "=>", then its orientation in the graph, for
# example "8->1 => 1->8". The tables of the paper in results/tables/ use two
# compact forms, which are also read: "a-b" for an undirected edge, and
# "a<->b" for a reversed edge that is a->b in the graph and b->a in D_hat.
# The sets of non-Gaussian nodes that a graph adds or removes are separated
# by ",".
#
# Output: <figure folder>/uq/, where the figure folder is the environment
# variable LEWIS_FIG_DIR or reproducibility/figures/ (see common/paths.R).
#
# Run:   Rscript reproducibility/04_uq/replot_uq_figures.R
#        (from the root of the repository, or from any other folder)
# =============================================================================

# ---- Folder of this script and shared code ----------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "..", "common")
source(file.path(COMMON_DIR, "plot_theme.R"))

# ---- Folders -----------------------------------------------------------------
TAB_DIR <- file.path(.here, "results", "tables")
FIG_DIR <- lewis_fig_dir("uq")

# ---- Read the tables -----------------------------------------------------------
uq <- read.csv(file.path(TAB_DIR, "posterior_uncertainty_table.csv"),
               stringsAsFactors = FALSE)
summ <- read.csv(file.path(TAB_DIR, "decision_rule_summary.csv"),
                 stringsAsFactors = FALSE)
need <- c("rank", "mass", "category", "shd_to_Dhat", "n_edges_added",
          "n_edges_removed", "n_edges_reversed", "added_edges_str",
          "removed_edges_str", "reversed_edges_str", "nG_added_str",
          "nG_removed_str")
miss <- setdiff(need, names(uq))
if (length(miss) > 0L) {
  stop("Columns missing in posterior_uncertainty_table.csv: ", paste(miss, collapse = ", "))
}
if (!all(c("eps_alpha", "alpha_ball") %in% names(summ))) {
  stop("decision_rule_summary.csv must have the columns eps_alpha and alpha_ball.")
}
eps_alpha <- as.integer(summ$eps_alpha[1])
alpha_ball <- as.numeric(summ$alpha_ball[1])
uq <- uq[order(uq$rank), ]
cat(sprintf("%d graphs read. Radius %d at alpha_ball = %g\n", nrow(uq), eps_alpha, alpha_ball))

# ---- Classes of Step 4 ----------------------------------------------------------
CLASS_LEVELS <- c("identical", "upper", "lower", "rearranged")
CLASS_COLORS <- c(identical = "#8B5CF6", upper = "#3B82F6",
                   rearranged = "#F59E0B", lower = "#EF4444")
if (!all(uq$category %in% CLASS_LEVELS)) stop("Unknown class in posterior_uncertainty_table.csv.")
CLASSES_SEEN <- CLASS_LEVELS[CLASS_LEVELS %in% unique(uq$category)]

# ---- Edit labels -----------------------------------------------------------------
# Empty string for a missing value.
.na0 <- function(x) ifelse(is.na(x), "", as.character(x))
# Split an edit list into its entries.
.entries <- function(x) {
  x <- .na0(x)
  if (!nzchar(x)) return(character(0))
  trimws(strsplit(x, "|", fixed = TRUE)[[1]])
}
# Format the entries of an added or removed list. A compact undirected edge
# "a-b" between numbered nodes becomes "a - b".
.fmt_edges <- function(x) {
  e <- .entries(x)
  e <- sub("^([0-9]+)-([0-9]+)$", "\\1 - \\2", e)
  paste(e, collapse = " | ")
}
# Format the entries of a reversed list. A compact entry "a<->b" becomes
# "b->a => a->b" (the orientation in D_hat first).
.fmt_reversed <- function(x) {
  e <- .entries(x)
  e <- ifelse(grepl("<->", e, fixed = TRUE),
              sub("^(.+)<->(.+)$", "\\2->\\1 => \\1->\\2", e), e)
  paste(e, collapse = " | ")
}
# Format a set of non-Gaussian nodes that a graph adds (sign "+") or
# removes (sign "-") relative to nG_hat.
.fmt_nodes <- function(x, sign) {
  v <- trimws(strsplit(.na0(x), "[|,]")[[1]])
  v <- v[nzchar(v)]
  paste0(sign, v, collapse = " | ")
}
# The label of one graph.
edit_label <- function(i) {
  bits <- character(0)
  if (uq$n_edges_added[i] > 0L) bits <- c(bits, paste("+e:", .fmt_edges(uq$added_edges_str[i])))
  if (uq$n_edges_removed[i] > 0L) bits <- c(bits, paste("-e:", .fmt_edges(uq$removed_edges_str[i])))
  if (uq$n_edges_reversed[i] > 0L) bits <- c(bits, paste("rev:", .fmt_reversed(uq$reversed_edges_str[i])))
  if (nzchar(.na0(uq$nG_removed_str[i]))) bits <- c(bits, paste("nG:", .fmt_nodes(uq$nG_removed_str[i], "-")))
  if (nzchar(.na0(uq$nG_added_str[i]))) bits <- c(bits, paste("nG:", .fmt_nodes(uq$nG_added_str[i], "+")))
  if (length(bits) == 0L) return("(identical)")
  paste(bits, collapse = "   ")
}
uq$edit_label <- vapply(seq_len(nrow(uq)), edit_label, character(1))

# ---- Line breaks in long labels --------------------------------------------------
# The width of a label is estimated as a share of the panel width from its
# characters: an upper-case letter or a digit takes W_WIDE, any other
# character W_NARROW. A label longer than the room beyond its bar is broken
# between its parts (+e, -e, rev, nG). The room is never less than MIN_ROOM.
W_WIDE <- 1.03 * 0.0109
W_NARROW <- 1.03 * 0.0056
MIN_ROOM <- 0.35
.label_width <- function(s) {
  ch <- strsplit(s, "")[[1]]
  sum(ifelse(grepl("[A-Z0-9]", ch), W_WIDE, W_NARROW))
}
.wrap_label <- function(s, room) {
  if (.label_width(s) <= room) return(s)
  parts <- strsplit(s, "   ", fixed = TRUE)[[1]]
  lines <- parts[1]
  for (q in parts[-1]) {
    cand <- paste(lines[length(lines)], q, sep = "   ")
    if (.label_width(cand) <= room) lines[length(lines)] <- cand else lines <- c(lines, q)
  }
  paste(lines, collapse = "\n")
}

# ---- Figure 17: posterior mass with the edits ------------------------------------
K <- min(15L, nrow(uq))
pm <- uq[seq_len(K), ]
# The mass axis runs to 1.5 times the largest mass.
room <- pmax(MIN_ROOM, 0.97 - pm$mass / (1.5 * max(pm$mass)))
pm$bar_label <- mapply(.wrap_label, sprintf("%.1f%%   %s", 100 * pm$mass, pm$edit_label),
                       room, USE.NAMES = FALSE)
pm$category <- factor(pm$category, levels = CLASSES_SEEN)
# Each graph gets a band of the vertical axis as tall as its label needs:
# one unit for the first line and one unit for each further line.
pm$n_lines <- 1L + lengths(regmatches(pm$bar_label, gregexpr("\n", pm$bar_label, fixed = TRUE)))
band <- pm$n_lines
pm$pos <- cumsum(band) - band / 2
total_band <- sum(band)

# Draw the bars and their labels.
g_marg <- ggplot(pm, aes(x = pos, y = mass, fill = category)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.3) +
  geom_text(aes(label = bar_label), hjust = -0.03, vjust = 0.5, size = 4.6,
            lineheight = 0.9, color = "grey20") +
  scale_fill_manual(values = CLASS_COLORS[CLASSES_SEEN], breaks = CLASSES_SEEN,
                    labels = CLASSES_SEEN,
                    name = expression(paste("Class relative to ", hat(D), ":  "))) +
  scale_x_continuous(breaks = pm$pos, labels = as.character(pm$rank),
                     expand = expansion(add = 0.4)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = scales::breaks_pretty(n = 4),
                     expand = expansion(mult = c(0, 0.5))) +
  coord_flip() +
  labs(x = "RPDAG, ranked by posterior mass", y = "Posterior mass") +
  theme_lewis() +
  theme(legend.position = "top",
        legend.title = element_text(size = rel(0.95)),
        legend.key.width = grid::unit(1.2, "lines"),
        panel.grid.minor.y = element_blank())

# The height is 5.2 inches for 15 one-line labels, of which 1.1 inches are
# axes and legend.
save_figure(g_marg, file.path(FIG_DIR, "posterior_marginal_annotated"),
            width = 12, height = 1.1 + 4.1 * total_band / 15)

# ---- Figure 11(b): SHD from D_hat and change in the number of edges ---------------
ht <- uq[seq_len(K), ]
# The vertical position is the number of edges of the graph minus that of
# D_hat. An added edge counts +1 and a removed edge -1. A reversal does not
# change the number of edges. An edge that changes between directed and
# undirected is one removal and one addition, so it counts 0.
ht$d_edges <- ht$n_edges_added - ht$n_edges_removed

# Check the change in the number of edges against the edge counts, when the
# tables give them.
e_hat <- if ("D_hat_edges" %in% names(summ)) summ$D_hat_edges[1] else NA
if ("n_edges" %in% names(ht) && !is.na(e_hat)) {
  if (!isTRUE(all(ht$n_edges - e_hat == ht$d_edges))) {
    stop("The edit lists disagree with the numbers of edges of the graphs.")
  }
  cat(sprintf("Edit lists agree with the numbers of edges of all %d graphs.\n", nrow(ht)))
} else {
  cat("The tables give no edge counts. The check of the edit lists is skipped.\n")
}

# The class and the position of each bubble.
ht$category <- factor(ht$category, levels = CLASSES_SEEN)
ht$xval <- ht$shd_to_Dhat
ht$x_pos <- ht$xval
ht$y_pos <- ht$d_edges
# Graphs at the same position are moved slightly apart along a short
# diagonal, in order of rank. The offset carries no meaning.
.cell <- paste(ht$xval, ht$d_edges)
for (cc in unique(.cell)) {
  k <- which(.cell == cc)
  if (length(k) > 1L) {
    k <- k[order(ht$rank[k])]
    dy <- min(0.15 + 0.075 * (length(k) - 2L), 0.30)
    dx <- min(0.24 + 0.18 * (length(k) - 2L), 0.60)
    ht$y_pos[k] <- ht$d_edges[k] + seq(-dy, dy, length.out = length(k))
    ht$x_pos[k] <- ht$xval[k] + seq(-dx, dx, length.out = length(k))
  }
}
# A moved bubble gets its rank beside it, the others above it.
ht$lab_h <- 0.5
ht$lab_v <- -1.2
for (cc in unique(.cell)) {
  k <- which(.cell == cc)
  if (length(k) > 1L) {
    ht$lab_h[k] <- -1.0
    ht$lab_v[k] <- 0.5
  }
}
# A bubble with at least a quarter of the largest mass carries its rank
# inside, in white.
ht$inside <- ht$mass >= 0.25 * max(ht$mass)
y_lo <- min(ht$y_pos) - 0.6
y_hi <- max(ht$y_pos) + 0.6
band_max <- max(ht$xval, eps_alpha) + 1

# Draw the bubbles, the radius of the ball and the region outside it.
g_hasse <- ggplot(ht, aes(x = x_pos, y = y_pos, size = mass, color = category)) +
  annotate("rect", xmin = eps_alpha + 0.5, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#FEE2E2", alpha = 0.4) +
  geom_vline(xintercept = eps_alpha + 0.5, linetype = "dashed",
             color = "#7C3AED", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.85) +
  geom_text(data = ht[!ht$inside, ],
            aes(label = sprintf("%d", rank), hjust = lab_h, vjust = lab_v),
            size = 4.8, color = "black", show.legend = FALSE) +
  geom_text(data = ht[ht$inside, ], aes(label = sprintf("%d", rank)),
            size = 5.2, color = "white", fontface = "bold", show.legend = FALSE) +
  scale_size_area(max_size = 14, labels = scales::percent_format(accuracy = 1),
                  name = "Posterior mass") +
  scale_color_manual(values = CLASS_COLORS[CLASSES_SEEN], breaks = CLASSES_SEEN,
                      labels = CLASSES_SEEN,
                      name = expression(paste("Class relative to ", hat(D)))) +
  scale_x_continuous(breaks = seq(0L, as.integer(band_max) + 1L, by = 1L),
                     limits = c(-0.6, band_max + 0.5)) +
  scale_y_continuous(breaks = seq(ceiling(y_lo), floor(y_hi), by = 1L),
                     limits = c(y_lo, y_hi)) +
  annotate("text", x = eps_alpha, y = y_hi,
           label = sprintf("epsilon[%g] == %d", alpha_ball, eps_alpha), parse = TRUE,
           vjust = 1.2, hjust = 0.5, color = "#5B21B6", size = 5.2, fontface = "bold") +
  labs(x = expression(SHD(D, hat(D))),
       y = expression(paste("Edges of ", D, " minus edges of ", hat(D)))) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_lewis() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size = rel(0.95)),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = grid::unit(0, "pt"))

# Write the figure.
save_figure(g_hasse, file.path(FIG_DIR, "hasse_partition_compact"),
            width = 8.2, height = 5.6)
