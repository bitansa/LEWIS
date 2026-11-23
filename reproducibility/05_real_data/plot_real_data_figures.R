# =============================================================================
# plot_real_data_figures.R
#
# The figures of Appendix I, drawn from the saved tables of the reported
# analyses (RD_REPORTED in real_data_config.R) only:
#
#   real_data_uq.pdf
#       Figure 18: the posterior mass of the RPDAGs by their SHD from
#       D_hat, split into the classes of Step 4 of the post-processing, with
#       the radius of the 95 percent credible ball. One panel per analysis.
#   posterior_marginal_annotated_<tag>.pdf
#       Figure 19: the posterior mass of the 15 most probable
#       RPDAGs of each analysis, colored by class and labeled with the edges
#       in which they differ from D_hat. The layout is that of Figure 17
#       of Appendix H.
#
# The figures go to lewis_fig_dir("real_data", "TCGA"), that is
# <LEWIS_FIG_DIR>/real_data/TCGA/ (default reproducibility/figures/). The
# masses of Figure 18 are also written to
# results/tables/real_data_uq_by_shd.csv. Any error stops the script.
#
# Run:  Rscript plot_real_data_figures.R
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))
options(warn = 1)

# ---- Folder of this script, settings and style -------------------------------
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) normalizePath(dirname(fa[1]), mustWork = FALSE)
  else getwd()
})
COMMON_DIR <- file.path(.here, "..", "common")
Sys.setenv(RD_OUT_DIR = file.path(.here, "results"))
source(file.path(.here, "real_data_config.R"))
source(file.path(COMMON_DIR, "plot_theme.R"))

# The folder of the figures.
FIG_DIR <- lewis_fig_dir("real_data", "TCGA")

# ---- The reported analyses ---------------------------------------------------
tags <- vapply(seq_len(nrow(RD_REPORTED)), function(k) {
  rd_tag(RD_REPORTED$network[k], RD_REPORTED$n_iter[k])
}, character(1))
# Panel titles: the label of the network without the number of genes.
panel_names <- sprintf("(%s) %s", letters[seq_along(tags)],
                       sub(" \\([0-9]+ genes.*\\)$", "",
                           vapply(RD_REPORTED$network,
                                  function(z) RD_NETWORKS[[z]]$label, character(1))))

# Read a table of an analysis. Stop if it is missing.
read_tab <- function(path) {
  if (!file.exists(path)) stop("Missing file ", path, ". Run run_reported_analyses.R first.")
  read.csv(path, stringsAsFactors = FALSE)
}

# Classes of Step 4, in the order and with the colors of Appendix H.
CLASS_LEVELS <- c("identical", "upper", "lower", "rearranged")
CLASS_COLORS <- c(identical = "#8B5CF6", upper = "#3B82F6",
                   rearranged = "#F59E0B", lower = "#EF4444")

# =============================================================================
# 1. Figure 18: posterior mass by SHD from D_hat
# =============================================================================
uq <- do.call(rbind, lapply(seq_along(tags), function(k) {
  post <- read_tab(file.path(RD_TAB_DIR, sprintf("posterior_rpdag_%s.csv", tags[k])))
  summ <- read_tab(file.path(RD_TAB_DIR, sprintf("uq_%s", tags[k]), "decision_rule_summary.csv"))
  # Sum the mass by SHD and class.
  a <- aggregate(mass ~ shd_to_Dhat + category, data = post, FUN = sum)
  a$panel <- panel_names[k]
  a$eps <- summ$eps_alpha[1]
  a
}))
uq$category <- factor(uq$category, levels = CLASS_LEVELS)
uq$panel <- factor(uq$panel, levels = panel_names)
write.csv(uq[order(uq$panel, uq$shd_to_Dhat, uq$category), ],
          file.path(RD_TAB_DIR, "real_data_uq_by_shd.csv"), row.names = FALSE)

# The axis ends at the largest SHD with at least 0.1 percent of the mass in
# some panel, and at least one step beyond the radius.
eps_lab <- unique(uq[, c("panel", "eps")])
tot <- aggregate(mass ~ shd_to_Dhat, data = uq, FUN = max)
x_max <- max(tot$shd_to_Dhat[tot$mass >= 0.001], eps_lab$eps + 1)
uq <- uq[uq$shd_to_Dhat <= x_max, ]
ymax <- max(aggregate(mass ~ panel + shd_to_Dhat, data = uq, FUN = sum)$mass)

# Draw the stacked bars, the radius of the ball and the region outside it.
g_uq <- ggplot(uq, aes(shd_to_Dhat, 100 * mass, fill = category)) +
  geom_rect(data = eps_lab, aes(xmin = eps + 0.5, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "#FEE2E2", alpha = 0.4, inherit.aes = FALSE) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  geom_vline(data = eps_lab, aes(xintercept = eps + 0.5), linetype = "dashed",
             color = "#7C3AED", linewidth = 0.6) +
  geom_text(data = eps_lab,
            aes(x = eps + 0.4, y = 100 * ymax * 0.98,
                label = sprintf("epsilon[0.05] == %d", eps)),
            parse = TRUE, hjust = 1, vjust = 1, size = 8 / .pt, color = "#5B21B6",
            fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = CLASS_COLORS, drop = FALSE,
                    name = expression(paste("Class relative to ", hat(D), ":  "))) +
  scale_x_continuous(breaks = 0:x_max, limits = c(-0.6, x_max + 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(SHD(D, hat(D))), y = "Posterior mass (%)") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 8),
        plot.margin = margin(2, 4, 2, 2))
save_figure(g_uq, file.path(FIG_DIR, "real_data_uq"), width = 6.5, height = 2.5)

# =============================================================================
# 2. Figure 19: the 15 most probable RPDAGs of each analysis
# =============================================================================
# Edges are stored with node numbers, "1->5" for a directed edge and "1-5"
# for an undirected edge. A reversal is stored "a<->b": the graph has a->b
# and D_hat has b->a. Changes of the non-Gaussian set are stored as node
# numbers separated by commas.

# Replace a missing value by an empty string.
na0 <- function(x) ifelse(is.na(x), "", as.character(x))

# Write an undirected edge "1-5" as "1 - 5" and separate edges by " | ".
space_undirected <- function(x) {
  gsub("|", " | ", gsub("([0-9])-([0-9])", "\\1 - \\2", x), fixed = TRUE)
}

# Write a reversal "a<->b" as "b->a => a->b", the orientation in D_hat first.
direct_rev <- function(x) {
  if (is.na(x) || !nzchar(x)) return("")
  parts <- strsplit(x, "|", fixed = TRUE)[[1]]
  out <- vapply(parts, function(e) {
    ij <- strsplit(e, "<->", fixed = TRUE)[[1]]
    if (length(ij) != 2L) stop("Unexpected reversal: ", e)
    sprintf("%s->%s => %s->%s", ij[2L], ij[1L], ij[1L], ij[2L])
  }, character(1))
  paste(out, collapse = " | ")
}

# Approximate width of a label as a share of the panel width. An upper-case
# letter or a digit takes W_WIDE, any other character W_NARROW. A label is
# broken between its parts when it would run past the right edge of the
# panel, and it is never wrapped narrower than MIN_ROOM of the panel.
W_WIDE <- 1.03 * 0.0109
W_NARROW <- 1.03 * 0.0056
MIN_ROOM <- 0.35
label_width <- function(s) {
  ch <- strsplit(s, "")[[1]]
  sum(ifelse(grepl("[A-Z0-9]", ch), W_WIDE, W_NARROW))
}
wrap_label <- function(s, room) {
  if (label_width(s) <= room) return(s)
  parts <- strsplit(s, "   ", fixed = TRUE)[[1]]
  lines <- parts[1]
  for (q in parts[-1]) {
    cand <- paste(lines[length(lines)], q, sep = "   ")
    if (label_width(cand) <= room) lines[length(lines)] <- cand else lines <- c(lines, q)
  }
  paste(lines, collapse = "\n")
}

# One figure for each reported analysis.
for (k in seq_along(tags)) {
  uq_dir <- file.path(RD_TAB_DIR, sprintf("uq_%s", tags[k]))
  tab <- read.csv(file.path(uq_dir, "posterior_uncertainty_table.csv"),
                  colClasses = c(D_key = "character"), stringsAsFactors = FALSE)
  labels <- read_tab(file.path(uq_dir, "node_labels.csv"))$symbol
  # Replace node numbers by gene names.
  relabel <- function(x) {
    if (is.na(x) || !nzchar(x)) return(x)
    m <- gregexpr("[0-9]+", x)
    regmatches(x, m) <- lapply(regmatches(x, m), function(v) labels[as.integer(v)])
    x
  }
  tab <- tab[order(tab$rank), ]

  # The label of each graph: the edges added (+e), removed (-e) and reversed
  # (rev), and the genes that leave (nG: -gene) or join (nG: +gene) the
  # non-Gaussian set of D_hat.
  edit_label <- function(i) {
    bits <- character(0)
    if (tab$n_edges_added[i] > 0L) {
      bits <- c(bits, sprintf("+e: %s", relabel(space_undirected(na0(tab$added_edges_str[i])))))
    }
    if (tab$n_edges_removed[i] > 0L) {
      bits <- c(bits, sprintf("-e: %s", relabel(space_undirected(na0(tab$removed_edges_str[i])))))
    }
    if (tab$n_edges_reversed[i] > 0L) {
      bits <- c(bits, sprintf("rev: %s", relabel(direct_rev(tab$reversed_edges_str[i]))))
    }
    if (nzchar(na0(tab$nG_removed_str[i]))) {
      bits <- c(bits, sprintf("nG: -%s", relabel(gsub("[|,]", " | -", na0(tab$nG_removed_str[i])))))
    }
    if (nzchar(na0(tab$nG_added_str[i]))) {
      bits <- c(bits, sprintf("nG: +%s", relabel(gsub("[|,]", " | +", na0(tab$nG_added_str[i])))))
    }
    if (length(bits) == 0L) return("(identical)")
    paste(bits, collapse = "   ")
  }
  tab$edit_label <- vapply(seq_len(nrow(tab)), edit_label, character(1))

  # The graphs shown, with their labels and the vertical band each needs.
  pm <- tab[seq_len(min(RD_TOP_K, nrow(tab))), ]
  seen <- CLASS_LEVELS[CLASS_LEVELS %in% unique(pm$category)]
  room <- pmax(MIN_ROOM, 0.97 - pm$mass / (1.5 * max(pm$mass)))
  pm$bar_label <- mapply(wrap_label, sprintf("%.1f%%   %s", 100 * pm$mass, pm$edit_label),
                         room, USE.NAMES = FALSE)
  pm$category <- factor(pm$category, levels = seen)
  pm$n_lines <- 1L + lengths(regmatches(pm$bar_label, gregexpr("\n", pm$bar_label, fixed = TRUE)))
  band <- pm$n_lines
  pm$pos <- cumsum(band) - band / 2

  g_marg <- ggplot(pm, aes(x = pos, y = mass, fill = category)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.3) +
    geom_text(aes(label = bar_label), hjust = -0.03, vjust = 0.5, size = 4.6,
              lineheight = 0.9, color = "grey20") +
    scale_fill_manual(values = CLASS_COLORS[seen], breaks = seen, labels = seen,
                      name = expression(paste("Class relative to ", hat(D), ":  "))) +
    scale_x_continuous(breaks = pm$pos, labels = as.character(pm$rank),
                       expand = expansion(add = 0.4)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = scales::breaks_pretty(n = 4),
                       expand = expansion(mult = c(0, 0.5))) +
    coord_flip() +
    labs(x = "RPDAG, ranked by posterior mass", y = "Posterior mass") +
    theme_lewis() +
    theme(legend.position = "top", legend.title = element_text(size = rel(0.95)),
          legend.key.width = grid::unit(1.2, "lines"),
          panel.grid.minor.y = element_blank())

  # The height grows with the number of label lines.
  h_marg <- 1.1 + 4.1 * sum(band) / 15
  save_figure(g_marg, file.path(FIG_DIR, paste0("posterior_marginal_annotated_", tags[k])),
              width = 12, height = h_marg)
}

cat(sprintf("Tables read from %s\nFigures written to %s\n", RD_TAB_DIR, FIG_DIR))
