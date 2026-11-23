# =============================================================================
# paper_numbers.R
#
# Prints the summary numbers that the paper quotes in its text, in its
# figure captions and in its typed tables. Each number is computed from the
# saved results of the studies. Nothing is run and nothing is written.
#
#   1. Comparative study (02_comparison): Section 7.2 and Appendix G.3.
#   2. Computing time of LEWIS (02_comparison): Appendix E.1.
#   3. Calibration of the complexity prior (03_calibration): the remark
#      "Choice of alpha and d" in Section 5 and Appendix G.2.
#   4. Figure 1 (01_motivation): the mean SHD at both ends of each row, as a
#      check of the figure. The text quotes no number of Figure 1.
#   5. Accuracy of the Laplace approximation (06_laplace_vs_exact):
#      Appendix F and Table 1.
#   6. Uncertainty quantification (04_uq): Section 7.3, Figure 11 and
#      Appendix H.
#   7. Real-data analysis (05_real_data): Section 7.4, Appendix I, Table 2,
#      Figure 18 and Table S15 of the online appendix.
#
# INPUTS. Most numbers come from the small summary files in the results/
# folders of this repository. A few numbers need the results of each
# dataset, which are in the Zenodo archive (see ZENODO.md):
#     02_comparison/results/results_<method>.csv   the paired standard
#         errors, the empty graphs of RCD and the computing times of the
#         existing methods (Section 7.2 and Appendix G.3)
#     04_uq/results/tables/samples_per_dag.csv     the mass outside the
#         credible ball (Section 7.3), if decision_rule_summary.csv does not
#         hold it
# Each block says where its input is. When a file is missing, the block
# prints a message that names the file and skips its numbers.
#
# The correlation "up to 0.999" of Section 7.2 is not computed here. It is
# printed by 02_comparison/variance_by_scale.R (column corr_max).
#
# ROUNDING. Values are printed with more digits than the paper, rounded
# half up. The paper rounds half up. All counts of settings use the
# unrounded means, as in the paper.
#
# Run from the root of the repository:
#     Rscript reproducibility/paper_numbers.R
# The script finds its own folder, so it also runs from any other folder.
# =============================================================================

# ---- Folder of this script --------------------------------------------------
# With Rscript, the folder comes from the --file argument. Otherwise, the
# working directory is used if it is the reproducibility folder or the root
# of the repository.
.here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)),
             fixed = TRUE)
  if (length(fa) > 0L && nzchar(fa[1])) {
    normalizePath(dirname(fa[1]), mustWork = FALSE)
  } else if (dir.exists(file.path(getwd(), "02_comparison"))) {
    getwd()
  } else {
    normalizePath(file.path(getwd(), "reproducibility"), mustWork = FALSE)
  }
})
if (!dir.exists(file.path(.here, "02_comparison"))) {
  stop("Cannot find the reproducibility folder. Run this script with ",
       "Rscript reproducibility/paper_numbers.R from the root of the repository.")
}

# ---- Helpers ------------------------------------------------------------------
LEWIS_C <- "lewis_complexity"
LEWIS_U <- "lewis_uniform"

# Print a section header.
header <- function(title) {
  cat(sprintf("\n%s\n %s\n%s\n", strrep("=", 76), title, strrep("=", 76)))
}

# Print where the input of a block is.
source_note <- function(text) cat(sprintf("  Input: %s\n", text))

# Print one labeled value.
line <- function(label, value) cat(sprintf("  %-62s %s\n", label, value))

# Path of a file of this folder.
path <- function(...) file.path(.here, ...)

# Return TRUE if the file exists. Otherwise print a message and return
# FALSE. Files of the Zenodo archive get a pointer to ZENODO.md.
need <- function(f, zenodo = FALSE) {
  if (file.exists(f)) return(TRUE)
  # The path relative to the reproducibility folder.
  rel <- if (startsWith(f, paste0(.here, "/"))) substring(f, nchar(.here) + 2L) else f
  if (zenodo) {
    cat(sprintf("  Skipped: %s is missing. It is in the Zenodo archive (see ZENODO.md).\n", rel))
  } else {
    cat(sprintf("  Skipped: %s is missing. It is part of this repository.\n", rel))
  }
  FALSE
}

# Read a CSV file.
rd <- function(f) read.csv(f, stringsAsFactors = FALSE)

# Format numbers with k decimals, rounded half up as in the paper.
fmt <- function(x, k = 2) {
  sprintf(paste0("%.", k, "f"), sign(x) * floor(abs(x) * 10^k + 0.5 + 1e-9) / 10^k)
}

# Format a range of numbers.
rng <- function(x, k = 2) paste(fmt(min(x), k), "to", fmt(max(x), k))

# Standard error of the mean of the paired differences x - y.
paired_se <- function(x, y) sd(x - y) / sqrt(length(x))

# =============================================================================
# 1. Comparative study
# =============================================================================
header("1. Comparative study (Section 7.2, Figures 9 and 10, Appendix G.3)")
source_note("02_comparison/results/comparison_cell_summary.csv (this folder)")
f_sum <- path("02_comparison", "results", "comparison_cell_summary.csv")
if (need(f_sum)) {
  s <- rd(f_sum)
  # One key per setting.
  key <- c("experiment", "graph_type", "p", "n", "rho", "effect_scale")
  s$cell <- do.call(paste, c(s[key], sep = "|"))
  cells <- unique(s$cell)
  methods <- sort(unique(s$method))
  compet <- setdiff(methods, c(LEWIS_C, LEWIS_U))
  stopifnot(length(cells) == 54L, length(methods) == 10L, all(s$n_reps == 100L))

  # Matrix of the mean of one column: one row per setting, one column per
  # method.
  wide <- function(col) {
    W <- matrix(NA_real_, length(cells), length(methods), dimnames = list(cells, methods))
    W[cbind(match(s$cell, cells), match(s$method, methods))] <- s[[col]]
    W
  }
  W <- wide("shd_mean")
  # Ranks of the methods in each setting (rank 1 is the smallest mean SHD).
  R <- t(apply(W, 1, rank, ties.method = "average"))
  info <- s[match(cells, s$cell), key]
  fam <- info$experiment
  gt <- info$graph_type
  # The best method and the best existing method of each setting.
  best <- colnames(W)[apply(W, 1, which.min)]
  best_ex_val <- apply(W[, compet], 1, min)
  best_ex <- compet[apply(W[, compet], 1, which.min)]
  is_lewis <- function(m) m %in% c(LEWIS_C, LEWIS_U)

  # ---- Wins and ranks ---------------------------------------------------------
  cat(" Wins and ranks (mean SHD)\n")
  line("settings: core / effect_scale / mixture",
       paste(table(fam)[c("core", "effect_scale", "mixture")], collapse = " / "))
  # The 26 settings with rho = 0.5 and s = 1 (experiment "core").
  line("rho = 0.5, s = 1: a version of LEWIS is best / complexity / uniform",
       sprintf("%d / %d / %d of %d", sum(is_lewis(best[fam == "core"])),
               sum(best[fam == "core"] == LEWIS_C), sum(best[fam == "core"] == LEWIS_U),
               sum(fam == "core")))
  line("rho = 0.5, s = 1: largest p where LEWIS (uniform) is best",
       max(info$p[fam == "core" & best == LEWIS_U]))
  line("rho != 0.5: a version of LEWIS is best (settings)",
       sprintf("%d of %d", sum(is_lewis(best[fam == "mixture"])), sum(fam == "mixture")))
  line("s != 1: a version of LEWIS is best (settings)",
       sprintf("%d of %d", sum(is_lewis(best[fam == "effect_scale"])),
               sum(fam == "effect_scale")))
  line("all: LEWIS (complexity) has the smallest mean SHD", sum(best == LEWIS_C))
  line("all: LEWIS (complexity) is among the best three", sum(R[, LEWIS_C] <= 3))
  mr <- sort(colMeans(R))
  cat("  mean rank over the 54 settings:\n")
  for (m in names(mr)) line(sprintf("    %s", m), fmt(mr[m], 3))
  cat("  settings in which each existing method is among the best three:\n")
  for (m in compet) line(sprintf("    %s", m), sum(R[, m] <= 3))

  # ---- Examples quoted in Section 7.2 ----------------------------------------
  cat(" Mean SHD quoted in Section 7.2\n")
  # Mean SHD of one method in one setting.
  pick <- function(e, g, p, n, rho, sc, m) {
    W[paste(e, g, p, n, format(rho), format(sc), sep = "|"), m]
  }
  # LEWIS (complexity) and the best existing method at rho = 0.5, s = 1.
  core_row <- function(g, p, n) {
    cl <- paste("core", g, p, n, "0.5", "1", sep = "|")
    sprintf("LEWIS (complexity) %s, best existing %s %s",
            fmt(W[cl, LEWIS_C]), best_ex[match(cl, cells)], fmt(best_ex_val[match(cl, cells)]))
  }
  line("ER, p = 20, n = 1000", core_row("ER", 20, 1000))
  line("ER, p = 100, n = 1000", core_row("ER", 100, 1000))
  line("SF, p = 100, n = 1000", core_row("SF", 100, 1000))
  line("ER, p = 50, rho = 1: LEWIS (complexity)",
       fmt(pick("mixture", "ER", 50, 500, 1, 1, LEWIS_C)))
  line("ER, p = 50, rho = 0: LEWIS (complexity), PC, LEWIS (uniform)",
       paste(fmt(c(pick("mixture", "ER", 50, 500, 0, 1, LEWIS_C),
                   pick("mixture", "ER", 50, 500, 0, 1, "pc"),
                   pick("mixture", "ER", 50, 500, 0, 1, LEWIS_U))), collapse = ", "))

  # ---- Figure 10 ----------------------------------------------------------------
  # The settings of Figure 10: n = 500, s = 1, p in {20, 50}, all rho.
  cat(" Figure 10 (n = 500, s = 1, p in {20, 50})\n")
  fig_rho <- info$n == 500 & info$effect_scale == 1 & info$p %in% c(20, 50) &
    fam %in% c("core", "mixture")
  ratio_ex <- best_ex_val / W[, LEWIS_C]
  gap <- best_ex_val - W[, LEWIS_C]
  for (r in sort(unique(info$rho[fig_rho]))) {
    idx <- fig_rho & info$rho == r
    line(sprintf("rho = %g: panels in which a version of LEWIS is best", r),
         sprintf("%d of %d", sum(is_lewis(best[idx])), sum(idx)))
    line(sprintf("rho = %g: best existing over LEWIS (complexity)", r), rng(ratio_ex[idx], 3))
    line(sprintf("rho = %g: best existing minus LEWIS (complexity)", r), rng(gap[idx]))
  }
  line("SF, p = 50, rho = 0.5: best existing minus LEWIS (complexity)",
       fmt(gap[match(paste("core", "SF", 50, 500, "0.5", "1", sep = "|"), cells)]))

  # ---- Against the best existing method ---------------------------------------
  cat(" LEWIS (complexity) against the existing methods\n")
  smaller <- W[, LEWIS_C] < best_ex_val
  line("smaller mean SHD than the best existing method (settings)", sum(smaller))
  for (m in compet) {
    line(sprintf("  smaller mean SHD than %s (settings)", m), sum(W[, LEWIS_C] < W[, m]))
  }
  # The 20 settings that the calibration does not use: rho in {0.25, 0.75}
  # or s != 1.
  uncal <- (fam == "mixture" & info$rho %in% c(0.25, 0.75)) | fam == "effect_scale"
  line("not in the calibration: settings", sum(uncal))
  line("not in the calibration: LEWIS (complexity) best, smaller than best existing",
       sprintf("%d, %d", sum(best[uncal] == LEWIS_C), sum(smaller[uncal])))

  # ---- The two priors -----------------------------------------------------------
  cat(" The two priors (Appendix G.3)\n")
  u_ahead <- W[, LEWIS_U] < W[, LEWIS_C]
  line("LEWIS (uniform) smaller mean SHD than LEWIS (complexity)", sum(u_ahead))
  line("  largest p of those settings", max(info$p[u_ahead]))
  line("  largest difference (edges)", fmt(max((W[, LEWIS_C] - W[, LEWIS_U])[u_ahead])))
  # All 12 settings on ER graphs with p >= 50, whatever rho and s.
  idx <- gt == "ER" & info$p >= 50
  line(sprintf("ER, p >= 50 (%d settings): LEWIS (uniform) over LEWIS (complexity)", sum(idx)),
       rng(W[idx, LEWIS_U] / W[idx, LEWIS_C]))

  # ---- SID ----------------------------------------------------------------------
  cat(" SID (Appendix G.3)\n")
  for (col in c("sid_lo_mean", "sid_up_mean")) {
    S <- wide(col)
    bs <- colnames(S)[apply(S, 1, which.min)]
    line(sprintf("%s: a version of LEWIS is smallest (settings)", col), sum(is_lewis(bs)))
    line(sprintf("%s: uniform smaller than complexity, SF and ER", col),
         sprintf("%d of %d, %d of %d",
                 sum(S[gt == "SF", LEWIS_U] < S[gt == "SF", LEWIS_C]), sum(gt == "SF"),
                 sum(S[gt == "ER", LEWIS_U] < S[gt == "ER", LEWIS_C]), sum(gt == "ER")))
  }

  # ---- Numbers that need the results of each dataset ----------------------------
  cat(" Numbers from the results of each dataset\n")
  source_note("02_comparison/results/results_<method>.csv (Zenodo archive)")
  res_files <- path("02_comparison", "results", sprintf("results_%s.csv", methods))
  if (all(vapply(res_files, need, logical(1), zenodo = TRUE))) {
    a <- do.call(rbind, lapply(res_files, rd))
    a$cell <- do.call(paste, c(a[key], sep = "|"))
    a <- a[order(a$cell, a$method, a$rep), ]
    # SHD of one method in one setting, ordered by dataset.
    shd_of <- function(cl, m) a$shd_rpdag[a$cell == cl & a$method == m]
    # Paired comparison with the best existing method of each setting.
    by2se <- vapply(seq_along(cells), function(i) {
      x <- shd_of(cells[i], LEWIS_C)
      y <- shd_of(cells[i], best_ex[i])
      stopifnot(length(x) == 100L, length(y) == 100L)
      mean(x - y) < -2 * paired_se(x, y)
    }, logical(1))
    line("smaller than the best existing method by more than 2 SE", sum(by2se))
    # RCD on the ER graph with p = 10: datasets with the empty graph.
    for (cl in cells[info$graph_type == "ER" & info$p == 10]) {
      z <- a[a$cell == cl & a$method == "rcd", ]
      line(sprintf("RCD, %s: empty graphs, mean SHD", cl),
           sprintf("%d of %d, %s", sum(z$n_edges_est == 0), nrow(z), fmt(mean(z$shd_rpdag))))
    }
    # Median computing time at (p, n) = (100, 1000) with rho = 0.5, s = 1.
    cat("  median computing time at (p, n) = (100, 1000), ER and SF, in seconds and minutes:\n")
    for (m in methods) {
      md <- vapply(c("ER", "SF"), function(g) {
        median(a$runtime_s[a$method == m & a$experiment == "core" & a$graph_type == g &
                             a$p == 100 & a$n == 1000])
      }, numeric(1))
      line(sprintf("    %s", m),
           sprintf("%s s, %s s (%s min, %s min)", fmt(md[1], 1), fmt(md[2], 1),
                   fmt(md[1] / 60), fmt(md[2] / 60)))
    }
    cat("  The run times of LEWIS in these files are not those of Appendix E.1.\n",
        "  Appendix E.1 uses lewis_runtime_complexity.csv (block 2).\n", sep = "")
  }
}

# =============================================================================
# 2. Computing time of LEWIS
# =============================================================================
header("2. Computing time of LEWIS (Appendix E.1)")
source_note("02_comparison/results/lewis_runtime_complexity.csv (this folder)")
f <- path("02_comparison", "results", "lewis_runtime_complexity.csv")
if (need(f)) {
  rt <- rd(f)
  # Medians over the datasets of each graph type and (p, n). The column
  # ms_per_iter is net of the TL initialization:
  # 1000 * (runtime_s - warm_start_s) / n_iter.
  agg <- aggregate(cbind(runtime_s, warm_start_s, warm_start_frac, ms_per_iter, n_iter) ~
                     graph + p + n, data = rt, FUN = median)
  agg <- agg[order(agg$p, agg$n, agg$graph), ]
  cat("  medians over the datasets:\n")
  cat(sprintf("  %-3s %4s %5s %9s %11s %9s %10s %9s\n", "", "p", "n", "iter",
              "total (s)", "TL (s)", "TL share", "ms/iter"))
  for (k in seq_len(nrow(agg))) {
    cat(sprintf("  %-3s %4d %5d %9d %11s %9s %9s%% %9s\n",
                agg$graph[k], as.integer(agg$p[k]), as.integer(agg$n[k]),
                as.integer(agg$n_iter[k]), fmt(agg$runtime_s[k]), fmt(agg$warm_start_s[k]),
                fmt(100 * agg$warm_start_frac[k], 0), fmt(agg$ms_per_iter[k], 3)))
  }
}

# =============================================================================
# 3. Calibration of the complexity prior
# =============================================================================
header("3. Calibration of the complexity prior (remark in Section 5, Appendix G.2)")
source_note("03_calibration/results/ (this folder)")
for (nm in c("d", "alpha")) {
  f <- path("03_calibration", "results", sprintf("%s_sweep_regret.csv", nm))
  if (need(f)) {
    r <- rd(f)
    # Mean regret in percent, the number of settings in which the value is
    # best, and the number in which its regret is at most 10 percent.
    cat(sprintf("  %s: mean regret (%%), settings in which best, settings within 10%% (of %d)\n",
                nm, r$n_scenarios[1]))
    for (k in seq_len(nrow(r))) {
      line(sprintf("    %s = %g", nm, r[[nm]][k]),
           sprintf("%8s   %3d   %3d", fmt(100 * r$mean_regret[k], 1), r$best_in[k],
                   r$within_10pct[k]))
    }
  }
}
f <- path("03_calibration", "results", "d_sweep_best_d_per_scenario.csv")
if (need(f)) {
  b <- rd(f)
  # The best value of d by the fraction of non-Gaussian errors.
  line("rho = 0: settings with best d >= 0.03",
       sprintf("%d of %d", sum(b$best_d[b$rho == 0] >= 0.03), sum(b$rho == 0)))
  line("rho = 1: settings with best d <= 0.03",
       sprintf("%d of %d", sum(b$best_d[b$rho == 1] <= 0.03), sum(b$rho == 1)))
}

# =============================================================================
# 4. Figure 1
# =============================================================================
header("4. Figure 1, n = 500 (check only, the text quotes no number)")
source_note("01_motivation/results/motivation_cell_summary.csv (this folder)")
f <- path("01_motivation", "results", "motivation_cell_summary.csv")
if (need(f)) {
  m <- rd(f)
  # The regimes of each row, from k = 12 (all t5) to k = 0.
  rows <- list(
    "top row (t5 and Gaussian)" = c("all_t5", sprintf("count_t5g_k%02d", 11:1), "gaussian"),
    "bottom row (t5 and t15)" = c("all_t5", sprintf("count_t5t15_k%02d", 11:1), "t15"))
  for (r in names(rows)) {
    z <- m[m$regime %in% rows[[r]] & m$n == 500, ]
    z$k <- (12:0)[match(z$regime, rows[[r]])]
    cat(sprintf("  %s: mean SHD at k = 12, at k = 0, largest over k\n", r))
    for (meth in intersect(c(LEWIS_C, LEWIS_U, sort(unique(z$method))), unique(z$method))) {
      zz <- z[z$method == meth, ]
      line(sprintf("    %s", meth),
           paste(fmt(c(zz$res_mean[zz$k == 12], zz$res_mean[zz$k == 0], max(zz$res_mean))),
                 collapse = "  "))
    }
  }
}

# =============================================================================
# 5. Accuracy of the Laplace approximation
# =============================================================================
header("5. Laplace approximation against the exact marginal likelihood (Appendix F, Table 1)")
source_note("06_laplace_vs_exact/results/ (this folder)")
f1 <- path("06_laplace_vs_exact", "results", "exact_vs_laplace_summary.csv")
f2 <- path("06_laplace_vs_exact", "results", "exact_vs_laplace_runtime_ratio.csv")
if (need(f1) && need(f2)) {
  sm <- rd(f1)
  ra <- rd(f2)
  # Table 1. The time ratio is the mean computing time of the exact chain
  # over the mean computing time of the Laplace chain.
  cat("  Table 1: p, n, then for each prior SHD (L), SHD (E), time ratio (E over L)\n")
  pn <- unique(sm[, c("p", "n")])
  pn <- pn[order(pn$p, pn$n), ]
  val <- function(pr, meth, p, n) {
    sm$shd_cpdag_mean[sm$prior == pr & sm$method == meth & sm$p == p & sm$n == n]
  }
  rat <- function(pr, p, n) {
    ra$runtime_ratio_exact_over_laplace[ra$prior == pr & ra$p == p & ra$n == n]
  }
  for (k in seq_len(nrow(pn))) {
    p <- pn$p[k]
    n <- pn$n[k]
    cat(sprintf("  %4d %5d | %7s %7s %5s | %7s %7s %5s\n", as.integer(p), as.integer(n),
                fmt(val("complexity", "Laplace", p, n)), fmt(val("complexity", "Exact", p, n)),
                fmt(rat("complexity", p, n), 1), fmt(val("uniform", "Laplace", p, n)),
                fmt(val("uniform", "Exact", p, n)), fmt(rat("uniform", p, n), 1)))
  }
  # Summaries quoted in the text of Appendix F.
  for (pr in c("complexity", "uniform")) {
    L <- sm[sm$prior == pr & sm$method == "Laplace", ]
    E <- sm[sm$prior == pr & sm$method == "Exact", ]
    E <- E[match(paste(L$p, L$n), paste(E$p, E$n)), ]
    d <- L$shd_cpdag_mean - E$shd_cpdag_mean
    line(sprintf("[%s] mean SHD over the 13 pairs, Laplace and exact", pr),
         paste(fmt(c(mean(L$shd_cpdag_mean), mean(E$shd_cpdag_mean))), collapse = ", "))
    line(sprintf("[%s] largest |Laplace - exact| over all pairs", pr), fmt(max(abs(d))))
    line(sprintf("[%s] largest |Laplace - exact|, p <= 50", pr), fmt(max(abs(d[L$p <= 50]))))
    line(sprintf("[%s] Laplace - exact at p = 100", pr), rng(d[L$p == 100]))
  }
  line("time ratio (exact over Laplace), range", rng(ra$runtime_ratio_exact_over_laplace, 1))
}

# =============================================================================
# 6. Uncertainty quantification
# =============================================================================
header("6. Uncertainty quantification on one dataset (Section 7.3, Figure 11, Appendix H)")
source_note("04_uq/results/tables/ and 04_uq/results/seed_sweep/ (this folder)")
d_uq <- path("04_uq", "results")
f1 <- file.path(d_uq, "tables", "decision_rule_summary.csv")
f2 <- file.path(d_uq, "tables", "posterior_uncertainty_table.csv")
f3 <- file.path(d_uq, "tables", "posterior_edges.csv")
f4 <- file.path(d_uq, "tables", "posterior_nodes.csv")
f5 <- file.path(d_uq, "seed_sweep", "seed_sweep_summary.csv")
f6 <- file.path(d_uq, "tables", "samples_per_dag.csv")
if (need(f1) && need(f2) && need(f3) && need(f4)) {
  # The run on the selected dataset.
  ds <- rd(f1)[1, ]
  for (v in c("seed", "acceptance_rate", "distinct_dags", "distinct_rpdags",
              "D_hat_edges", "nG_star", "nG_hat", "eps_alpha", "shd_Dhat_truth",
              "truth_in_ball_shd")) {
    line(v, as.character(ds[[v]]))
  }
  # Mass outside the credible ball. Tables written by uq_decision_rule_demo.R
  # hold it. For the stored tables, it is computed from the draws in the
  # Zenodo archive.
  if ("mass_outside_ball" %in% names(ds)) {
    line("mass outside the credible ball (%)", fmt(100 * ds$mass_outside_ball))
  } else {
    source_note("04_uq/results/tables/samples_per_dag.csv (Zenodo archive)")
    if (need(f6, zenodo = TRUE)) {
      sp <- rd(f6)
      line("mass outside the credible ball (%)", fmt(100 * mean(sp$shd_to_Dhat > ds$eps_alpha)))
    }
  }
  # Posterior mass of the most probable graphs.
  u <- rd(f2)
  line("mass of D_hat (%)", fmt(100 * u$mass[u$is_D_hat][1]))
  line("mass of graphs 2 and 3 (%)", paste(fmt(100 * u$mass[2:3]), collapse = ", "))
  line("mass of graphs 1 to 3 (%)", fmt(100 * sum(u$mass[1:3])))
  line("mass of graph 4 (%)", fmt(100 * u$mass[4]))
  line("largest mass of graphs 4 to 15 (%)", fmt(100 * max(u$mass[4:15])))
  # Figure 11(b): graphs drawn at the same position (SHD from D_hat, number
  # of edges minus that of D_hat).
  pos <- paste0("(", u$shd_to_Dhat, ", ", u$n_edges_added - u$n_edges_removed, ")")
  top <- seq_len(min(15L, nrow(u)))
  grp <- split(u$rank[top], factor(pos[top], levels = unique(pos[top])))
  grp <- grp[lengths(grp) > 1L]
  line("Figure 11(b): graphs at the same position",
       paste(sprintf("%s at %s", vapply(grp, paste, character(1), collapse = ", "), names(grp)),
             collapse = ". "))
  # Edge probabilities of D_hat.
  e <- rd(f3)
  line("edges of D_hat with probability >= 0.99", sum(e$prob >= 0.99))
  low <- e[e$prob < 0.99, ]
  line("the other edges and their probabilities",
       paste(sprintf("%s%s%s %s", low$from, ifelse(low$type == "directed", "->", " - "),
                     low$to, fmt(low$prob, 3)), collapse = ", "))
  # Proportion of draws in which a node is declared non-Gaussian.
  nd <- rd(f4)
  nd$node <- as.character(nd$node)
  line("node 1 declared non-Gaussian (% of draws)",
       fmt(100 * nd$prob_nongaussian[nd$node == "1"]))
  line("node 8 declared non-Gaussian (% of draws)",
       fmt(100 * nd$prob_nongaussian[nd$node == "8"]))
  line("nodes classified the same way in all draws",
       paste(nd$node[nd$prob_nongaussian %in% c(0, 1)], collapse = ", "))
}
if (need(f5)) {
  # The 10 datasets of the selection rule (Appendix H).
  sw <- rd(f5)
  line("10 datasets: SHD of D_hat from the truth", rng(sw$shd_to_truth, 0))
  line("10 datasets: radius of the credible ball", rng(sw$eps_alpha, 0))
  line("10 datasets: the ball contains the true RPDAG",
       sprintf("%d of %d", sum(as.logical(sw$truth_in_ball)), nrow(sw)))
}
f <- path("02_comparison", "results", "comparison_cell_summary.csv")
if (need(f)) {
  # The benchmark of the selection rule: the median SHD of LEWIS
  # (complexity) in the same setting of the comparative study.
  cs <- rd(f)
  z <- cs[cs$method == LEWIS_C & cs$experiment == "core" & cs$graph_type == "ER" &
            cs$p == 10 & cs$n == 500, ]
  line("median SHD of LEWIS (complexity), ER, p = 10, n = 500", z$shd_med)
}

# =============================================================================
# 7. Real-data analysis
# =============================================================================
header("7. Real-data analysis (Section 7.4, Appendix I, Table 2, Figure 18, Table S15)")
source_note("05_real_data/results/tables/ (this folder)")
rd_tab <- path("05_real_data", "results", "tables")
analyses <- c(prostate = "PRAD_ar8", liver = "LIHC_chol8")
# The paper calls an edge between two regulators or two targets "same level".
order_label <- function(x) ifelse(x == "within tier", "same level", x)
f_screen <- file.path(rd_tab, "screen_summary.csv")
for (an in names(analyses)) {
  net <- analyses[[an]]
  tag <- paste0(net, "_it10000k")
  files <- file.path(rd_tab, c(sprintf("summary_%s.csv", tag),
                               sprintf("point_estimate_%s.csv", tag),
                               sprintf("posterior_rpdag_%s.csv", tag),
                               sprintf("convergence_%s.csv", tag),
                               sprintf("stability_%s.csv", tag),
                               sprintf("node_posterior_%s.csv", tag),
                               sprintf("qc_summary_%s.csv", net),
                               file.path(sprintf("uq_%s", tag), "decision_rule_summary.csv")))
  if (!all(vapply(files, need, logical(1)))) next
  sm <- rd(files[1])
  pe <- rd(files[2])
  po <- rd(files[3])
  cv <- rd(files[4])
  sb <- rd(files[5])
  np <- rd(files[6])
  qc <- rd(files[7])
  uq <- rd(files[8])
  qv <- function(k) qc$value[qc$key == k]
  cat(sprintf(" [%s analysis, %s]\n", an, tag))
  # Data and preprocessing.
  line("primary tumors, removed FFPE, removed flagged",
       sprintf("%s, %s, %s", qv("n_primary_tumor"), qv("removed_ffpe"), qv("removed_do_not_use")))
  line("removed by library size", qv("removed_library_size"))
  line("fraction with a known subtype", qv("subtype_known_frac"))
  line("n (patients), p (genes)", sprintf("%d, %d", sm$n, sm$p))
  # Point estimate.
  line("genes declared Gaussian", if (nzchar(sm$gaussian)) sm$gaussian else "none")
  line("edges of D_hat, edges with probability >= 0.9",
       sprintf("%d, %d", nrow(pe), sum(pe$p_status_in_D_hat >= 0.9)))
  line("regulator-target edges that agree, that point against",
       sprintf("%d, %d", sm$agrees, sm$against))
  # Uncertainty.
  line("mass of D_hat (%) and its rank",
       sprintf("%s, %d", fmt(100 * uq$mass_D_hat), po$rank[po$is_D_hat][1]))
  line("radius of the credible ball", uq$eps_alpha)
  line("distinct RPDAGs", uq$n_distinct_rpdags)
  line("mass outside the ball (%)", fmt(100 * uq$mass_outside_ball))
  cl <- tapply(po$mass, po$category, sum)
  line("class mass (%): upper, lower, rearranged",
       paste(fmt(100 * cl[c("upper", "lower", "rearranged")], 1), collapse = ", "))
  line("mass of the two most probable graphs (%)", paste(fmt(100 * po$mass[1:2]), collapse = ", "))
  line("mass beyond SHD 10 from D_hat (%)", fmt(100 * sum(po$mass[po$shd_to_Dhat > 10]), 3))
  # Convergence.
  line("PSRF of the score, of the number of edges",
       paste(fmt(c(max(cv$psrf_log_score), max(cv$psrf_edges_dag)), 4), collapse = ", "))
  line("gap between the highest scores of the chains",
       fmt(max(cv$max_log_score) - min(cv$max_log_score), 3))
  line("largest difference between chains, adjacency and orientation",
       paste(fmt(c(max(cv$max_diff_p_adjacent), max(cv$max_diff_p_directed)), 4),
             collapse = ", "))
  if (file.exists(f_screen)) {
    # The same difference with 10^6 iterations, from the screen.
    sc <- rd(f_screen)
    v6 <- sc$max_pair_diff[sc$network == net & sc$length == "1000k"]
    if (length(v6) == 1L) line("largest difference between chains with 10^6 iterations", fmt(v6, 4))
  }
  line("effective sample size of the score", fmt(uq$ess_log_score, 0))
  line("repetitions with other seeds that return D_hat",
       sprintf("%d of %d", sum(sb$shd_to_reported == 0), nrow(sb)))
  line("computing time of one chain (minutes)", fmt(uq$runtime_s / 60, 1))
  # Table 2: Shapiro-Wilk p-values and non-Gaussian proportions.
  cat("  Table 2: gene, level, Shapiro-Wilk p-value, non-Gaussian proportion\n")
  for (k in seq_len(nrow(np))) {
    cat(sprintf("    %-8s %-10s %9.1e %6s\n", np$gene[k], np$role[k],
                np$sw_p_value[k], fmt(np$p_non_gaussian[k])))
  }
  line("level alpha_sw / p", unique(np$sw_level))
  # Table S15: edges of D_hat.
  cat("  Table S15: edge, P(edge), P(adjacent), order\n")
  for (k in seq_len(nrow(pe))) {
    cat(sprintf("    %-22s %6s %6s  %s\n", pe$edge[k], fmt(pe$p_status_in_D_hat[k], 3),
                fmt(pe$p_adjacent[k], 3), order_label(pe$known_order[k])))
  }
}
if (need(f_screen)) {
  # The screen of the 18 gene sets (Appendix I, "Choice of the gene sets").
  sc <- rd(f_screen)
  cat(" Screen of the 18 gene sets (Appendix I)\n")
  line("gene sets, with 14 genes, with 8 genes",
       sprintf("%d, %d, %d", length(unique(sc$network)),
               length(unique(sc$network[sc$p == 14])), length(unique(sc$network[sc$p == 8]))))
  cat(sprintf("  %-14s %-8s %4s %8s %7s %6s %4s  %s\n", "set", "length", "p", "PSRF",
              "gap", "diff", "rep", "pass"))
  for (k in seq_len(nrow(sc))) {
    cat(sprintf("  %-14s %-8s %4d %8s %7s %6s %4s  %s\n", sc$network[k], sc$length[k],
                as.integer(sc$p[k]), fmt(sc$psrf_log_score[k], 4), fmt(sc$best_score_gap[k]),
                fmt(sc$max_pair_diff[k], 3), as.character(sc$n_reproduced[k]),
                as.character(sc$pass[k])))
  }
  line("sets that pass", paste(unique(sc$network[as.logical(sc$pass)]), collapse = ", "))
}
cat("\n")
