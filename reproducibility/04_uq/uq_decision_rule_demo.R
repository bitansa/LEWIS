# =============================================================================
# uq_decision_rule_demo.R
#
# The uncertainty quantification study of Section 7.3 and Appendix H. Runs
# LEWIS (complexity) once on the selected dataset and writes the tables
# behind the text, Figure 11(b) and Figure 17.
#
# The dataset is the one selected by the rule of Appendix H. The seed is
# read from results/seed_sweep/seed_sweep_summary.csv, which
# choose_seed_for_demo.R writes. The rule gives seed 2.
#
# Output in results/tables/:
#   decision_rule_summary.csv        one row: acceptance rate, number of
#       distinct DAGs and RPDAGs, edges of D_hat, nG_hat, radius of the
#       credible ball (eps_alpha), mass outside the ball, SHD of D_hat from
#       the true RPDAG, and whether the ball contains the true RPDAG
#   top_rpdag_posterior.csv          every distinct RPDAG with its mass,
#       its SHD from D_hat, its number of edges and its class
#   posterior_uncertainty_table.csv  the 15 most probable RPDAGs with the
#       edges added, removed or reversed relative to D_hat, and the changes
#       of the most frequent non-Gaussian set among their draws
#   posterior_edges.csv              the edges of D_hat with their posterior
#       probabilities (fit$edges)
#   posterior_nodes.csv              the proportion of draws in which each
#       node is declared non-Gaussian (fit$nodes)
#
# The figures are drawn by replot_uq_figures.R from these tables.
#
# Run:   Rscript reproducibility/04_uq/uq_decision_rule_demo.R
#        (from the root of the repository, or from any other folder)
# Time:  about ten seconds on one core.
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
source(file.path(.here, "uq_common.R"))

# ---- The selected seed -------------------------------------------------------
# Apply the rule of Appendix H to the seed sweep: the smallest distance to
# the benchmark median, then the smaller seed.
f_sweep <- file.path(.here, "results", "seed_sweep", "seed_sweep_summary.csv")
if (!file.exists(f_sweep)) {
  stop("File not found: ", f_sweep, ". Run choose_seed_for_demo.R first.")
}
sweep <- read.csv(f_sweep, stringsAsFactors = FALSE)
SEED <- sweep$seed[order(sweep$dist_to_bench_median, sweep$seed)][1]
cat(sprintf("Selected seed: %d\n", SEED))

# ---- Run LEWIS (complexity) ----------------------------------------------------
res <- run_uq(SEED)
fit <- res$fit
s <- res$summary

# ---- Report the numbers of Section 7.3 and Appendix H ---------------------------
cat(sprintf("Schedule: %d iterations, burn-in %d, thinning %d, %d draws\n",
            s$n_iter, s$burn_in, s$thin, s$M))
cat(sprintf("Acceptance rate: %.4f\n", s$acceptance_rate))
cat(sprintf("Distinct DAGs after the burn-in: %d. Distinct RPDAGs: %d\n",
            s$distinct_dags, s$distinct_rpdags))
cat(sprintf("D_hat: %d edges, posterior mass %.4f, SHD from the true RPDAG %d\n",
            s$D_hat_edges, res$top$mass[res$top$is_D_hat], s$shd_Dhat_truth))
cat(sprintf("nG_hat = {%s}. True nG = {%s}\n", s$nG_hat, s$nG_star))
cat(sprintf("Credible ball: radius %d, mass outside %.4f, contains the true RPDAG: %s\n",
            s$eps_alpha, s$mass_outside_ball, s$truth_in_ball_shd))
cat("\nEdges of D_hat:\n")
print(res$edges, row.names = FALSE)
cat("\nNodes:\n")
print(res$nodes, row.names = FALSE)
cat("\nThe 15 most probable RPDAGs:\n")
print(res$uq[, c("rank", "mass", "shd_to_Dhat", "category", "added_edges_str",
                 "removed_edges_str", "reversed_edges_str", "nongaussian")],
      row.names = FALSE)

# ---- Check that the seed sweep gave the same run ------------------------------
row <- sweep[sweep$seed == SEED, ]
if (!isTRUE(row$shd_to_truth == s$shd_Dhat_truth && row$eps_alpha == s$eps_alpha &&
            row$n_dags == s$distinct_dags)) {
  stop("The run differs from the row of seed ", SEED, " in ", f_sweep, ".")
}

# ---- Write the tables ----------------------------------------------------------
write_uq_tables(res, file.path(.here, "results", "tables"))
