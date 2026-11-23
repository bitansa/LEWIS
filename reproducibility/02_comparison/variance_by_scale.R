# =============================================================================
# variance_by_scale.R
#
# How the effect-size scale s changes the variables of the fixed graphs of
# the comparative study. Section 7.2 quotes the result: on SF graphs at
# s = 2, parents and children become nearly collinear, with correlations up
# to 0.999 (column corr_max).
#
# No data are drawn. The covariance of each instance is computed exactly
# from its coefficients and error scales:
#     X = B' X + e,  Sigma = A D A',  A = (I - B')^{-1},  D = diag(sigma^2),
# with B[k, j] the coefficient of k -> j (view$betas, as in draw_dataset()).
#
# Printed for each graph type, p and s:
#   own_share_med, own_share_min  share of the variance of a node that comes
#                                 from its own error (median and smallest
#                                 over the nodes)
#   var_max                       largest variance of a node
#   corr_max                      largest absolute correlation between a
#                                 node and one of its parents
#
# Usage, from any folder, for example from the root of the repository:
#     Rscript reproducibility/02_comparison/variance_by_scale.R
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
suppressPackageStartupMessages(source(file.path(COMMON_DIR, "scenario_instances.R")))

# One row for each graph type, p and s.
rows <- list()
for (gt in c("ER", "SF")) for (p in c(20L, 50L)) {
  inst <- build_instance(gt, p)
  for (s in c(0.25, 0.5, 1, 2)) {
    # Implied covariance of the view.
    v <- instance_view(inst, rho = 0.5, effect_scale = s)
    A <- solve(diag(p) - t(v$betas))
    S <- A %*% diag(v$sigmas^2) %*% t(A)
    # Share of the own error, and correlations along the edges.
    share <- v$sigmas^2 / diag(S)
    R <- cov2cor(S)
    pc <- which(v$dag == 1, arr.ind = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      graph = gt, p = p, s = s, edges = sum(v$dag),
      own_share_med = round(median(share), 3),
      own_share_min = signif(min(share), 3),
      var_max = signif(max(diag(S)), 3),
      corr_max = round(max(abs(R[pc])), 3))
  }
}
# Print the table.
print(do.call(rbind, rows), row.names = FALSE)
