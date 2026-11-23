# =============================================================================
# settings.R
#
# Settings shared by every study: the constants of the complexity prior,
# the level of the Shapiro-Wilk tests, and the names, order and labels of
# the methods.
#
# The values of the prior are those of the paper (Remark "Choice of alpha
# and d" and Appendix G.2). They are also the defaults of lewis().
#
# This file has no dependencies.
# =============================================================================

# Constants of the complexity prior pi(gamma) proportional to
# exp(-n^alpha * d * |gamma|).
LEWIS_ALPHA <- 0.75
LEWIS_D <- 0.03

# Level of the Shapiro-Wilk tests of Step 1, before the Bonferroni division
# by p.
LEWIS_ALPHA_SW <- 0.05

# Method keys in the order used by every table and figure.
METHOD_ORDER <- c("lewis_complexity", "lewis_uniform", "pc", "ges",
                  "lingam", "direct_lingam", "rcd", "tl_zhao",
                  "notears", "pairwise_lingam")

# Display labels of the methods, as in the paper.
METHOD_LABELS <- c(
  lewis_complexity = "LEWIS (complexity)",
  lewis_uniform    = "LEWIS (uniform)",
  pc               = "PC",
  ges              = "GES",
  lingam           = "ICA-LiNGAM",
  direct_lingam    = "DirectLiNGAM",
  rcd              = "RCD",
  tl_zhao          = "TL",
  notears          = "NOTEARS",
  pairwise_lingam  = "Pairwise-LiNGAM"
)

# Convert a vector of method keys to an ordered factor of display labels.
method_factor <- function(x) {
  factor(METHOD_LABELS[x], levels = unname(METHOD_LABELS[METHOD_ORDER]))
}
