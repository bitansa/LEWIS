# LEWIS 0.1.0

* First release.
* `lewis()` fits LEWIS: a Metropolis-Hastings sampler over DAGs under the
  approximate posterior given by the Laplace approximation of the marginal
  likelihood, followed by the four post-processing steps of the paper
  (estimating the non-Gaussian nodes, projection onto the space of RPDAGs,
  point estimate and credible ball, and the shape of the posterior).
* The defaults are the settings of the paper: the complexity prior with
  `d = 0.03` and `alpha = 0.75`, a start at the estimate of the TL method,
  10,000 max(1, (p/20)^2) iterations, a burn-in of 20 percent, at most 8000
  kept draws (`max_draws = Inf` keeps every draw), the Shapiro-Wilk level
  0.05 / p and the 95 percent credible ball.
* Convergence checks: the effective sample size of the score, and with
  `n_chains > 1` the potential scale reduction factor and the differences
  between chains.
* Building blocks: `lewis_schedule()`, `lewis_score()`, `nongaussian_nodes()`,
  `rpdag()`, `shd()` and `tl_lingam()`. `tl_lingam()` is a copy of the TL
  function of the package TransGraph 1.1.0, which is no longer on CRAN.
* Methods `print()`, `summary()` and `plot()` for fitted objects.
* The data set `tcga_prostate` with the gene expression data of the prostate
  cancer analysis of the paper.
