#' LEWIS: Consistent DAG Selection for Bayesian Causal Discovery
#'
#' LEWIS (Laplace Error Working model for Inference and Selection of DAG)
#' learns a directed acyclic graph (DAG) from observational data generated
#' by a linear structural equation model. The errors may be Gaussian for
#' some variables and non-Gaussian for others. LEWIS uses a Laplace error
#' working model, a complexity prior on the DAG and a Laplace approximation
#' of the marginal likelihood. Under mild conditions, the posterior
#' probability of the distribution equivalence class of the true DAG tends
#' to one as the sample size grows. A Metropolis-Hastings sampler draws DAGs,
#' and the draws are summarized by a restricted partially directed acyclic
#' graph (RPDAG), a credible ball and posterior probabilities.
#'
#' The main function is [lewis()]. The building blocks [lewis_schedule()],
#' [lewis_score()], [nongaussian_nodes()], [rpdag()], [shd()] and
#' [tl_lingam()] are also exported. The data set [tcga_prostate] contains the gene expression data
#' of the prostate cancer analysis of the paper.
#'
#' @references
#' Chaudhuri, A., Bhattacharya, A., Sarkar, B. and Ni, Y. (2026). Consistent
#' DAG Selection for Bayesian Causal Discovery under Model
#' Misspecification. Manuscript.
#'
#' @keywords internal
#' @useDynLib LEWIS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
