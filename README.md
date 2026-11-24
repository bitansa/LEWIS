
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# LEWIS

<!-- badges: start -->

[![R-CMD-check](https://github.com/bitansa/LEWIS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bitansa/LEWIS/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/LEWIS)](https://CRAN.R-project.org/package=LEWIS)
[![License:
GPL-2](https://img.shields.io/badge/License-GPL--2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
<!-- badges: end -->

LEWIS (Laplace Error Working model for Inference and Selection of DAG)
is a Bayesian method for causal discovery. It learns a directed acyclic
graph (DAG) from observational data that follow a linear structural
equation model, in which the errors of some variables are Gaussian and
the errors of others are not. LEWIS samples DAGs from an approximate
posterior and summarizes the draws by an estimate of the distribution
equivalence class, a credible ball and posterior probabilities. In the
estimate, an edge is directed when the data determine its direction and
undirected otherwise.

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("bitansa/LEWIS")
```

To also build the vignette, use

``` r
remotes::install_github("bitansa/LEWIS", build_vignettes = TRUE)
```

Once the package is on CRAN, install it with:

``` r
install.packages("LEWIS")
```

The package contains C++ and Fortran code. Installing from GitHub needs
compilers for both: the Xcode command line tools and the GNU Fortran
compiler of <https://mac.r-project.org/tools/> on macOS, Rtools on
Windows, and `g++` and `gfortran` on Linux.

## Quick start

All you need is a numeric matrix or data frame `X` with one row per
observation and one column per variable.

``` r
library(LEWIS)
fit <- lewis(X)
fit
```

This runs the sampler and the four post-processing steps with the
settings of the paper. `print(fit)` gives a short description: the data
size, the prior, the schedule of the chain, the effective sample size of
the score, the edges and the non-Gaussian nodes of the point estimate,
the radius of the credible ball and any warnings.

## What LEWIS does

1.  **Model.** Each variable is a linear function of its parents in the
    DAG plus an error. The working model takes the errors to be Laplace.
    The true errors need not be Laplace. Each true error may be Gaussian
    or a Gaussian scale mixture, such as a Student’s t error. The
    coefficients have the prior N(0, 100 I) and the error scales have
    the prior Inv-Gamma(1, 1). The DAG has the complexity prior,
    proportional to exp(-n^alpha d \|gamma\|), where \|gamma\| is the
    number of edges. With `d = 0` this is the uniform prior.
2.  **Sampling.** The marginal likelihood is replaced by its Laplace
    approximation, which gives the score s_n of each DAG. A
    Metropolis-Hastings sampler proposes the addition, deletion or
    reversal of one edge in each iteration.
3.  **Post-processing.** The draws kept after the burn-in are summarized
    in four steps.
    - Step 1 (estimating the non-Gaussian nodes): for each draw, the
      residuals of each node are tested for normality with the
      Shapiro-Wilk test at the level 0.05 / p.
    - Step 2 (projection onto the space of RPDAGs): each draw is mapped
      to the restricted partially directed acyclic graph (RPDAG) of its
      distribution equivalence class.
    - Step 3 (point estimate and credible ball): the point estimate
      `D_hat` is the projected draw with the smallest mean structural
      Hamming distance (SHD) to all projected draws. The credible ball
      holds the graphs within a given SHD of `D_hat`.
    - Step 4 (the shape of the posterior): the distinct RPDAGs are
      listed by posterior mass and compared with `D_hat`.

## The defaults and the paper

The defaults of `lewis()` are the settings of the paper.

| Setting | Default | Argument |
|:---|:---|:---|
| Complexity prior | `d = 0.03`, `alpha = 0.75`. Use `d = 0` for the uniform prior. | `d`, `alpha` |
| Start of the chain | The estimate of the TL method of Zhao, He and Wang (2022), computed by `tl_lingam()`. Use `"empty"` for the empty graph, or give a DAG as a 0/1 matrix with the parents in the rows. | `start` |
| Number of iterations | 10,000 max(1, (p/20)^2) | `n_iter` |
| Burn-in | The first 20 percent of the iterations | `burn_in` |
| Thinning | At most 8000 draws are kept. `max_draws = Inf` keeps every draw after the burn-in. | `max_draws` |
| Shapiro-Wilk level | 0.05 / p (Bonferroni correction over the p nodes) | `alpha_sw` |
| Credible ball | 95 percent | `alpha_ball` |
| Centering | The columns of `X` are centered, because the model has no intercept. | `center` |

The schedule of the paper is 10,000 iterations for up to 20 variables,
62,500 iterations for 50 variables and 250,000 iterations for 100
variables. `lewis_schedule()` shows the iterations, the burn-in, the
thinning interval and the number of kept draws for any number of
variables:

|   p | Iterations | Burn-in | Thinning | Kept draws |
|----:|-----------:|--------:|---------:|-----------:|
|  10 |     10,000 |   2,000 |        1 |      8,000 |
|  50 |     62,500 |  12,500 |        7 |      7,142 |
| 100 |    250,000 |  50,000 |       25 |      8,000 |

``` r
lewis_schedule(50)
```

## Checking convergence

A fit is only useful if the chain has converged. LEWIS helps in three
ways.

- **The effective sample size.** LEWIS computes the effective sample
  size (ESS) of the score of the kept draws and gives a warning when it
  is below `ess_min = 100`. If you see this warning, increase `n_iter`.
- **The trace.** `plot(fit, type = "trace")` shows the score of the kept
  draws. It should show no trend after the burn-in.
- **Several chains.** With `n_chains = 2` (or more), LEWIS runs further
  chains, from the empty graph and from the TL estimate. It reports the
  potential scale reduction factor (PSRF) of the score and of the number
  of edges, and the largest difference between the chains in the
  posterior probabilities of adjacencies and orientations. It gives a
  warning when a PSRF is above 1.05 or a difference is above 0.1. All
  estimates come from chain 1.

``` r
fit <- lewis(X, n_iter = 100000, n_chains = 2, seed = 1)
fit$diagnostics[c("ess_log_score", "psrf_log_score", "psrf_n_edges",
                  "max_diff_adjacency", "max_diff_orientation")]
```

The default number of iterations works well in the simulations of the
paper. Real data can need many more. The real-data analysis of the paper
uses 10^7 iterations.

## What the output contains

`lewis()` returns an object of class `"lewis"`, a list with these
elements.

| Element | Content |
|:---|:---|
| `D_hat` | The point estimate, an RPDAG, as a p x p matrix of 0 and 1 |
| `nG_hat` | The names of the nodes whose errors are declared non-Gaussian |
| `edges` | One row per edge of `D_hat`: `from`, `to`, `type` (directed or undirected), `prob` (the posterior probability that the edge has the same status as in `D_hat`) and `prob_adjacent` (the posterior probability that the two nodes are adjacent) |
| `nodes` | One row per node: `prob_nongaussian` (the proportion of draws in which its error is declared non-Gaussian) and `in_nG_hat` |
| `graphs` | One row per distinct RPDAG, in decreasing order of posterior mass: its `mass`, its `class` (identical, upper, lower or rearranged), its SHD to `D_hat`, the edges `added`, `removed` and `reversed` relative to `D_hat`, its most frequent set of non-Gaussian nodes and whether it is in the credible ball |
| `ball` | The credible ball: its `radius` (an SHD), `alpha_ball` and the posterior mass outside it |
| `trace` | One row per kept draw: the iteration, the score (`log_score`) and the number of edges of the DAG |
| `diagnostics` | The acceptance rate, the schedule, the number of distinct DAGs and RPDAGs, the ESS of the score, the start used, the run time and the warnings. With several chains also the PSRF and the differences between chains. |
| `settings`, `call` | The arguments used and the call |

`summary(fit)` prints the tables `edges`, `nodes` and the 10 most
probable graphs, the posterior mass of each class and the diagnostics.

### How to read an RPDAG matrix

An RPDAG is a graph with directed and undirected edges. It is stored as
a p x p matrix `D` of 0 and 1, with the parent in the row.

- `D[i, j] = 1` and `D[j, i] = 0`: the directed edge i -\> j.
- `D[i, j] = 1` and `D[j, i] = 1`: the undirected edge i - j. The data
  do not determine its direction.
- `D[i, j] = 0` and `D[j, i] = 0`: no edge between i and j.

A directed edge has this direction in every DAG of the distribution
equivalence class. An undirected edge has both directions within the
class. For example, take the DAG a -\> b -\> c. If all errors are
Gaussian, `rpdag()` returns its CPDAG, in which both edges are
undirected. If the error of `a` is non-Gaussian, both edges are
directed.

``` r
library(LEWIS)
A <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
A["a", "b"] <- 1
A["b", "c"] <- 1
rpdag(A)            # all errors Gaussian: a - b - c
#>   a b c
#> a 0 1 0
#> b 1 0 1
#> c 0 1 0
rpdag(A, nG = "a")  # error of a non-Gaussian: a -> b -> c
#>   a b c
#> a 0 1 0
#> b 0 0 1
#> c 0 0 0
```

The classes of Step 4 compare a graph `D` with `D_hat`. `D` is identical
if it equals `D_hat`. Otherwise, `D` is upper if every edge of `D_hat`
is an edge of `D` with the same status. `D` is lower if every edge of
`D` is an edge of `D_hat` with the same status. In all other cases `D`
is rearranged. In the columns `added`, `removed` and `reversed`, a
directed edge is written `a->b` and an undirected edge `a - b`. A
reversed edge is written as its orientation in `D_hat`, then `=>`, then
its orientation in `D`. An edge that changes between directed and
undirected is listed as one removal and one addition. Several edges are
separated by `|`.

## Plots

``` r
plot(fit)                  # the point estimate
plot(fit, type = "ball")   # posterior mass by SHD from the point estimate
plot(fit, type = "trace")  # the score of the kept draws
plot(fit, type = "edges")  # posterior probabilities of the edges
```

- `type = "graph"` draws `D_hat` with the nodes on a circle. Nodes whose
  errors are declared non-Gaussian are black. An edge is solid if its
  posterior probability is at least 0.9 and dashed otherwise.
- `type = "ball"` draws the posterior mass by the SHD from `D_hat`,
  split by class. The dashed line marks the radius of the credible ball.
- `type = "trace"` draws the score of the kept draws against the
  iteration.
- `type = "edges"` draws the posterior probability of each edge of
  `D_hat` as a bar and its adjacency probability as a point.

## A simulated example

We simulate 500 observations from the DAG x1 -\> x2 -\> x3 -\> x4. The
errors of x1 and x3 are Laplace and the errors of x2 and x4 are
Gaussian.

``` r
library(LEWIS)
set.seed(1)
n <- 500
laplace <- function(n) (rexp(n) - rexp(n)) / sqrt(2)
x1 <- laplace(n)
x2 <- 0.8 * x1 + rnorm(n)
x3 <- 0.8 * x2 + laplace(n)
x4 <- 0.8 * x3 + rnorm(n)
X <- cbind(x1, x2, x3, x4)

fit <- lewis(X, seed = 1)
fit
fit$D_hat
fit$edges
summary(fit)
plot(fit)
```

Every edge of this DAG has an endpoint with a non-Gaussian error (x1 or
x3). So the true RPDAG is the DAG itself, with all three edges directed.
Compare the estimate with it by the SHD:

``` r
dag <- matrix(0, 4, 4, dimnames = list(colnames(X), colnames(X)))
dag["x1", "x2"] <- 1
dag["x2", "x3"] <- 1
dag["x3", "x4"] <- 1
truth <- rpdag(dag, nG = c("x1", "x3"))
shd(fit$D_hat, truth)
```

## A real-data example

The package contains the data of the prostate cancer analysis of the
paper, `tcga_prostate`. It holds the expression of 8 genes of the
androgen-receptor program (AR, FOXA1, HOXB13, KLK3, KLK2, TMPRSS2,
NKX3-1 and FKBP5) in 462 primary prostate tumors from The Cancer Genome
Atlas. The data are already adjusted for covariates, centered and
scaled. See `?tcga_prostate` for the preprocessing.

``` r
dim(tcga_prostate)

# A short chain, for illustration only.
fit <- lewis(tcga_prostate, n_iter = 20000, center = FALSE, seed = 1)
fit
fit$nodes
plot(fit)
```

A short chain like this one may give a warning about the effective
sample size. The paper uses 10^7 iterations (`n_iter = 1e7`) and checks
convergence with several chains, because shorter chains started from
different graphs did not agree on these data. One chain of this length
takes about 13 minutes on a laptop. The call below runs two chains.

``` r
fit <- lewis(tcga_prostate, n_iter = 1e7, center = FALSE, n_chains = 2,
             seed = 1)
```

## Building blocks

| Function | Purpose |
|:---|:---|
| `lewis()` | Fit LEWIS: the sampler and the four post-processing steps |
| `lewis_schedule()` | The iterations, burn-in and thinning for `p` variables |
| `lewis_score()` | The score s_n of one DAG |
| `nongaussian_nodes()` | Step 1: the nodes of a DAG with non-Gaussian errors |
| `rpdag()` | Step 2: the RPDAG of a DAG and a set of non-Gaussian nodes |
| `shd()` | The structural Hamming distance between two graphs (Step 3) |
| `tl_lingam()` | The TL method, used for the default start of the chain |
| `print()`, `summary()`, `plot()` | Methods for fitted objects |

The vignette gives a longer introduction: `vignette("LEWIS")`.

## Reproducing the paper

The folder `reproducibility/` of the GitHub repository contains the
scripts that reproduce the simulations, the real-data analysis, the
figures and the tables of the paper. Its
[README](https://github.com/bitansa/LEWIS/blob/main/reproducibility/README.md)
explains what each script reproduces and the order in which to run them.
The scripts call `lewis()` with `center = FALSE`. The folder is not part
of the R package. The results of each dataset of the simulations, from
which the scripts regenerate the figures and tables without rerunning
the simulations, are archived on Zenodo,
<https://doi.org/10.5281/zenodo.24514618>.

## Citation

If you use LEWIS, please cite the paper:

Chaudhuri, A., Bhattacharya, A., Sarkar, B. and Ni, Y. (2026).
Consistent DAG Selection for Bayesian Causal Discovery under Model
Misspecification. Manuscript.

In R, `citation("LEWIS")` gives this reference.

## License

GPL-2. The package contains copies of the Fortran routines `rqbr` and
`rqfnb` of the R package quantreg (Roger Koenker, GPL (\>= 2)), of the
Shapiro-Wilk algorithm of R (R Core Team, GPL (\>= 2)), and of the TL
functions of the R package TransGraph 1.1.0 (Mingyang Ren, Ruixuan Zhao,
Xin He and Junhui Wang, GPL-2).
