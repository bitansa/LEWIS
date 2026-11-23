# common

Code shared by the studies. A script of a study sets `COMMON_DIR` to this
folder and then sources the files it needs.

| File | Contents |
| --- | --- |
| `settings.R` | The constants of the complexity prior (alpha = 0.75, d = 0.03), the level of the Shapiro-Wilk tests (0.05), and the keys, order and labels of the methods. |
| `paths.R` | `lewis_fig_dir()`: the folder of the figures and LaTeX tables. It is the environment variable `LEWIS_FIG_DIR`, or `reproducibility/figures/` if the variable is not set. |
| `scenario_instances.R` | Data generation of the simulations: `stable_seed()`, the ER and SF graphs, the fixed instance of each graph type and p, the view of a setting, and `draw_dataset()`. |
| `paper12_instance.R` | The 12-node DAG of Figure 7(a) with its coefficients and error scales, used by Figure 1. |
| `graph_tools.R` | `to_cpdag()` (CPDAG of a DAG with `pcalg`), `sid_bounds()` (SID bounds with the package `SID`), `n_edges_partial()` and `edge_list()`. |
| `method_wrappers.R` | One wrapper for each method of Section 7.2, `require_method()`, `run_method()`, and `est_partial_dag()`, which maps the output of a method to the graph that is compared with the true RPDAG. |
| `scoring.R` | `true_rpdag_of()` and `score_estimate()`: the SHD from the true RPDAG and the SID bounds against the true DAG. |
| `plot_theme.R` | The style of the figures and `save_figure()`, which writes a PDF and a 300-dpi PNG. |
| `laplace_exact.R`, `laplace_exact.cpp` | The two chains of Appendix F: the chain with the Laplace approximation and the pseudo-marginal chain with importance sampling estimates of the marginal likelihood. |
| `pairwise_lingam/run_pwling_batch.m` | The MATLAB stage of Pairwise-LiNGAM. It runs the authors' function `pwling` on the exported datasets. |
| `record_versions.R` | Writes `versions.txt` with the versions of R, of the R and Python packages, of the C++ compiler and of MATLAB. |
| `check_orientation_conventions.R` | Checks that every wrapper returns its graph with the parent in the row. |

## Conventions

- **Adjacency matrices.** `A[k, j] = 1` means the directed edge k -> j (the
  row is the parent). `A[k, j] = A[j, k] = 1` means the undirected edge
  k - j.
- **Fixed instances.** For each graph type and p, the graph, its
  coefficients and its error scales are drawn once and kept for all
  settings. The datasets of a setting differ only in the errors.
- **Seeds.** Every seed is computed by `stable_seed()` from a readable
  string. The string of a dataset contains a regime tag (for example
  `"mixture"` in the comparative study). Renaming a tag changes the
  datasets.
- **LEWIS.** Every study calls `LEWIS::lewis()` with `center = FALSE`, as in
  the paper. The complexity prior uses `d = 0.03` and `alpha = 0.75`. The
  uniform prior uses `d = 0`. The chain starts at the TL estimate, and the
  number of iterations follows `lewis_schedule(p)`.
- **Failures.** A wrapper that fails prints the message and flags the row.
  The scripts of the studies exit with a non-zero status when a row is
  missing or flagged.

## Before a run on a new machine

```sh
Rscript reproducibility/common/check_orientation_conventions.R
```

It must end with "No wrapper returns reversed edges."

## Versions

To record the software versions of a machine, run, in the same environment
as the jobs:

```sh
export RETICULATE_PYTHON="$(command -v python3)"
Rscript reproducibility/common/record_versions.R
```

This writes `common/versions.txt`. A file name given as the first argument
replaces `versions.txt`.
