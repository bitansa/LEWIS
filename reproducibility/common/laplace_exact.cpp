// =============================================================================
// laplace_exact.cpp
//
// The exact chain of Appendix F of the paper: the Metropolis-Hastings
// sampler of LEWIS in which the Laplace approximation of each node is
// replaced by an unbiased importance-sampling estimate of the marginal
// likelihood of the node (a pseudo-marginal chain). Its stationary
// distribution is the exact posterior over DAGs.
//
// This file is the part of the compiled engine that produced Table 1 of
// the paper and that the package LEWIS does not contain:
//   * The least absolute deviation fit of a node on its parents, with the
//     Barrodale-Roberts routine rqbr of the R package quantreg (and the
//     fallbacks rqfnb and ridge least squares, as in the runs of the paper).
//   * The importance-sampling estimate of log m_j(D | pa) of one node.
//   * The sampler, with the per-node estimates as the score.
//
// SAME RANDOM NUMBERS AS THE RUNS OF THE PAPER. The random numbers come from
// the generator of R (unif_rand, norm_rand, R_unif_index, Rf_rnorm,
// Rf_rchisq), drawn in the same order as the pure R version of the method:
//   * Proposal: sample.int(3, 1) for the move, sample.int(p, 2) with up to
//     64 tries for an addition, sample.int(|gamma|, 1) for a deletion or a
//     reversal, and log(runif(1)) for the acceptance.
//   * Importance sampling of one node, for each of the M draws: rnorm(1)
//     for log theta_j, then (when the node has parents) rnorm(d_j) and
//     rchisq(1, 5) for the multivariate t draw of mvtnorm::rmvt.
// The linear algebra uses the BLAS and LAPACK routines that R uses for the
// same expressions (dsyrk, dgemv, dgemm, dgesv, dsyevr, dpotrf, dtrsm,
// dgesdd), and sums are accumulated in long double as in R.
//
// The calling R code (laplace_exact.R) compiles this file with
// Rcpp::sourceCpp() and hands over the two Fortran routines of quantreg
// with exact_set_quantreg_symbols().
// =============================================================================

#include <vector>
#include <string>
#include <cmath>
#include <cfloat>
#include <climits>
#include <algorithm>
#include <unordered_map>
#include <cstring>
#include <cstdio>

#ifndef USE_FC_LEN_T
#  define USE_FC_LEN_T
#endif
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Random.h>
#include <Rmath.h>
#ifndef FCONE
#  define FCONE
#endif
#define LEWIS_ERROR(msg) Rcpp::stop(msg)
#define LEWIS_PRINTF Rprintf

// =============================================================================
// 0.  SUMS AND MEANS AS IN R
// =============================================================================
// R accumulates sums in long double. mean() adds a second refinement pass.

static double r_sum(const double* x, int n) {
  long double s = 0.0L;
  for (int i = 0; i < n; ++i) s += x[i];
  return (double) s;
}

static double r_mean(const double* x, int n) {
  long double s = 0.0L;
  for (int i = 0; i < n; ++i) s += x[i];
  s /= n;
  if (std::isfinite((double) s)) {
    long double t = 0.0L;
    for (int i = 0; i < n; ++i) t += (x[i] - s);
    s += t / n;
  }
  return (double) s;
}

// =============================================================================
// 1.  DAG UTILITIES
// =============================================================================
// Adjacency A is p x p, column-major, A[i + p*j] = 1 iff i -> j (row = parent).

static bool is_acyclic_core(const std::vector<int>& A, int p) {
  std::vector<int> indeg(p, 0);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i) if (A[i + p * j]) indeg[j]++;
  std::vector<int> stack; stack.reserve(p);
  for (int j = 0; j < p; ++j) if (indeg[j] == 0) stack.push_back(j);
  int done = 0;
  while (!stack.empty()) {
    int v = stack.back(); stack.pop_back(); ++done;
    for (int j = 0; j < p; ++j)
      if (A[v + p * j] && --indeg[j] == 0) stack.push_back(j);
  }
  return done == p;
}

// which(A[, j] == 1L): parents of j in increasing order (0-based).
static std::vector<int> parents_of(const std::vector<int>& A, int p, int j) {
  std::vector<int> pa;
  for (int i = 0; i < p; ++i) if (A[i + p * j] == 1) pa.push_back(i);
  return pa;
}

// =============================================================================
// 2.  L1 (MEDIAN) REGRESSION  ==  quantreg::rq.fit(x, y, tau = 0.5, "br")
// =============================================================================

typedef void (*rqbr_fn)(int* m, int* nn, int* m5, int* n3, int* n4,
                        double* a, double* b, double* t, double* toler,
                        int* ift, double* x, double* e, int* s, double* wa,
                        double* wb, int* nsol, int* ndsol, double* sol,
                        double* dsol, int* lsol, int* h, double* qn,
                        double* cutoff, double* ci, double* tnmat,
                        double* big, int* lci1);

typedef void (*rqfnb_fn)(int* n, int* p, double* a, double* c, double* rhs,
                         double* d, double* u, double* beta, double* eps,
                         double* wn, double* wp, int* nit, int* info);

// The two entry points are handed over by laplace_exact.R through
// exact_set_quantreg_symbols(). R_FindSymbol is the fallback.
static rqbr_fn  RQBR_PTR  = NULL;
static rqfnb_fn RQFNB_PTR = NULL;

static rqbr_fn  get_rqbr() {
  if (RQBR_PTR == NULL) {
    RQBR_PTR = (rqbr_fn) R_FindSymbol("rqbr", "quantreg", NULL);
    if (RQBR_PTR == NULL) LEWIS_ERROR("laplace_exact: quantreg's Fortran routine 'rqbr' is not "
                                      "available -- loadNamespace('quantreg') before use");
  }
  return RQBR_PTR;
}
static rqfnb_fn get_rqfnb() {
  if (RQFNB_PTR == NULL) {
    RQFNB_PTR = (rqfnb_fn) R_FindSymbol("rqfnb", "quantreg", NULL);
    if (RQFNB_PTR == NULL) LEWIS_ERROR("laplace_exact: quantreg's Fortran routine 'rqfnb' is not "
                                       "available -- loadNamespace('quantreg') before use");
  }
  return RQFNB_PTR;
}

// qr(x)$rank of base R (LINPACK dqrdc2, tol = 1e-7): Householder reduction
// with limited column pivoting. A column whose remaining norm falls below
// tol times its original norm is moved to the end and not counted. Only
// the rank is needed here, so the pivoted columns are skipped.
static int linpack_qr_rank(std::vector<double> x, int n, int p, double tol) {
  std::vector<double> qraux(p), work2(p);
  for (int j = 0; j < p; ++j) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += x[i + n * j] * x[i + n * j];
    qraux[j] = std::sqrt(s);
    work2[j] = qraux[j] == 0.0 ? 1.0 : qraux[j];
  }
  int lup = std::min(n, p), k = p + 1;   // 1-based k as in dqrdc2
  for (int l = 1; l <= lup; ++l) {
    // cycle columns l..p left to right until one with non-negligible norm
    while (!(l >= k || qraux[l - 1] >= work2[l - 1] * tol)) {
      // move column l to the end, shift the others left
      std::vector<double> col(x.begin() + n * (l - 1), x.begin() + n * l);
      double q = qraux[l - 1], w = work2[l - 1];
      for (int j = l; j < p; ++j) {
        std::copy(x.begin() + n * j, x.begin() + n * (j + 1), x.begin() + n * (j - 1));
        qraux[j - 1] = qraux[j]; work2[j - 1] = work2[j];
      }
      std::copy(col.begin(), col.end(), x.begin() + n * (p - 1));
      qraux[p - 1] = q; work2[p - 1] = w;
      --k;
    }
    if (l == n) break;
    double nrmxl = 0.0;
    for (int i = l - 1; i < n; ++i) nrmxl += x[i + n * (l - 1)] * x[i + n * (l - 1)];
    nrmxl = std::sqrt(nrmxl);
    if (nrmxl == 0.0) continue;
    double* xl = &x[(l - 1) + n * (l - 1)];
    if (xl[0] != 0.0) nrmxl = (xl[0] < 0) ? -std::fabs(nrmxl) : std::fabs(nrmxl);
    for (int i = 0; i < n - l + 1; ++i) xl[i] *= 1.0 / nrmxl;
    xl[0] = 1.0 + xl[0];
    for (int j = l + 1; j <= p; ++j) {
      double* xj = &x[(l - 1) + n * (j - 1)];
      double t = 0.0;
      for (int i = 0; i < n - l + 1; ++i) t += xl[i] * xj[i];
      t = -t / xl[0];
      for (int i = 0; i < n - l + 1; ++i) xj[i] += t * xl[i];
      if (qraux[j - 1] == 0.0) continue;
      double tt = 1.0 - (std::fabs(xj[0]) / qraux[j - 1]) * (std::fabs(xj[0]) / qraux[j - 1]);
      tt = std::max(tt, 0.0);
      if (std::fabs(tt) < 1e-6) {
        double s = 0.0;
        for (int i = 1; i < n - l + 1; ++i) s += xj[i] * xj[i];
        qraux[j - 1] = std::sqrt(s);
      } else {
        qraux[j - 1] = qraux[j - 1] * std::sqrt(tt);
      }
    }
    qraux[l - 1] = xl[0];
    xl[0] = -nrmxl;
  }
  return std::min(k - 1, n);
}

struct L1Fit {
  std::vector<double> b;       // coefficients (length d_j)
  double theta;                // mean absolute residual (floored at 1e-12)
  std::vector<double> resid;   // finite residuals
  std::vector<double> resid_raw; // all n residuals (before the finite filter)
};

// The residual y - Xpa %*% b as R computes it (dgemv, then elementwise).
static void residual_dgemv(const std::vector<double>& Xpa, int n, int d,
                           const std::vector<double>& y,
                           const std::vector<double>& b,
                           std::vector<double>& res) {
  std::vector<double> fitted(n, 0.0);
  const char* tN = "N"; const double one = 1.0, zero = 0.0; const int ione = 1;
  F77_CALL(dgemv)(tN, &n, &d, &one, Xpa.data(), &n, b.data(), &ione, &zero,
                  fitted.data(), &ione FCONE);
  res.resize(n);
  for (int i = 0; i < n; ++i) res[i] = y[i] - fitted[i];
}

// crossprod(Xpa): dsyrk on the upper triangle, then mirrored (R's symcrossprod).
static void crossprod_sym(const std::vector<double>& Xpa, int n, int d,
                          std::vector<double>& G) {
  G.assign((size_t) d * d, 0.0);
  const char* up = "U"; const char* tr = "T"; const double one = 1.0, zero = 0.0;
  F77_CALL(dsyrk)(up, tr, &d, &n, &one, Xpa.data(), &n, &zero, G.data(), &d FCONE FCONE);
  for (int i = 1; i < d; ++i)
    for (int j = 0; j < i; ++j) G[i + d * j] = G[j + d * i];
}

// solve(A, B) of base R: dgesv, then the reciprocal-condition check with
// dgecon on the 1-norm. Returns false when R would have raised an error.
static bool r_solve(std::vector<double> A, int n, std::vector<double>& B, int nrhs) {
  std::vector<double> Aorig(A);
  std::vector<int> ipiv(n);
  int info = 0;
  F77_CALL(dgesv)(&n, &nrhs, A.data(), &n, ipiv.data(), B.data(), &n, &info);
  if (info != 0) return false;
  const char* one = "1";
  double anorm = F77_CALL(dlange)(one, &n, &n, Aorig.data(), &n, (double*) NULL FCONE);
  std::vector<double> work(4 * (size_t) n);
  double rcond = 0.0;
  F77_CALL(dgecon)(one, &n, A.data(), &n, &anorm, &rcond, work.data(), ipiv.data(), &info FCONE);
  if (rcond < DBL_EPSILON) return false;
  return true;
}

// Ridge least-squares fallback of the L1 fit:  solve(crossprod(Xpa) + 1e-8 I, crossprod(Xpa, y)).
static std::vector<double> ridge_ols(const std::vector<double>& Xpa, int n, int d,
                                     const std::vector<double>& y) {
  std::vector<double> G; crossprod_sym(Xpa, n, d, G);
  for (int k = 0; k < d; ++k) G[k + d * k] += 1e-8;
  std::vector<double> Xty(d, 0.0);
  const char* tT = "T"; const double one = 1.0, zero = 0.0; const int ione = 1;
  F77_CALL(dgemv)(tT, &n, &d, &one, Xpa.data(), &n, y.data(), &ione, &zero, Xty.data(), &ione FCONE);
  if (!r_solve(G, d, Xty, 1))
    LEWIS_ERROR("laplace_exact: ridge OLS fallback failed (singular Gram matrix)");
  return Xty;
}

// Coefficients of quantreg::rq.fit.br(x, y, tau = 0.5). Returns false where
// the R function would stop with "Singular design matrix".
static bool rq_fit_br(const std::vector<double>& Xpa, int n, int p,
                      const std::vector<double>& y, std::vector<double>& coef) {
  if (linpack_qr_rank(Xpa, n, p, 1e-7) < p) return false;
  double tol = R_pow(DBL_EPSILON, 2.0 / 3.0);
  double big = DBL_MAX;
  int nsol = 2, ndsol = 2;
  int m = n, nn = p, m5 = n + 5, n3 = p + 3, n4 = p + 4;
  double tau = 0.5;
  int ift = 1;
  std::vector<double> a(Xpa), b(y);                 // .Fortran copies its arguments
  std::vector<double> x(p, 0.0), e(n, 0.0);
  std::vector<int> s(n, 0);
  std::vector<double> wa((size_t) (n + 5) * (p + 4), 0.0), wb(n, 0.0);
  std::vector<double> sol((size_t) (p + 3) * nsol, 0.0), dsol((size_t) n * ndsol, 0.0);
  int lsol = 0;
  std::vector<int> h((size_t) p * nsol, 0);
  std::vector<double> qn(p, 0.0);
  double cutoff = 0.0;
  std::vector<double> ci(4 * (size_t) p, 0.0), tnmat(4 * (size_t) p, 0.0);
  int lci1 = 0;
  get_rqbr()(&m, &nn, &m5, &n3, &n4, a.data(), b.data(), &tau, &tol, &ift,
             x.data(), e.data(), s.data(), wa.data(), wb.data(), &nsol, &ndsol,
             sol.data(), dsol.data(), &lsol, h.data(), qn.data(), &cutoff,
             ci.data(), tnmat.data(), &big, &lci1);
  coef = x;
  return true;
}

// quantreg::rq.fit.fnb(x, y, tau = 0.5) coefficients (the "fn" fallback).
static void rq_fit_fnb(const std::vector<double>& Xpa, int n, int p,
                       const std::vector<double>& y, std::vector<double>& coef) {
  double tau = 0.5, beta = 0.99995, eps = 1e-06;
  std::vector<double> a((size_t) n * p);                 // t(x): p x n
  for (int i = 0; i < n; ++i) for (int j = 0; j < p; ++j) a[j + (size_t) p * i] = Xpa[i + (size_t) n * j];
  std::vector<double> c(n); for (int i = 0; i < n; ++i) c[i] = -y[i];
  std::vector<double> rhs(p);
  for (int j = 0; j < p; ++j) rhs[j] = (1 - tau) * r_sum(&Xpa[(size_t) n * j], n);
  std::vector<double> d(n, 1.0), u(n, 1.0), wn(10 * (size_t) n, 0.0);
  for (int i = 0; i < n; ++i) wn[i] = 1 - tau;
  std::vector<double> wp((size_t) (p + 3) * p, 0.0);
  int nit[3] = {0, 0, 0}, info = 0, nn = n, pp = p;
  get_rqfnb()(&nn, &pp, a.data(), c.data(), rhs.data(), d.data(), u.data(),
              &beta, &eps, wn.data(), wp.data(), nit, &info);
  coef.resize(p);
  for (int j = 0; j < p; ++j) coef[j] = -wp[j];
}

// L1 fit of node j on its parents pa: the coefficients b and the scale
// theta (the mean absolute residual). A node without parents has b empty
// and theta equal to the mean absolute value of its data.
static L1Fit fit_node_l1_core(const double* X, int n, int p, int j,
                              const std::vector<int>& pa) {
  L1Fit f;
  int d = (int) pa.size();
  if (d == 0) {
    f.resid.assign(X + (size_t) n * j, X + (size_t) n * (j + 1));
    f.resid_raw = f.resid;
    std::vector<double> ar(n);
    for (int i = 0; i < n; ++i) ar[i] = std::fabs(f.resid[i]);
    f.theta = std::max(r_mean(ar.data(), n), 1e-12);
    return f;
  }
  std::vector<double> y(X + (size_t) n * j, X + (size_t) n * (j + 1));
  std::vector<double> Xpa((size_t) n * d);
  for (int k = 0; k < d; ++k)
    std::copy(X + (size_t) n * pa[k], X + (size_t) n * (pa[k] + 1), Xpa.begin() + (size_t) n * k);

  std::vector<double> coef;
  bool ok = rq_fit_br(Xpa, n, d, y, coef);
  if (!ok) rq_fit_fnb(Xpa, n, d, y, coef);          // the Frisch-Newton fallback
  bool finite = true;
  for (int k = 0; k < d; ++k) if (!std::isfinite(coef[k])) { finite = false; break; }
  if (!finite) coef = ridge_ols(Xpa, n, d, y);
  for (int k = 0; k < d; ++k) if (!std::isfinite(coef[k])) coef[k] = 0.0;

  std::vector<double> res;
  residual_dgemv(Xpa, n, d, y, coef, res);
  f.resid_raw = res;
  f.resid.clear(); f.resid.reserve(n);
  for (int i = 0; i < n; ++i) if (std::isfinite(res[i])) f.resid.push_back(res[i]);
  int nr = (int) f.resid.size();
  if (nr == 0) f.theta = 1e-12;
  else {
    std::vector<double> ar(nr);
    for (int i = 0; i < nr; ++i) ar[i] = std::fabs(f.resid[i]);
    f.theta = std::max(r_mean(ar.data(), nr), 1e-12);
  }
  f.b = coef;
  return f;
}

// =============================================================================
// 3.  IMPORTANCE-SAMPLING ESTIMATE OF log m_j
// =============================================================================
// Priors: b_j ~ N(0, sigma_b2 I) and theta_j ~ Inv-Gamma(ig_alpha, ig_beta).
// Importance distributions: theta_j ~ LogNormal(log theta_hat - c_n / 2, c_n)
// with c_n = log(1 + 1/n), and b_j ~ multivariate t with nu degrees of
// freedom, center b_hat and scale ((nu - 2) / nu) (theta_hat / 2)
// (D_pa^T D_pa)^{-1}. The draws follow mvtnorm::rmvt (eigen square root of
// the scale matrix, rnorm(d) then rchisq(1, nu)) and the density follows
// mvtnorm::dmvt (Cholesky factor and backsolve).

// eigen(sigma, symmetric = TRUE): dsyevr as in La_rs, then R reverses to
// descending order.  Returns values (descending) and vectors (columns).
static void r_eigen_sym(const std::vector<double>& sigma, int n,
                        std::vector<double>& values, std::vector<double>& vectors) {
  std::vector<double> A(sigma);
  const char* jobv = "V"; const char* range = "A"; const char* uplo = "L";
  double vl = 0.0, vu = 0.0, abstol = 0.0;
  int il = 0, iu = 0, m = 0, info = 0;
  std::vector<double> w(n), z((size_t) n * n);
  std::vector<int> isuppz(2 * (size_t) n);
  double tmp = 0.0; int itmp = 0, lwork = -1, liwork = -1;
  F77_CALL(dsyevr)(jobv, range, uplo, &n, A.data(), &n, &vl, &vu, &il, &iu, &abstol,
                   &m, w.data(), z.data(), &n, isuppz.data(), &tmp, &lwork, &itmp,
                   &liwork, &info FCONE FCONE FCONE);
  if (info != 0) LEWIS_ERROR("laplace_exact: dsyevr workspace query failed");
  lwork = (int) tmp; liwork = itmp;
  std::vector<double> work(lwork); std::vector<int> iwork(liwork);
  F77_CALL(dsyevr)(jobv, range, uplo, &n, A.data(), &n, &vl, &vu, &il, &iu, &abstol,
                   &m, w.data(), z.data(), &n, isuppz.data(), work.data(), &lwork,
                   iwork.data(), &liwork, &info FCONE FCONE FCONE);
  if (info != 0) LEWIS_ERROR("laplace_exact: dsyevr failed");
  values.resize(n); vectors.resize((size_t) n * n);
  for (int k = 0; k < n; ++k) {
    values[k] = w[n - 1 - k];
    std::copy(z.begin() + (size_t) n * (n - 1 - k), z.begin() + (size_t) n * (n - k),
              vectors.begin() + (size_t) n * k);
  }
}

// MASS::ginv(G): SVD (dgesdd as R's svd), tolerance sqrt(eps).
static std::vector<double> mass_ginv(const std::vector<double>& G, int n) {
  std::vector<double> A(G), s(n), u((size_t) n * n), vt((size_t) n * n);
  const char* jobz = "S";
  int info = 0, lwork = -1; double tmp = 0.0;
  std::vector<int> iwork(8 * (size_t) n);
  F77_CALL(dgesdd)(jobz, &n, &n, A.data(), &n, s.data(), u.data(), &n, vt.data(), &n,
                   &tmp, &lwork, iwork.data(), &info FCONE);
  lwork = (int) tmp; std::vector<double> work(lwork);
  F77_CALL(dgesdd)(jobz, &n, &n, A.data(), &n, s.data(), u.data(), &n, vt.data(), &n,
                   work.data(), &lwork, iwork.data(), &info FCONE);
  if (info != 0) LEWIS_ERROR("laplace_exact: dgesdd failed in ginv");
  double tol = std::sqrt(DBL_EPSILON);
  std::vector<double> out((size_t) n * n, 0.0);
  // v %*% (1/d[Positive] * t(u[, Positive]))
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      long double acc = 0.0L;
      for (int k = 0; k < n; ++k)
        if (s[k] > std::max(tol * s[0], 0.0)) acc += vt[k + (size_t) n * i] * ((1.0 / s[k]) * u[j + (size_t) n * k]);
      out[i + (size_t) n * j] = (double) acc;
    }
  return out;
}

struct ImportanceResult { double log_marg; double ess; };

static ImportanceResult log_marg_lik_importance_node_core(
    const double* X, int n, int p, int j, const std::vector<int>& pa, int M,
    double sigma_b2, double ig_alpha, double ig_beta, double nu) {
  double cn = std::log(1 + 1.0 / n);
  L1Fit fit = fit_node_l1_core(X, n, p, j, pa);
  int d_j = (int) pa.size();
  const double* xj = X + (size_t) n * j;
  double mu_t = std::log(fit.theta) - cn / 2;
  double sd_t = std::sqrt(cn);

  std::vector<double> Xpa, Sigma, Rsqrt, chol;
  bool chol_ok = true;
  if (d_j > 0) {
    Xpa.resize((size_t) n * d_j);
    for (int k = 0; k < d_j; ++k)
      std::copy(X + (size_t) n * pa[k], X + (size_t) n * (pa[k] + 1), Xpa.begin() + (size_t) n * k);
    std::vector<double> G; crossprod_sym(Xpa, n, d_j, G);
    std::vector<double> Gr(G);
    for (int k = 0; k < d_j; ++k) Gr[k + d_j * k] += 1e-8;
    std::vector<double> Ginv((size_t) d_j * d_j, 0.0);
    for (int k = 0; k < d_j; ++k) Ginv[k + d_j * k] = 1.0;
    if (!r_solve(Gr, d_j, Ginv, d_j)) Ginv = mass_ginv(G, d_j);
    double sc = ((nu - 2) / nu) * (fit.theta / 2);
    Sigma.resize(Ginv.size());
    for (size_t k = 0; k < Ginv.size(); ++k) Sigma[k] = sc * Ginv[k];
    // rmvnorm "eigen" square root:  R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
    std::vector<double> ev, V; r_eigen_sym(Sigma, d_j, ev, V);
    std::vector<double> W((size_t) d_j * d_j);      // t(V) * sqrt(ev): W[i, k] = V[k, i] * s[i]
    for (int i = 0; i < d_j; ++i) {
      double s = std::sqrt(std::max(ev[i], 0.0));
      for (int k = 0; k < d_j; ++k) W[i + d_j * k] = V[k + d_j * i] * s;
    }
    std::vector<double> P((size_t) d_j * d_j, 0.0);
    if (d_j == 1) {
      P[0] = V[0] * W[0];
    } else {
      const char* tN = "N"; const double one = 1.0, zero = 0.0;
      F77_CALL(dgemm)(tN, tN, &d_j, &d_j, &d_j, &one, V.data(), &d_j, W.data(), &d_j,
                      &zero, P.data(), &d_j FCONE FCONE);
    }
    Rsqrt.resize(P.size());
    for (int i = 0; i < d_j; ++i) for (int k = 0; k < d_j; ++k) Rsqrt[i + d_j * k] = P[k + d_j * i];
    // dmvt: chol(sigma) (dpotrf "U"), lower triangle zeroed
    chol = Sigma;
    const char* up = "U"; int info = 0;
    F77_CALL(dpotrf)(up, &d_j, chol.data(), &d_j, &info FCONE);
    if (info != 0) chol_ok = false;
    else for (int i = 0; i < d_j; ++i) for (int k = 0; k < i; ++k) chol[i + d_j * k] = 0.0;
  }
  double log_det_chol = 0.0;
  if (d_j > 0 && chol_ok) {
    std::vector<double> ld(d_j);
    for (int k = 0; k < d_j; ++k) ld[k] = std::log(chol[k + d_j * k]);
    log_det_chol = r_sum(ld.data(), d_j);
  }
  double sum_abs_xj = 0.0;
  if (d_j == 0) {
    std::vector<double> ax(n);
    for (int i = 0; i < n; ++i) ax[i] = std::fabs(xj[i]);
    sum_abs_xj = r_sum(ax.data(), n);
  }

  std::vector<double> log_w(M);
  std::vector<double> z(d_j), bj(d_j), diff(d_j), res(n), ares(n), fitted(n), b2(d_j);
  for (int m = 0; m < M; ++m) {
    double log_theta = Rf_rnorm(mu_t, sd_t);
    double theta_j = std::exp(log_theta);
    double log_prop = Rf_dnorm4(log_theta, mu_t, sd_t, 1) - log_theta;
    double log_prior = ig_alpha * std::log(ig_beta) - Rf_lgammafn(ig_alpha) -
                       (ig_alpha + 1) * std::log(theta_j) - ig_beta / theta_j;
    double log_lik;
    if (d_j == 0) {
      log_lik = -(double) n * std::log(2.0) - (double) n * std::log(theta_j) - sum_abs_xj / theta_j;
    } else {
      for (int k = 0; k < d_j; ++k) z[k] = norm_rand();                      // rnorm(d_j)
      // (1 x d) %*% R  ->  R's vector-matrix branch: dgemv("T", R, z)
      std::vector<double> zr(d_j, 0.0);
      const char* tT = "T"; const double one = 1.0, zero = 0.0; const int ione = 1;
      F77_CALL(dgemv)(tT, &d_j, &d_j, &one, Rsqrt.data(), &d_j, z.data(), &ione, &zero, zr.data(), &ione FCONE);
      double chi = Rf_rchisq(nu);                                           // rchisq(1, nu)
      double denom = std::sqrt(chi / nu);
      for (int k = 0; k < d_j; ++k) bj[k] = zr[k] / denom + fit.b[k];        // sims/denom, + delta(0), + b_hat
      // dmvt(b_j - b_hat, sigma = Sigma_hat, df = nu, log = TRUE)
      for (int k = 0; k < d_j; ++k) diff[k] = bj[k] - fit.b[k];
      double log_dmvt;
      if (!chol_ok) {
        bool all_zero = true;
        for (int k = 0; k < d_j; ++k) if (diff[k] != 0.0) { all_zero = false; break; }
        log_dmvt = all_zero ? INFINITY : -INFINITY;
      } else {
        std::vector<double> rx(diff);
        const char* L = "L"; const char* U = "U"; const char* T = "T"; const char* N = "N";
        F77_CALL(dtrsm)(L, U, T, N, &d_j, &ione, &one, chol.data(), &d_j, rx.data(), &d_j FCONE FCONE FCONE FCONE);
        long double rss = 0.0L;
        for (int k = 0; k < d_j; ++k) rss += rx[k] * rx[k];
        double rssd = (double) rss;
        log_dmvt = Rf_lgammafn((d_j + nu) / 2) -
                   (Rf_lgammafn(nu / 2) + log_det_chol + d_j / 2.0 * std::log(M_PI * nu)) -
                   0.5 * (nu + d_j) * std::log1p(rssd / nu);
      }
      log_prop = log_prop + log_dmvt;
      for (int k = 0; k < d_j; ++k) b2[k] = bj[k] * bj[k];
      log_prior = log_prior - 0.5 * d_j * std::log(2 * M_PI * sigma_b2) -
                  r_sum(b2.data(), d_j) / (2 * sigma_b2);
      F77_CALL(dgemv)("N", &n, &d_j, &one, Xpa.data(), &n, bj.data(), &ione, &zero, fitted.data(), &ione FCONE);
      for (int i = 0; i < n; ++i) { res[i] = xj[i] - fitted[i]; ares[i] = std::fabs(res[i]); }
      log_lik = -(double) n * std::log(2.0) - (double) n * std::log(theta_j) - r_sum(ares.data(), n) / theta_j;
    }
    log_w[m] = log_lik + log_prior - log_prop;
  }
  double M_log = log_w[0];
  for (int m = 1; m < M; ++m) if (log_w[m] > M_log) M_log = log_w[m];
  std::vector<double> w(M), w2(M);
  for (int m = 0; m < M; ++m) { w[m] = std::exp(log_w[m] - M_log); w2[m] = w[m] * w[m]; }
  ImportanceResult out;
  out.log_marg = M_log + std::log(r_mean(w.data(), M));
  double sw = r_sum(w.data(), M);
  out.ess = sw * sw / r_sum(w2.data(), M);
  return out;
}

// =============================================================================
// 4.  THE SAMPLER WITH THE IMPORTANCE-SAMPLING SCORE
// =============================================================================
// The proposal adds an edge on a uniformly chosen non-adjacent ordered pair,
// deletes a uniformly chosen edge, or reverses a uniformly chosen edge, each
// with probability 1/3. An impossible or cyclic proposal is rejected. A
// valid move is accepted with probability min(1, exp(score(new) -
// score(old)) r), where r is the Hastings ratio. The score of a node is its
// importance-sampling estimate plus its share of the log prior. A move draws
// new estimates for the nodes whose parent sets change and keeps the
// estimates of the other nodes, as the pseudo-marginal chain requires.

struct ScoreSpec {
  int prior;       // 0 = uniform, 1 = complexity
  double alpha, d; // constants of the complexity prior
  int M;           // number of importance draws per node
};

// Score of node k with parent set pa. The constants 100, 1, 1 and 5 are the
// prior variance of the coefficients, the Inv-Gamma(1, 1) prior of the scale
// and the degrees of freedom of the multivariate t importance distribution.
static double score_node(const double* X, int n, int p, int k,
                         const std::vector<int>& pa, const ScoreSpec& sp) {
  ImportanceResult r = log_marg_lik_importance_node_core(X, n, p, k, pa, sp.M,
                                                         100.0, 1.0, 1.0, 5.0);
  double node_prior = (sp.prior == 1) ?
    -R_pow((double) n, sp.alpha) * sp.d * (double) pa.size() : 0.0;
  return r.log_marg + node_prior;
}

struct MCMCResult {
  std::vector<int> best_dag;
  double best_score;
  double acceptance_rate, invalid_rate;
};

static MCMCResult exact_mcmc_core(const double* X, int n, int p,
                                  std::vector<int> current, int n_iter,
                                  const ScoreSpec& sp, bool verbose) {
  if (!is_acyclic_core(current, p)) LEWIS_ERROR("The start graph has a directed cycle.");
  if (n_iter < 1) LEWIS_ERROR("The number of iterations must be positive.");
  int report_every = std::max(1, n_iter / 10);

  // Scores of the start graph, node by node in index order.
  int n_pairs = p * (p - 1);
  std::vector<double> per(p);
  for (int j = 0; j < p; ++j) per[j] = score_node(X, n, p, j, parents_of(current, p, j), sp);
  double total = r_sum(per.data(), p);
  int n_edges = 0;
  for (size_t k = 0; k < current.size(); ++k) n_edges += current[k];

  MCMCResult out_;
  out_.best_dag = current; out_.best_score = total;
  int n_accept = 0, n_invalid = 0;

  std::vector<int> new_dag(current.size());
  std::vector<double> new_per(p);

  for (int iter = 1; iter <= n_iter; ++iter) {
    // The move type: sample.int(3, 1).
    int move = (int) (R_unif_index(3.0) + 1);
    int n_absent = n_pairs - 2 * n_edges;
    new_dag = current;
    int changed[2]; int n_changed = 0;
    bool valid = false; double log_q = 0.0;

    if (move == 1) {
      // Addition on a non-adjacent ordered pair.
      int i1 = -1, i2 = -1;
      if (n_absent > 0) {
        bool found = false;
        for (int tries = 0; tries < 64 && !found; ++tries) {
          // sample.int(p, 2): the loop of R's do_sample without replacement.
          int nn = p; int x_last = nn - 1;
          int j1 = (int) R_unif_index(nn);
          int a1 = j1;
          int replaced = x_last; nn--;
          int j2 = (int) R_unif_index(nn);
          int a2 = (j2 == j1) ? replaced : j2;
          if (current[a1 + p * a2] == 0 && current[a2 + p * a1] == 0) { i1 = a1; i2 = a2; found = true; }
        }
        if (!found) {
          // All non-adjacent ordered pairs in column-major order, and one of
          // them chosen uniformly.
          std::vector<int> ci, cj;
          for (int j = 0; j < p; ++j) for (int i = 0; i < p; ++i)
            if (i != j && current[i + p * j] == 0 && current[j + p * i] == 0) { ci.push_back(i); cj.push_back(j); }
          int pick = (int) (R_unif_index((double) ci.size()) + 1) - 1;
          i1 = ci[pick]; i2 = cj[pick];
        }
        new_dag[i1 + p * i2] = 1;
        if (is_acyclic_core(new_dag, p)) {
          changed[0] = i2; n_changed = 1; valid = true;
          log_q = std::log((double) n_absent) - std::log((double) (n_edges + 1));
        }
      }
    } else {
      // Deletion or reversal of an edge chosen uniformly, as
      // which(A == 1, arr.ind = TRUE)[sample.int(|gamma|, 1), ].
      if (n_edges > 0) {
        int pick = (int) (R_unif_index((double) n_edges) + 1);
        int cnt = 0, i1 = -1, i2 = -1;
        for (int j = 0; j < p && i1 < 0; ++j) for (int i = 0; i < p; ++i)
          if (current[i + p * j] == 1 && ++cnt == pick) { i1 = i; i2 = j; break; }
        if (move == 2) {
          new_dag[i1 + p * i2] = 0;
          changed[0] = i2; n_changed = 1; valid = true;
          log_q = std::log((double) n_edges) - std::log((double) (n_absent + 2));
        } else {
          new_dag[i1 + p * i2] = 0; new_dag[i2 + p * i1] = 1;
          if (is_acyclic_core(new_dag, p)) {
            changed[0] = i1; changed[1] = i2; n_changed = 2; valid = true;
            log_q = 0.0;
          }
        }
      }
    }

    if (valid) {
      // New estimates for the changed nodes only.
      new_per = per;
      for (int c = 0; c < n_changed; ++c) {
        int k = changed[c];
        new_per[k] = score_node(X, n, p, k, parents_of(new_dag, p, k), sp);
      }
      double new_total = r_sum(new_per.data(), p);
      double log_acc = (new_total - total) + log_q;
      // Acceptance with log(runif(1)).
      if (!std::isnan(log_acc) && std::log(Rf_runif(0.0, 1.0)) < log_acc) {
        current = new_dag; per = new_per; total = new_total;
        n_edges = (move == 1) ? n_edges + 1 : (move == 2) ? n_edges - 1 : n_edges;
        n_accept++;
        if (total > out_.best_score) { out_.best_dag = current; out_.best_score = total; }
      }
    } else {
      n_invalid++;
    }

    if (verbose && iter % report_every == 0) {
      LEWIS_PRINTF("  iteration %d of %d, acceptance %.3f, edges %d, best score %.2f\n",
                   iter, n_iter, (double) n_accept / iter, n_edges, out_.best_score);
    }
  }
  out_.acceptance_rate = (double) n_accept / n_iter;
  out_.invalid_rate = (double) n_invalid / n_iter;
  return out_;
}

// =============================================================================
// 5.  FUNCTIONS CALLED FROM R
// =============================================================================
using namespace Rcpp;

// p x p integer matrix to a column-major vector, and back.
static std::vector<int> mat_to_vec(const IntegerMatrix& A) {
  std::vector<int> v(A.begin(), A.end());
  return v;
}
static IntegerMatrix vec_to_mat(const std::vector<int>& v, int p) {
  IntegerMatrix A(p, p);
  std::copy(v.begin(), v.end(), A.begin());
  return A;
}

// Hand over the Fortran routines rqbr and rqfnb of quantreg, as external
// pointers from getNativeSymbolInfo(..., PACKAGE = "quantreg")$address.
// [[Rcpp::export]]
void exact_set_quantreg_symbols(SEXP rqbr_ptr, SEXP rqfnb_ptr) {
  if (TYPEOF(rqbr_ptr) != EXTPTRSXP || TYPEOF(rqfnb_ptr) != EXTPTRSXP)
    stop("Expected two external pointers.");
  RQBR_PTR  = (rqbr_fn)  R_ExternalPtrAddr(rqbr_ptr);
  RQFNB_PTR = (rqfnb_fn) R_ExternalPtrAddr(rqfnb_ptr);
  if (RQBR_PTR == NULL || RQFNB_PTR == NULL) stop("A quantreg routine was not found.");
}

// The exact chain. prior is "uniform" or "complexity". Returns the
// highest-scoring DAG visited over all iterations and its score.
// [[Rcpp::export]]
List exact_mcmc_cpp(NumericMatrix X, IntegerMatrix start, int n_iter,
                    std::string prior, double alpha, double d, int M,
                    bool verbose) {
  int n = X.nrow(), p = X.ncol();
  if (start.nrow() != p || start.ncol() != p) stop("start must be a p x p matrix.");
  if (prior != "uniform" && prior != "complexity") stop("prior must be uniform or complexity.");
  if (M < 1) stop("M must be positive.");
  ScoreSpec sp;
  sp.prior = (prior == "complexity") ? 1 : 0;
  sp.alpha = alpha; sp.d = d; sp.M = M;
  MCMCResult res = exact_mcmc_core(X.begin(), n, p, mat_to_vec(start), n_iter, sp, verbose);
  return List::create(_["best_dag"] = vec_to_mat(res.best_dag, p),
                      _["best_score"] = res.best_score,
                      _["acceptance_rate"] = res.acceptance_rate,
                      _["invalid_rate"] = res.invalid_rate);
}
