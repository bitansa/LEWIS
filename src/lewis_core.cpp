// =============================================================================
// lewis_core.cpp
//
// The compiled engine of the LEWIS package.
//
// Contents
//   1. Summaries that follow R's own arithmetic (sum, mean, var, quantile).
//   2. DAG utilities.
//   3. Least absolute deviation (L1) regression of one node on its parents.
//   4. The Laplace-approximated score s_n of a DAG, node by node.
//   5. The Metropolis-Hastings sampler over DAGs.
//   6. The Shapiro-Wilk test.
//   7. Step 1 of the post-processing: the non-Gaussian nodes of a DAG.
//   8. Step 2: CPDAG, Meek's rules and the RPDAG. Step 3: the SHD.
//   9. Step 3: the point estimate (SHD medoid) and its non-Gaussian set.
//  10. Functions called from R.
//
// Numerical conventions
//   Random numbers come from R's generator through the R API (unif_rand,
//   R_unif_index, Rf_runif). The chain therefore depends on set.seed() in
//   R. Sums use long double accumulators where R does. Matrix products,
//   determinants and linear solves call the same BLAS and LAPACK routines
//   that R calls for the same expressions. This keeps the results equal to
//   those of a pure R implementation of the same steps.
//
// Third-party code
//   The L1 regressions call the Fortran routines rqbr and rqfnb of the R
//   package quantreg (version 6.1, Roger Koenker, GPL (>= 2)). Copies of
//   these routines are in src/rqbr.f and src/rqfnb.f.
//   Section 6 contains the Shapiro-Wilk algorithm of R's file swilk.c
//   (algorithm AS R94), copyright (C) the R Core Team, GPL (>= 2).
//
// Adjacency matrices are p x p and stored by column. A[i + p * j] = 1 means
// the edge i -> j (the row is the parent).
// =============================================================================

#ifndef USE_FC_LEN_T
#  define USE_FC_LEN_T
#endif
#include <Rcpp.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>
#include <Rmath.h>
#ifndef FCONE
#  define FCONE
#endif

#include <vector>
#include <string>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstring>
#include <algorithm>
#include <unordered_map>

// The Fortran L1 solvers of quantreg (src/rqbr.f and src/rqfnb.f).
// The routine rqbr has a LOGICAL argument (lci1). It is called through the
// wrapper lwrqbr of src/rqbr_wrap.f, which takes this flag as an INTEGER.
// So every argument below is a Fortran INTEGER or DOUBLE PRECISION.
extern "C" {
void F77_NAME(lwrqbr)(int* m, int* nn, int* m5, int* n3, int* n4,
                      double* a, double* b, double* t, double* toler,
                      int* ift, double* x, double* e, int* s, double* wa,
                      double* wb, int* nsol, int* ndsol, double* sol,
                      double* dsol, int* lsol, int* h, double* qn,
                      double* cutoff, double* ci, double* tnmat,
                      double* big, int* ilci1);
void F77_NAME(rqfnb)(int* n, int* p, double* a, double* c, double* rhs,
                     double* d, double* u, double* beta, double* eps,
                     double* wn, double* wp, int* nit, int* info);
}

// =============================================================================
// 1.  SUMMARIES THAT FOLLOW R'S ARITHMETIC
// =============================================================================
// R accumulates sums in long double. mean() adds a second refinement pass.
// var() centers with that mean and also accumulates in long double.

// sum(x) of R.
static double r_sum(const double* x, int n) {
  long double s = 0.0L;
  for (int i = 0; i < n; ++i) s += x[i];
  return (double) s;
}

// mean(x) of R, with its refinement pass.
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

// var(x) of R for one column: two-pass mean, then the sum of squared
// deviations divided by n - 1.
static double r_var(const double* x, int n) {
  long double sum = 0.0L, tmp;
  for (int k = 0; k < n; ++k) sum += x[k];
  tmp = sum / n;
  if (std::isfinite((double) tmp)) {
    sum = 0.0L;
    for (int k = 0; k < n; ++k) sum += (x[k] - tmp);
    tmp = tmp + sum / n;
  }
  long double xxm = (double) tmp;
  sum = 0.0L;
  for (int k = 0; k < n; ++k) sum += (x[k] - xxm) * (x[k] - xxm);
  return (double) (sum / (n - 1));
}

// quantile(x, prob, type = 7) of R on a sorted vector.
static double r_quantile7_sorted(const std::vector<double>& xs, double prob) {
  int n = (int) xs.size();
  double index = 1.0 + std::max(n - 1, 0) * prob;
  int lo = (int) std::floor(index);
  int hi = (int) std::ceil(index);
  double qs = xs[lo - 1];
  if (index > lo && xs[hi - 1] != qs) {
    double h = index - lo;
    qs = (1 - h) * qs + h * xs[hi - 1];
  }
  return qs;
}

// =============================================================================
// 2.  DAG UTILITIES
// =============================================================================

// TRUE if the directed graph A has no directed cycle (Kahn's algorithm).
static bool is_acyclic_core(const std::vector<int>& A, int p) {
  std::vector<int> indeg(p, 0);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i) if (A[i + p * j]) indeg[j]++;
  std::vector<int> stack;
  stack.reserve(p);
  for (int j = 0; j < p; ++j) if (indeg[j] == 0) stack.push_back(j);
  int done = 0;
  while (!stack.empty()) {
    int v = stack.back();
    stack.pop_back();
    ++done;
    for (int j = 0; j < p; ++j)
      if (A[v + p * j] && --indeg[j] == 0) stack.push_back(j);
  }
  return done == p;
}

// The parents of node j in increasing order (0-based).
static std::vector<int> parents_of(const std::vector<int>& A, int p, int j) {
  std::vector<int> pa;
  for (int i = 0; i < p; ++i) if (A[i + p * j] == 1) pa.push_back(i);
  return pa;
}

// The positions (0-based, by column) of the nonzero entries of A.
static std::vector<int> positions_of(const std::vector<int>& A) {
  std::vector<int> pos;
  for (size_t k = 0; k < A.size(); ++k) if (A[k] != 0) pos.push_back((int) k);
  return pos;
}

// =============================================================================
// 3.  LEAST ABSOLUTE DEVIATION (L1) REGRESSION
// =============================================================================
// The fit equals quantreg::rq.fit(x, y, tau = 0.5, method = "br"). If that
// function would stop with "Singular design matrix", the Frisch-Newton
// method ("fn") is used. If the coefficients are not finite, a ridge least
// squares fit is used.

// qr(x)$rank of R (LINPACK dqrdc2 with tol = 1e-7). Only the rank is needed,
// so pivoted columns are skipped.
static int linpack_qr_rank(std::vector<double> x, int n, int p, double tol) {
  std::vector<double> qraux(p), work2(p);
  for (int j = 0; j < p; ++j) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += x[i + n * j] * x[i + n * j];
    qraux[j] = std::sqrt(s);
    work2[j] = qraux[j] == 0.0 ? 1.0 : qraux[j];
  }
  int lup = std::min(n, p), k = p + 1;
  for (int l = 1; l <= lup; ++l) {
    // Move columns with a negligible remaining norm to the end.
    while (!(l >= k || qraux[l - 1] >= work2[l - 1] * tol)) {
      std::vector<double> col(x.begin() + n * (l - 1), x.begin() + n * l);
      double q = qraux[l - 1], w = work2[l - 1];
      for (int j = l; j < p; ++j) {
        std::copy(x.begin() + n * j, x.begin() + n * (j + 1), x.begin() + n * (j - 1));
        qraux[j - 1] = qraux[j];
        work2[j - 1] = work2[j];
      }
      std::copy(col.begin(), col.end(), x.begin() + n * (p - 1));
      qraux[p - 1] = q;
      work2[p - 1] = w;
      --k;
    }
    if (l == n) break;
    // Householder transformation for column l.
    double nrmxl = 0.0;
    for (int i = l - 1; i < n; ++i) nrmxl += x[i + n * (l - 1)] * x[i + n * (l - 1)];
    nrmxl = std::sqrt(nrmxl);
    if (nrmxl == 0.0) continue;
    double* xl = &x[(l - 1) + n * (l - 1)];
    if (xl[0] != 0.0) nrmxl = (xl[0] < 0) ? -std::fabs(nrmxl) : std::fabs(nrmxl);
    for (int i = 0; i < n - l + 1; ++i) xl[i] *= 1.0 / nrmxl;
    xl[0] = 1.0 + xl[0];
    // Apply it to the remaining columns and update their norms.
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

// The result of one L1 fit.
struct L1Fit {
  std::vector<double> b;          // coefficients, one per parent
  double theta;                   // mean absolute residual, floored at 1e-12
  std::vector<double> resid;      // the finite residuals
};

// The residual y - Xpa %*% b as R computes it (dgemv, then elementwise).
static void residual_dgemv(const std::vector<double>& Xpa, int n, int d,
                           const std::vector<double>& y,
                           const std::vector<double>& b,
                           std::vector<double>& res) {
  std::vector<double> fitted(n, 0.0);
  const char* tN = "N";
  const double one = 1.0, zero = 0.0;
  const int ione = 1;
  F77_CALL(dgemv)(tN, &n, &d, &one, Xpa.data(), &n, b.data(), &ione, &zero,
                  fitted.data(), &ione FCONE);
  res.resize(n);
  for (int i = 0; i < n; ++i) res[i] = y[i] - fitted[i];
}

// crossprod(Xpa) of R: dsyrk on the upper triangle, then mirrored.
static void crossprod_sym(const std::vector<double>& Xpa, int n, int d,
                          std::vector<double>& G) {
  G.assign((size_t) d * d, 0.0);
  const char* up = "U";
  const char* tr = "T";
  const double one = 1.0, zero = 0.0;
  F77_CALL(dsyrk)(up, tr, &d, &n, &one, Xpa.data(), &n, &zero, G.data(), &d FCONE FCONE);
  for (int i = 1; i < d; ++i)
    for (int j = 0; j < i; ++j) G[i + d * j] = G[j + d * i];
}

// solve(A, B) of R: dgesv, then the reciprocal condition number check with
// dgecon. Returns false where R would stop with an error.
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
  // The pivot vector is not needed any more. It serves as integer workspace.
  F77_CALL(dgecon)(one, &n, A.data(), &n, &anorm, &rcond, work.data(), ipiv.data(), &info FCONE);
  if (rcond < DBL_EPSILON) return false;
  return true;
}

// Ridge least squares: solve(crossprod(Xpa) + 1e-8 I, crossprod(Xpa, y)).
static std::vector<double> ridge_ols(const std::vector<double>& Xpa, int n, int d,
                                     const std::vector<double>& y) {
  std::vector<double> G;
  crossprod_sym(Xpa, n, d, G);
  for (int k = 0; k < d; ++k) G[k + d * k] += 1e-8;
  std::vector<double> Xty(d, 0.0);
  const char* tT = "T";
  const double one = 1.0, zero = 0.0;
  const int ione = 1;
  F77_CALL(dgemv)(tT, &n, &d, &one, Xpa.data(), &n, y.data(), &ione, &zero, Xty.data(), &ione FCONE);
  if (!r_solve(G, d, Xty, 1))
    Rcpp::stop("The least squares fit of a node on its parents failed. "
               "The data may contain linearly dependent columns.");
  return Xty;
}

// Coefficients of quantreg::rq.fit.br(x, y, tau = 0.5). Returns false where
// that function would stop with "Singular design matrix".
static bool rq_fit_br(const std::vector<double>& Xpa, int n, int p,
                      const std::vector<double>& y, std::vector<double>& coef) {
  if (linpack_qr_rank(Xpa, n, p, 1e-7) < p) return false;
  // The arguments that rq.fit.br() passes to rqbr.
  double tol = R_pow(DBL_EPSILON, 2.0 / 3.0);
  double big = DBL_MAX;
  int nsol = 2, ndsol = 2;
  int m = n, nn = p, m5 = n + 5, n3 = p + 3, n4 = p + 4;
  double tau = 0.5;
  int ift = 1;
  std::vector<double> a(Xpa), b(y);
  std::vector<double> x(p, 0.0), e(n, 0.0);
  std::vector<int> s(n, 0);
  std::vector<double> wa((size_t) (n + 5) * (p + 4), 0.0), wb(n, 0.0);
  std::vector<double> sol((size_t) (p + 3) * nsol, 0.0), dsol((size_t) n * ndsol, 0.0);
  int lsol = 0;
  std::vector<int> h((size_t) p * nsol, 0);
  std::vector<double> qn(p, 0.0);
  double cutoff = 0.0;
  std::vector<double> ci(4 * (size_t) p, 0.0), tnmat(4 * (size_t) p, 0.0);
  // No confidence intervals are computed (lci1 = FALSE in rq.fit.br()).
  int lci1 = 0;
  // The Barrodale-Roberts simplex solver, through the wrapper of rqbr.
  F77_CALL(lwrqbr)(&m, &nn, &m5, &n3, &n4, a.data(), b.data(), &tau, &tol, &ift,
                   x.data(), e.data(), s.data(), wa.data(), wb.data(), &nsol, &ndsol,
                   sol.data(), dsol.data(), &lsol, h.data(), qn.data(), &cutoff,
                   ci.data(), tnmat.data(), &big, &lci1);
  coef = x;
  return true;
}

// Coefficients of quantreg::rq.fit.fnb(x, y, tau = 0.5).
static void rq_fit_fnb(const std::vector<double>& Xpa, int n, int p,
                       const std::vector<double>& y, std::vector<double>& coef) {
  // The arguments that rq.fit.fnb() passes to rqfnb.
  double tau = 0.5, beta = 0.99995, eps = 1e-06;
  std::vector<double> a((size_t) n * p);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < p; ++j) a[j + (size_t) p * i] = Xpa[i + (size_t) n * j];
  std::vector<double> c(n);
  for (int i = 0; i < n; ++i) c[i] = -y[i];
  std::vector<double> rhs(p);
  for (int j = 0; j < p; ++j) rhs[j] = (1 - tau) * r_sum(&Xpa[(size_t) n * j], n);
  std::vector<double> d(n, 1.0), u(n, 1.0), wn(10 * (size_t) n, 0.0);
  for (int i = 0; i < n; ++i) wn[i] = 1 - tau;
  std::vector<double> wp((size_t) (p + 3) * p, 0.0);
  int nit[3] = {0, 0, 0}, info = 0, nn = n, pp = p;
  // The Frisch-Newton interior point solver.
  F77_CALL(rqfnb)(&nn, &pp, a.data(), c.data(), rhs.data(), d.data(), u.data(),
                  &beta, &eps, wn.data(), wp.data(), nit, &info);
  coef.resize(p);
  for (int j = 0; j < p; ++j) coef[j] = -wp[j];
}

// The L1 fit of node j on its parents pa. With no parents, the residuals
// are the data of node j, because the model has no intercept.
static L1Fit fit_node_l1_core(const double* X, int n, int p, int j,
                              const std::vector<int>& pa) {
  L1Fit f;
  int d = (int) pa.size();
  // No parents: only the scale is estimated.
  if (d == 0) {
    f.resid.assign(X + (size_t) n * j, X + (size_t) n * (j + 1));
    std::vector<double> ar(n);
    for (int i = 0; i < n; ++i) ar[i] = std::fabs(f.resid[i]);
    f.theta = std::max(r_mean(ar.data(), n), 1e-12);
    return f;
  }
  // Response and design matrix.
  std::vector<double> y(X + (size_t) n * j, X + (size_t) n * (j + 1));
  std::vector<double> Xpa((size_t) n * d);
  for (int k = 0; k < d; ++k)
    std::copy(X + (size_t) n * pa[k], X + (size_t) n * (pa[k] + 1), Xpa.begin() + (size_t) n * k);

  // Barrodale-Roberts fit, with the two fallbacks.
  std::vector<double> coef;
  bool ok = rq_fit_br(Xpa, n, d, y, coef);
  if (!ok) rq_fit_fnb(Xpa, n, d, y, coef);
  bool finite = true;
  for (int k = 0; k < d; ++k) if (!std::isfinite(coef[k])) { finite = false; break; }
  if (!finite) coef = ridge_ols(Xpa, n, d, y);
  for (int k = 0; k < d; ++k) if (!std::isfinite(coef[k])) coef[k] = 0.0;

  // Residuals and the scale estimate.
  std::vector<double> res;
  residual_dgemv(Xpa, n, d, y, coef, res);
  f.resid.clear();
  f.resid.reserve(n);
  for (int i = 0; i < n; ++i) if (std::isfinite(res[i])) f.resid.push_back(res[i]);
  int nr = (int) f.resid.size();
  if (nr == 0) {
    f.theta = 1e-12;
  } else {
    std::vector<double> ar(nr);
    for (int i = 0; i < nr; ++i) ar[i] = std::fabs(f.resid[i]);
    f.theta = std::max(r_mean(ar.data(), nr), 1e-12);
  }
  f.b = coef;
  return f;
}

// =============================================================================
// 4.  THE LAPLACE-APPROXIMATED SCORE
// =============================================================================
// The score of node j with parent set S is
//   s_{n,j}(S) = -n log 2 - n log theta_hat - n - (|S| + 1)/2 log n
//                + c_hat_j(S) - n^alpha d |S|,
// and s_n of a DAG is the sum over the nodes. The term c_hat_j(S) is the
// estimate of the constant of the Laplace approximation theorem of the
// paper. It is the log prior density of the reparametrized parameters
// (eta, w) = (1 / theta, b / theta) at the estimate, plus
// (|S| + 1)/2 log(2 pi), minus half the log determinant of the information
// matrix J_j.

// The priors are fixed: b_j ~ N(0, 100 I) and theta_j ~ Inv-Gamma(1, 1).
static const double PRIOR_B_VARIANCE = 100.0;
static const double PRIOR_THETA_SHAPE = 1.0;
static const double PRIOR_THETA_RATE = 1.0;

// Gaussian kernel estimate of the residual density at zero, with the
// bandwidth 0.9 min(sd, IQR / 1.349) n^(-1/5). The neutral value 1 is used
// when fewer than 5 residuals are available or the spread is zero.
static double resid_density_at_zero_core(const std::vector<double>& resid) {
  std::vector<double> r;
  r.reserve(resid.size());
  for (size_t i = 0; i < resid.size(); ++i) if (std::isfinite(resid[i])) r.push_back(resid[i]);
  int n = (int) r.size();
  if (n < 5) return 1.0;
  // Robust scale and bandwidth.
  double s = std::sqrt(r_var(r.data(), n));
  std::vector<double> rs(r);
  std::sort(rs.begin(), rs.end());
  double iqr = (r_quantile7_sorted(rs, 0.75) - r_quantile7_sorted(rs, 0.25)) / 1.349;
  double scale_r = std::min(s, (std::isfinite(iqr) && iqr > 0) ? iqr : s);
  if (!std::isfinite(scale_r) || scale_r <= 0) return 1.0;
  double h = 0.9 * scale_r * R_pow((double) n, -1.0 / 5.0);
  if (!std::isfinite(h) || h <= 0) return 1.0;
  // Kernel average, floored at 1e-8.
  std::vector<double> dn(n);
  for (int i = 0; i < n; ++i) dn[i] = Rf_dnorm4(r[i] / h, 0.0, 1.0, 0);
  return std::max(r_mean(dn.data(), n) / h, 1e-8);
}

// The information matrix J_j of node j, of size (|S| + 1) x (|S| + 1).
static std::vector<double> build_J_core(const std::vector<double>& b, double theta,
                                        const double* X, int n, int p,
                                        const std::vector<int>& pa, double f0) {
  int d = (int) pa.size();
  double eta = 1.0 / theta;
  // No parents: J_j is the scalar theta_hat^2.
  if (d == 0) return std::vector<double>(1, theta * theta);
  // Gram matrix M = X_pa' X_pa / n, v = M w and S = w' M w.
  std::vector<double> Xpa((size_t) n * d);
  for (int k = 0; k < d; ++k)
    std::copy(X + (size_t) n * pa[k], X + (size_t) n * (pa[k] + 1), Xpa.begin() + (size_t) n * k);
  std::vector<double> w(d);
  for (int k = 0; k < d; ++k) w[k] = b[k] / theta;
  std::vector<double> M;
  crossprod_sym(Xpa, n, d, M);
  for (size_t k = 0; k < M.size(); ++k) M[k] = M[k] / n;
  std::vector<double> v(d, 0.0);
  const char* tN = "N";
  const double one = 1.0, zero = 0.0;
  const int ione = 1;
  F77_CALL(dgemv)(tN, &d, &d, &one, M.data(), &d, w.data(), &ione, &zero, v.data(), &ione FCONE);
  std::vector<double> wv(d);
  for (int k = 0; k < d; ++k) wv[k] = w[k] * v[k];
  double S = r_sum(wv.data(), d);
  // Fill the blocks of J_j.
  int D = d + 1;
  std::vector<double> J((size_t) D * D, 0.0);
  J[0] = 2 * f0 / R_pow(eta, 3.0) * S + 1 / (eta * eta);
  for (int k = 0; k < d; ++k) {
    double off = -2 * f0 / (eta * eta) * v[k];
    J[0 + D * (k + 1)] = off;
    J[(k + 1) + D * 0] = off;
  }
  for (int l = 0; l < d; ++l)
    for (int k = 0; k < d; ++k)
      J[(k + 1) + D * (l + 1)] = 2 * f0 / eta * M[k + d * l];
  return J;
}

// determinant(A, logarithm = TRUE)$modulus of R (LU factorization by dgetrf).
// A singular matrix gives -Inf. The score is then not finite, and the
// caller stops with an error.
static double log_det_modulus(std::vector<double> A, int n) {
  std::vector<int> jpvt(n);
  int info = 0;
  F77_CALL(dgetrf)(&n, &n, A.data(), &n, jpvt.data(), &info);
  if (info < 0) Rcpp::stop("Internal error: invalid argument in the LU factorization.");
  if (info > 0) return -INFINITY;
  double modulus = 0.0;
  for (int i = 0; i < n; ++i) {
    double dii = A[i * (n + 1)];
    modulus += std::log(dii < 0 ? -dii : dii);
  }
  return modulus;
}

// The estimate c_hat_j of node j.
static double node_c_gamma_core(const std::vector<double>& b, double theta,
                                const double* X, int n, int p,
                                const std::vector<int>& pa,
                                const std::vector<double>& resid) {
  const double sigma_b2 = PRIOR_B_VARIANCE;
  const double ig_alpha = PRIOR_THETA_SHAPE;
  const double ig_beta = PRIOR_THETA_RATE;
  int d_j = (int) pa.size();
  int d_total = d_j + 1;
  // Normalizing constant of the Gaussian integral.
  double gauss_norm = d_total / 2.0 * std::log(2 * M_PI);
  // Log prior density of b and theta at the estimate.
  double log_prior_b = 0.0;
  if (d_j > 0) {
    std::vector<double> b2(d_j);
    for (int k = 0; k < d_j; ++k) b2[k] = b[k] * b[k];
    log_prior_b = -(double) d_j / 2.0 * std::log(2 * M_PI * sigma_b2) -
                  r_sum(b2.data(), d_j) / (2 * sigma_b2);
  }
  double log_prior_theta = ig_alpha * std::log(ig_beta) - Rf_lgammafn(ig_alpha) -
                           (ig_alpha + 1) * std::log(theta) - ig_beta / theta;
  // Jacobian of the map (eta, w) -> (theta, b).
  double log_jacobian = (double) (d_j + 2) * std::log(theta);
  double log_prior_etaw = log_prior_b + log_prior_theta + log_jacobian;
  // Information matrix and its log determinant.
  double f0 = (d_j == 0) ? 1.0 : resid_density_at_zero_core(resid);
  std::vector<double> J = build_J_core(b, theta, X, n, p, pa, f0);
  double log_det_J = log_det_modulus(J, d_j + 1);
  return gauss_norm + log_prior_etaw - 0.5 * log_det_J;
}

// The score s_{n,j}(S) of node j. With d = 0 (uniform prior) the prior term
// is skipped.
static double node_score_core(const double* X, int n, int p, int j,
                              const std::vector<int>& pa,
                              double alpha, double d) {
  // L1 fit of node j on its parents.
  L1Fit fit = fit_node_l1_core(X, n, p, j, pa);
  int d_j = (int) pa.size();
  double nn = (double) n;
  // Maximized log-likelihood and the log n penalty.
  double base = -nn * std::log(2.0) - nn * std::log(fit.theta) - nn -
                (1.0 + d_j) / 2.0 * std::log(nn);
  // The share of the complexity prior.
  if (d != 0.0) base = base - R_pow(nn, alpha) * d * d_j;
  // The constant c_hat_j.
  double c = node_c_gamma_core(fit.b, fit.theta, X, n, p, pa, fit.resid);
  return base + c;
}

// =============================================================================
// 5.  THE METROPOLIS-HASTINGS SAMPLER
// =============================================================================
// Each iteration proposes one of three moves with probability 1/3 each:
// add an edge on a uniformly chosen non-adjacent ordered pair, delete a
// uniformly chosen edge, or reverse a uniformly chosen edge. A move that is
// impossible or creates a cycle is rejected. Otherwise the move is accepted
// with probability min(1, exp(s_n(new) - s_n(old)) r), where r is the
// Hastings ratio. Only the nodes whose parent sets change are scored again.
//
// The retained draws are stored compactly. Each distinct DAG is stored once,
// as the positions of its edges. Each draw stores the index of its DAG.

struct ChainResult {
  std::vector< std::vector<int> > dags;   // distinct DAGs, edge positions
  std::vector<int> draw_dag;              // DAG index of each draw (1-based)
  std::vector<int> draw_iter;             // iteration of each draw
  std::vector<double> draw_score;         // s_n of each draw
  std::vector<int> draw_edges;            // number of edges of each draw
  std::vector<int> best_dag;
  double best_score;
  int n_accept, n_invalid;
  int moves[3][2];
};

// Turn a list of edge positions into a key for the hash table.
static std::string key_of_positions(const std::vector<int>& pos) {
  std::string key(pos.size() * sizeof(int), '\0');
  if (!pos.empty()) std::memcpy(&key[0], pos.data(), pos.size() * sizeof(int));
  return key;
}

static ChainResult mcmc_core(const double* X, int n, int p,
                             std::vector<int> current,
                             int n_iter, int burn_in, int thin,
                             double alpha, double d, bool verbose) {
  ChainResult out;
  int report_every = std::max(1, n_iter / 10);
  int n_pairs = p * (p - 1);

  // Scores of the start graph.
  std::vector<double> per(p);
  for (int j = 0; j < p; ++j) {
    per[j] = node_score_core(X, n, p, j, parents_of(current, p, j), alpha, d);
    if (!std::isfinite(per[j]))
      Rcpp::stop("The score of the start graph is not finite. "
                 "Check the data for linearly dependent columns.");
  }
  double total = r_sum(per.data(), p);
  int n_edges = 0;
  for (size_t k = 0; k < current.size(); ++k) n_edges += current[k];

  // Bookkeeping.
  out.best_dag = current;
  out.best_score = total;
  int n_accept = 0, n_invalid = 0;
  for (int a = 0; a < 3; ++a) for (int b = 0; b < 2; ++b) out.moves[a][b] = 0;
  int n_ret = (n_iter > burn_in) ? (n_iter - burn_in) / thin : 0;
  out.draw_dag.reserve(n_ret);
  out.draw_iter.reserve(n_ret);
  out.draw_score.reserve(n_ret);
  out.draw_edges.reserve(n_ret);
  std::unordered_map<std::string, int> seen;

  std::vector<int> new_dag(current.size());
  std::vector<double> new_per(p);

  for (int iter = 1; iter <= n_iter; ++iter) {
    // Allow the user to interrupt long runs.
    if (iter % 1000 == 0) Rcpp::checkUserInterrupt();

    // Choose the move type, as sample.int(3, 1).
    int move = (int) (R_unif_index(3.0) + 1);
    out.moves[move - 1][0]++;
    int n_absent = n_pairs - 2 * n_edges;
    new_dag = current;
    int changed[2] = {0, 0};
    int n_changed = 0;
    bool valid = false;
    double log_q = 0.0;

    if (move == 1) {
      // Addition. Draw a non-adjacent ordered pair by rejection, with
      // sample.int(p, 2) of R. After 64 failures, draw uniformly from the
      // list of non-adjacent pairs.
      int i1 = -1, i2 = -1;
      if (n_absent > 0) {
        bool found = false;
        for (int tries = 0; tries < 64 && !found; ++tries) {
          int nn = p;
          int x_last = nn - 1;
          int j1 = (int) R_unif_index(nn);
          int a1 = j1;
          int replaced = x_last;
          nn--;
          int j2 = (int) R_unif_index(nn);
          int a2 = (j2 == j1) ? replaced : j2;
          if (current[a1 + p * a2] == 0 && current[a2 + p * a1] == 0) {
            i1 = a1;
            i2 = a2;
            found = true;
          }
        }
        if (!found) {
          std::vector<int> ci, cj;
          for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i)
              if (i != j && current[i + p * j] == 0 && current[j + p * i] == 0) {
                ci.push_back(i);
                cj.push_back(j);
              }
          int pick = (int) (R_unif_index((double) ci.size()) + 1) - 1;
          i1 = ci[pick];
          i2 = cj[pick];
        }
        new_dag[i1 + p * i2] = 1;
        if (is_acyclic_core(new_dag, p)) {
          changed[0] = i2;
          n_changed = 1;
          valid = true;
          log_q = std::log((double) n_absent) - std::log((double) (n_edges + 1));
        }
      }
    } else {
      // Deletion or reversal of a uniformly chosen edge, as
      // which(A == 1, arr.ind = TRUE)[sample.int(n_edges, 1), ].
      if (n_edges > 0) {
        int pick = (int) (R_unif_index((double) n_edges) + 1);
        int cnt = 0, i1 = -1, i2 = -1;
        for (int j = 0; j < p && i1 < 0; ++j)
          for (int i = 0; i < p; ++i)
            if (current[i + p * j] == 1 && ++cnt == pick) {
              i1 = i;
              i2 = j;
              break;
            }
        if (move == 2) {
          new_dag[i1 + p * i2] = 0;
          changed[0] = i2;
          n_changed = 1;
          valid = true;
          log_q = std::log((double) n_edges) - std::log((double) (n_absent + 2));
        } else {
          new_dag[i1 + p * i2] = 0;
          new_dag[i2 + p * i1] = 1;
          if (is_acyclic_core(new_dag, p)) {
            changed[0] = i1;
            changed[1] = i2;
            n_changed = 2;
            valid = true;
            log_q = 0.0;
          }
        }
      }
    }

    if (valid) {
      // Score the changed nodes and apply the acceptance rule.
      new_per = per;
      for (int c = 0; c < n_changed; ++c) {
        int k = changed[c];
        new_per[k] = node_score_core(X, n, p, k, parents_of(new_dag, p, k), alpha, d);
        if (!std::isfinite(new_per[k]))
          Rcpp::stop("A score is not finite. "
                     "Check the data for linearly dependent columns.");
      }
      double new_total = r_sum(new_per.data(), p);
      double log_acc = (new_total - total) + log_q;
      if (!std::isnan(log_acc) && std::log(Rf_runif(0.0, 1.0)) < log_acc) {
        current = new_dag;
        per = new_per;
        total = new_total;
        n_edges = (move == 1) ? n_edges + 1 : (move == 2) ? n_edges - 1 : n_edges;
        n_accept++;
        out.moves[move - 1][1]++;
        if (total > out.best_score) {
          out.best_dag = current;
          out.best_score = total;
        }
      }
    } else {
      n_invalid++;
    }

    // Keep every thin-th state after the burn-in.
    if (iter > burn_in && ((iter - burn_in) % thin) == 0 &&
        (int) out.draw_dag.size() < n_ret) {
      std::vector<int> pos = positions_of(current);
      std::string key = key_of_positions(pos);
      std::unordered_map<std::string, int>::iterator it = seen.find(key);
      int u;
      if (it == seen.end()) {
        u = (int) out.dags.size();
        seen[key] = u;
        out.dags.push_back(pos);
      } else {
        u = it->second;
      }
      out.draw_dag.push_back(u + 1);
      out.draw_iter.push_back(iter);
      out.draw_score.push_back(total);
      out.draw_edges.push_back(n_edges);
    }

    // Progress report.
    if (verbose && iter % report_every == 0) {
      Rprintf("  iteration %d of %d: acceptance %.3f, edges %d, score %.2f\n",
              iter, n_iter, (double) n_accept / iter, n_edges, total);
    }
  }
  out.n_accept = n_accept;
  out.n_invalid = n_invalid;
  return out;
}

// =============================================================================
// 6.  SHAPIRO-WILK TEST
// =============================================================================
// Algorithm AS R94 as in R's file src/library/stats/src/swilk.c.
// Copyright (C) the R Core Team. Distributed under the GNU General Public
// License, version 2 or later.

static double sw_poly(const double* cc, int nord, double x) {
  double p, ret_val;
  ret_val = cc[0];
  if (nord > 1) {
    p = x * cc[nord - 1];
    for (int j = nord - 2; j > 0; j--) p = (p + cc[j]) * x;
    ret_val += p;
  }
  return ret_val;
}

static void swilk(double* x, int n, double* w, double* pw, int* ifault) {
  int nn2 = n / 2;
  std::vector<double> a(nn2 + 1); /* 1-based */
  double small = 1e-19;
  double g[2] = { -2.273, .459 };
  double c1[6] = { 0., .221157, -.147981, -2.07119, 4.434685, -2.706056 };
  double c2[6] = { 0., .042981, -.293762, -1.752461, 5.682633, -3.582633 };
  double c3[4] = { .544, -.39978, .025054, -6.714e-4 };
  double c4[4] = { 1.3822, -.77857, .062767, -.0020322 };
  double c5[4] = { -1.5861, -.31082, -.083751, .0038915 };
  double c6[3] = { -.4803, -.082676, .0030302 };
  int i, j, i1;
  double ssassx, summ2, ssumm2, gamma, range;
  double a1, a2, an, m, s, sa, xi, sx, xx, y, w1;
  double fac, asa, an25, ssa, sax, rsn, ssx, xsx;

  *pw = 1.;
  if (n < 3) { *ifault = 1; return; }
  an = (double) n;
  if (n == 3) {
    a[1] = 0.70710678;
  } else {
    an25 = an + .25;
    summ2 = 0.;
    for (i = 1; i <= nn2; i++) {
      a[i] = Rf_qnorm5((i - 0.375) / an25, 0., 1., 1, 0);
      double r__1 = a[i];
      summ2 += r__1 * r__1;
    }
    summ2 *= 2.;
    ssumm2 = std::sqrt(summ2);
    rsn = 1. / std::sqrt(an);
    a1 = sw_poly(c1, 6, rsn) - a[1] / ssumm2;
    if (n > 5) {
      i1 = 3;
      a2 = -a[2] / ssumm2 + sw_poly(c2, 6, rsn);
      fac = std::sqrt((summ2 - 2. * (a[1] * a[1]) - 2. * (a[2] * a[2]))
                      / (1. - 2. * (a1 * a1) - 2. * (a2 * a2)));
      a[2] = a2;
    } else {
      i1 = 2;
      fac = std::sqrt((summ2 - 2. * (a[1] * a[1])) / (1. - 2. * (a1 * a1)));
    }
    a[1] = a1;
    for (i = i1; i <= nn2; i++) a[i] /= -fac;
  }
  range = x[n - 1] - x[0];
  if (range < small) { *ifault = 6; return; }
  *ifault = 0;
  xx = x[0] / range;
  sx = xx;
  sa = -a[1];
  for (i = 1, j = n - 1; i < n; j--) {
    xi = x[i] / range;
    if (xx - xi > small) *ifault = 7;
    sx += xi;
    i++;
    if (i != j) sa += Rf_sign((double) (i - j)) * a[std::min(i, j)];
    xx = xi;
  }
  if (n > 5000) *ifault = 2;
  sa /= n;
  sx /= n;
  ssa = ssx = sax = 0.;
  for (i = 0, j = n - 1; i < n; i++, j--) {
    if (i != j) asa = Rf_sign((double) (i - j)) * a[1 + std::min(i, j)] - sa; else asa = -sa;
    xsx = x[i] / range - sx;
    ssa += asa * asa;
    ssx += xsx * xsx;
    sax += asa * xsx;
  }
  ssassx = std::sqrt(ssa * ssx);
  w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
  *w = 1. - w1;
  if (n == 3) {
    double pi6 = 1.90985931710274, stqr = 1.04719755119660;
    *pw = pi6 * (std::asin(std::sqrt(*w)) - stqr);
    if (*pw < 0.) *pw = 0.;
    return;
  }
  y = std::log(w1);
  xx = std::log(an);
  if (n <= 11) {
    gamma = sw_poly(g, 2, an);
    if (y >= gamma) { *pw = 1e-99; return; }
    y = -std::log(gamma - y);
    m = sw_poly(c3, 4, an);
    s = std::exp(sw_poly(c4, 4, an));
  } else {
    m = sw_poly(c5, 4, xx);
    s = std::exp(sw_poly(c6, 3, xx));
  }
  *pw = Rf_pnorm5(y, m, s, 0, 0);
}

// The p-value of stats::shapiro.test(x). The value 1 is returned in the
// cases where shapiro.test() stops with an error (fewer than 3 or more
// than 5000 values, or a range of zero).
static double shapiro_pvalue_core(std::vector<double> x) {
  int n = (int) x.size();
  if (n < 3 || n > 5000) return 1.0;
  std::sort(x.begin(), x.end());
  double rng = x[n - 1] - x[0];
  if (rng < 1e-10) return 1.0;
  double W = 0.0, pw = 1.0;
  int ifault = 0;
  swilk(x.data(), n, &W, &pw, &ifault);
  if (ifault > 0 && ifault != 7) return 1.0;
  return pw;
}

// =============================================================================
// 7.  STEP 1: THE NON-GAUSSIAN NODES OF A DAG
// =============================================================================
// Node j is declared non-Gaussian when the Shapiro-Wilk p-value of its L1
// residuals is below alpha_sw / p. With more than 5000 residuals, the test
// uses the deterministic subsample res[round(seq(1, n, length.out = 5000))],
// because shapiro.test() accepts at most 5000 values.

static std::vector<int> nongaussian_core(const double* X, int n, int p,
                                         const std::vector<int>& dag, double alpha_sw,
                                         std::vector<double>* pvalues) {
  std::vector<int> out;
  double thr = alpha_sw / p;
  if (pvalues != NULL) pvalues->assign(p, 1.0);
  for (int j = 0; j < p; ++j) {
    // L1 residuals of node j.
    std::vector<int> pa = parents_of(dag, p, j);
    L1Fit f = fit_node_l1_core(X, n, p, j, pa);
    std::vector<double> res = f.resid;
    // Deterministic subsample of 5000 residuals.
    if ((int) res.size() > 5000) {
      int len = (int) res.size();
      std::vector<double> sub(5000);
      double by = (double) (len - 1) / 4999.0;
      for (int k = 0; k < 5000; ++k) {
        double pos = (k == 4999) ? (double) len : 1.0 + k * by;
        int idx = (int) std::nearbyint(pos);
        sub[k] = res[idx - 1];
      }
      res = sub;
    }
    // Test against the Bonferroni level.
    double pv = ((int) res.size() < 3) ? 1.0 : shapiro_pvalue_core(res);
    if (pvalues != NULL) (*pvalues)[j] = pv;
    if (pv < thr) out.push_back(j);
  }
  return out;
}

// =============================================================================
// 8.  STEP 2: CPDAG, MEEK'S RULES AND THE RPDAG.  STEP 3: THE SHD.
// =============================================================================
// A partially directed graph g uses the same storage as a DAG.
// g[a, b] = 1 and g[b, a] = 0 is the edge a -> b. g[a, b] = g[b, a] = 1 is
// the undirected edge a - b.

// Meek's rules R1 to R4, applied in this order until no rule applies.
static void meek_closure(std::vector<int>& g, int p) {
  std::vector<int> old;
#define G(a, b) g[(a) + p * (b)]
  while (old != g) {
    old = g;
    // Rule 1: a -> b - c with a and c not adjacent gives b -> c.
    {
      std::vector<int> ia, ib;
      for (int b = 0; b < p; ++b)
        for (int a = 0; a < p; ++a)
          if (G(a, b) == 1 && G(b, a) == 0) { ia.push_back(a); ib.push_back(b); }
      for (size_t t = 0; t < ia.size(); ++t) {
        int a = ia[t], b = ib[t];
        for (int c = 0; c < p; ++c)
          if (G(b, c) == 1 && G(c, b) == 1 && G(a, c) == 0 && G(c, a) == 0) {
            G(b, c) = 1;
            G(c, b) = 0;
          }
      }
    }
    // Rule 2: a -> c -> b with a - b gives a -> b.
    {
      std::vector<int> ia, ib;
      for (int b = 0; b < p; ++b)
        for (int a = 0; a < p; ++a)
          if (G(a, b) == 1 && G(b, a) == 1) { ia.push_back(a); ib.push_back(b); }
      for (size_t t = 0; t < ia.size(); ++t) {
        int a = ia[t], b = ib[t];
        bool any = false;
        for (int c = 0; c < p; ++c)
          if (G(a, c) == 1 && G(c, a) == 0 && G(c, b) == 1 && G(b, c) == 0) { any = true; break; }
        if (any) {
          G(a, b) = 1;
          G(b, a) = 0;
        }
      }
    }
    // Rule 3: a - c1 -> b and a - c2 -> b with c1 and c2 not adjacent, and
    // a - b, gives a -> b.
    {
      std::vector<int> ia, ib;
      for (int b = 0; b < p; ++b)
        for (int a = 0; a < p; ++a)
          if (G(a, b) == 1 && G(b, a) == 1) { ia.push_back(a); ib.push_back(b); }
      for (size_t t = 0; t < ia.size(); ++t) {
        int a = ia[t], b = ib[t];
        std::vector<int> C;
        for (int c = 0; c < p; ++c)
          if (G(a, c) == 1 && G(c, a) == 1 && G(c, b) == 1 && G(b, c) == 0) C.push_back(c);
        if (C.size() >= 2) {
          bool nonadj = false;
          for (size_t u = 0; u < C.size() && !nonadj; ++u)
            for (size_t v = 0; v < C.size(); ++v)
              if (u != v && G(C[u], C[v]) + G(C[v], C[u]) == 0) { nonadj = true; break; }
          if (nonadj) {
            G(a, b) = 1;
            G(b, a) = 0;
          }
        }
      }
    }
    // Rule 4: a - c -> d -> b with c and b not adjacent, and a - b, gives
    // a -> b. Meek's rule also requires a and d to be adjacent. This
    // condition is not checked here. The orientation is still implied:
    // b -> a would give either the cycle a -> c -> d -> b -> a or the new
    // v-structure c -> a <- b.
    {
      std::vector<int> ia, ib;
      for (int b = 0; b < p; ++b)
        for (int a = 0; a < p; ++a)
          if (G(a, b) == 1 && G(b, a) == 1) { ia.push_back(a); ib.push_back(b); }
      for (size_t t = 0; t < ia.size(); ++t) {
        int a = ia[t], b = ib[t];
        bool found = false;
        for (int c = 0; c < p && !found; ++c) {
          if (!(G(a, c) == 1 && G(c, a) == 1 && G(c, b) == 0 && G(b, c) == 0)) continue;
          for (int d = 0; d < p; ++d)
            if (G(c, d) == 1 && G(d, c) == 0 && G(d, b) == 1 && G(b, d) == 0) { found = true; break; }
        }
        if (found) G(b, a) = 0;
      }
    }
  }
#undef G
}

// The CPDAG of a DAG: the skeleton, the v-structures kept directed, and
// Meek's rules.
static std::vector<int> cpdag_core(const std::vector<int>& dag, int p) {
  // Skeleton.
  std::vector<int> g((size_t) p * p, 0);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i)
      if (dag[i + p * j] == 1) { g[i + p * j] = 1; g[j + p * i] = 1; }
  // V-structures i -> k <- j with i and j not adjacent.
  for (int k = 0; k < p; ++k) {
    std::vector<int> pa = parents_of(dag, p, k);
    for (size_t u = 0; u < pa.size(); ++u)
      for (size_t v = u + 1; v < pa.size(); ++v) {
        int i = pa[u], j = pa[v];
        if (dag[i + p * j] == 0 && dag[j + p * i] == 0) {
          g[i + p * k] = 1; g[k + p * i] = 0;
          g[j + p * k] = 1; g[k + p * j] = 0;
        }
      }
  }
  // Meek's rules.
  meek_closure(g, p);
  return g;
}

// The RPDAG of (dag, nG), by steps (A) to (C) of the paper. nG is 0-based.
static std::vector<int> rpdag_core(const std::vector<int>& dag, int p,
                                   const std::vector<int>& nG) {
  // Step (A): the CPDAG.
  std::vector<int> cp = cpdag_core(dag, p);
  if (nG.empty()) return cp;
  std::vector<char> inG(p, 0);
  for (size_t k = 0; k < nG.size(); ++k) inG[nG[k]] = 1;
  // Step (B): restore the direction of the DAG on every undirected edge
  // with an endpoint in nG.
  std::vector<int> g(cp);
  int n_orient = 0;
  for (int i = 0; i < p - 1; ++i)
    for (int j = i + 1; j < p; ++j)
      if (cp[i + p * j] == 1 && cp[j + p * i] == 1 && (inG[i] || inG[j])) {
        if (dag[i + p * j] == 1 && dag[j + p * i] == 0) {
          g[i + p * j] = 1; g[j + p * i] = 0; ++n_orient;
        } else if (dag[j + p * i] == 1 && dag[i + p * j] == 0) {
          g[j + p * i] = 1; g[i + p * j] = 0; ++n_orient;
        }
      }
  if (n_orient == 0) return cp;
  // Step (C): Meek's rules.
  meek_closure(g, p);
  return g;
}

// The status of every unordered pair i < j, in the column order of
// upper.tri(): 0 no edge, 1 j -> i, 2 i -> j, 3 undirected.
static std::vector<int> pair_state_core(const std::vector<int>& A, int p) {
  std::vector<int> s;
  s.reserve((size_t) p * (p - 1) / 2);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < j; ++i) s.push_back(2 * A[i + p * j] + A[j + p * i]);
  return s;
}

// SHD: the number of unordered pairs whose status differs.
static int shd_core(const std::vector<int>& A, const std::vector<int>& B, int p) {
  int s = 0;
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < j; ++i)
      if (2 * A[i + p * j] + A[j + p * i] != 2 * B[i + p * j] + B[j + p * i]) ++s;
  return s;
}

// The number of edges of a partially directed graph. An undirected edge
// counts once.
static double n_edges_pdag_core(const std::vector<int>& D, int p) {
  int und = 0, dir = 0;
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i) {
      if (D[i + p * j] == 1 && D[j + p * i] == 1) ++und;
      if (D[i + p * j] == 1 && D[j + p * i] == 0) ++dir;
    }
  return (double) (und / 2) + dir;
}

// =============================================================================
// 9.  STEP 3: THE POINT ESTIMATE
// =============================================================================
// D_hat is the distinct RPDAG with the smallest mean SHD to the M projected
// draws. Ties go to fewer edges, then to the graph visited first. The mean
// SHD of every candidate u is computed from a table of pair statuses,
//   sum_m SHD(D_u, D^(m)) = L M - sum_l W[s_u[l], l],
// where W[k, l] counts the draws whose pair l has status k. The estimate of
// the non-Gaussian set is the most frequent set among the draws with
// D^(m) = D_hat. Ties go to the smaller set, then to the set visited first.

static void medoid_core(const std::vector<int>& D_idx,
                        const std::vector< std::vector<int> >& D_unique, int p,
                        const std::vector<int>& nG_idx, const std::vector<int>& nG_size,
                        int& best_out, int& best_g_out,
                        std::vector<double>& mean_shd, std::vector<double>& edge_cnt) {
  int M = (int) D_idx.size(), U = (int) D_unique.size();
  int L = p * (p - 1) / 2;
  // Number of draws of each distinct RPDAG.
  std::vector<double> counts(U, 0.0);
  for (int m = 0; m < M; ++m) counts[D_idx[m] - 1] += 1.0;
  // Pair statuses and edge counts of the distinct RPDAGs.
  std::vector< std::vector<int> > S(U);
  edge_cnt.assign(U, 0.0);
  for (int u = 0; u < U; ++u) {
    S[u] = pair_state_core(D_unique[u], p);
    edge_cnt[u] = n_edges_pdag_core(D_unique[u], p);
  }
  // Table W of pair statuses over the draws.
  std::vector<double> W(4 * (size_t) L, 0.0);
  for (int u = 0; u < U; ++u)
    for (int l = 0; l < L; ++l) W[S[u][l] + 4 * (size_t) l] += counts[u];
  // Mean SHD of every candidate.
  mean_shd.assign(U, 0.0);
  for (int u = 0; u < U; ++u) {
    long double agree = 0.0L;
    for (int l = 0; l < L; ++l) agree += W[S[u][l] + 4 * (size_t) l];
    mean_shd[u] = ((double) L * (double) M - (double) agree) / (double) M;
  }
  // Smallest mean SHD, then fewest edges, then first visited.
  int best = 0;
  for (int u = 1; u < U; ++u) {
    if (mean_shd[u] < mean_shd[best] ||
        (mean_shd[u] == mean_shd[best] && edge_cnt[u] < edge_cnt[best])) best = u;
  }
  // Most frequent non-Gaussian set among the draws with D^(m) = D_hat.
  int G = 0;
  for (int m = 0; m < M; ++m) if (nG_idx[m] > G) G = nG_idx[m];
  std::vector<int> tab(G + 1, 0), first_pos(G + 1, INT_MAX);
  for (int m = 0; m < M; ++m) {
    if (D_idx[m] - 1 != best) continue;
    int g = nG_idx[m];
    tab[g]++;
    if (first_pos[g] == INT_MAX) first_pos[g] = m;
  }
  int max_tab = 0;
  for (int g = 1; g <= G; ++g) if (tab[g] > max_tab) max_tab = tab[g];
  int best_g = -1;
  for (int g = 1; g <= G; ++g) {
    if (tab[g] != max_tab) continue;
    if (best_g < 0 ||
        nG_size[g - 1] < nG_size[best_g - 1] ||
        (nG_size[g - 1] == nG_size[best_g - 1] && first_pos[g] < first_pos[best_g])) best_g = g;
  }
  best_out = best + 1;
  best_g_out = best_g;
}

// =============================================================================
// 10.  FUNCTIONS CALLED FROM R
// =============================================================================
// Every function checks its arguments. Invalid input gives an error and
// never an access outside the arrays. Only the sampler uses
// random numbers. The other functions are exported with rng = false, so
// they neither read nor change the state of the random number generator.

using namespace Rcpp;

// Check a data matrix: at least 2 columns, more rows than columns, finite.
static void check_data_matrix(const NumericMatrix& X) {
  if (X.ncol() < 2) stop("The data must have at least 2 columns.");
  if (X.nrow() < 3 || X.nrow() <= X.ncol())
    stop("The data must have at least 3 rows and more rows than columns.");
  for (R_xlen_t k = 0; k < X.size(); ++k)
    if (!std::isfinite(X[k])) stop("The data must contain only finite values.");
}

// Check a p x p 0/1 matrix with a zero diagonal and copy it.
static std::vector<int> check_adjacency(const IntegerMatrix& A, int p, const char* what) {
  if (A.nrow() != p || A.ncol() != p)
    stop("The %s must be a %d x %d matrix.", what, p, p);
  std::vector<int> v(A.begin(), A.end());
  for (size_t k = 0; k < v.size(); ++k)
    if (v[k] != 0 && v[k] != 1) stop("The %s must contain only 0 and 1.", what);
  for (int i = 0; i < p; ++i)
    if (v[i + (size_t) p * i] != 0) stop("The %s must have a zero diagonal.", what);
  return v;
}

// Check a DAG: a 0/1 matrix with a zero diagonal and no directed cycle.
static std::vector<int> check_dag(const IntegerMatrix& A, int p) {
  std::vector<int> v = check_adjacency(A, p, "DAG");
  if (!is_acyclic_core(v, p)) stop("The DAG must not contain a directed cycle.");
  return v;
}

// Turn 1-based edge positions into a p x p 0/1 vector.
static std::vector<int> graph_from_positions(const IntegerVector& pos, int p) {
  std::vector<int> A((size_t) p * p, 0);
  for (R_xlen_t k = 0; k < pos.size(); ++k) {
    int q = pos[k];
    if (q == NA_INTEGER || q < 1 || q > p * p) stop("Internal error: edge position out of range.");
    if ((q - 1) % (p + 1) == 0) stop("Internal error: edge on the diagonal.");
    A[q - 1] = 1;
  }
  return A;
}

// Turn a p x p vector into 1-based edge positions.
static IntegerVector positions_to_r(const std::vector<int>& A) {
  std::vector<int> pos = positions_of(A);
  IntegerVector out(pos.size());
  for (size_t k = 0; k < pos.size(); ++k) out[k] = pos[k] + 1;
  return out;
}

static IntegerMatrix vec_to_mat(const std::vector<int>& v, int p) {
  IntegerMatrix A(p, p);
  std::copy(v.begin(), v.end(), A.begin());
  return A;
}

// TRUE if the 0/1 matrix A has no directed cycle.
// [[Rcpp::export(rng = false)]]
bool lewis_is_acyclic_cpp(IntegerMatrix A) {
  int p = A.nrow();
  std::vector<int> v = check_adjacency(A, p, "graph");
  return is_acyclic_core(v, p);
}

// The score s_n of one DAG and its node terms.
// [[Rcpp::export(rng = false)]]
List lewis_score_dag_cpp(NumericMatrix X, IntegerMatrix dag, double alpha, double d) {
  check_data_matrix(X);
  int n = X.nrow(), p = X.ncol();
  if (!std::isfinite(alpha) || !std::isfinite(d) || d < 0) stop("Invalid prior constants.");
  std::vector<int> A = check_dag(dag, p);
  NumericVector per(p);
  for (int j = 0; j < p; ++j)
    per[j] = node_score_core(X.begin(), n, p, j, parents_of(A, p, j), alpha, d);
  std::vector<double> pv(per.begin(), per.end());
  return List::create(_["total"] = r_sum(pv.data(), p), _["per_node"] = per);
}

// The sampler. The distinct DAGs are returned as 1-based edge positions.
// [[Rcpp::export]]
List lewis_mcmc_cpp(NumericMatrix X, IntegerMatrix start, int n_iter, int burn_in,
                    int thin, double alpha, double d, bool verbose) {
  check_data_matrix(X);
  int n = X.nrow(), p = X.ncol();
  std::vector<int> A = check_dag(start, p);
  if (n_iter < 1) stop("The number of iterations must be positive.");
  if (burn_in < 0 || burn_in >= n_iter) stop("The burn-in must be smaller than the number of iterations.");
  if (thin < 1) stop("The thinning interval must be positive.");
  if (!std::isfinite(alpha) || !std::isfinite(d) || d < 0) stop("Invalid prior constants.");

  ChainResult res = mcmc_core(X.begin(), n, p, A, n_iter, burn_in, thin, alpha, d, verbose);

  // Distinct DAGs as 1-based edge positions.
  List dags(res.dags.size());
  for (size_t u = 0; u < res.dags.size(); ++u) {
    IntegerVector pos(res.dags[u].size());
    for (size_t k = 0; k < res.dags[u].size(); ++k) pos[k] = res.dags[u][k] + 1;
    dags[u] = pos;
  }
  // Move counts.
  IntegerMatrix moves(3, 2);
  for (int a = 0; a < 3; ++a) for (int b = 0; b < 2; ++b) moves(a, b) = res.moves[a][b];
  moves.attr("dimnames") = List::create(CharacterVector::create("add", "delete", "reverse"),
                                        CharacterVector::create("proposed", "accepted"));
  return List::create(_["dags"] = dags,
                      _["draw_dag"] = IntegerVector(res.draw_dag.begin(), res.draw_dag.end()),
                      _["iteration"] = IntegerVector(res.draw_iter.begin(), res.draw_iter.end()),
                      _["log_score"] = NumericVector(res.draw_score.begin(), res.draw_score.end()),
                      _["n_edges"] = IntegerVector(res.draw_edges.begin(), res.draw_edges.end()),
                      _["best_dag"] = vec_to_mat(res.best_dag, p),
                      _["best_score"] = res.best_score,
                      _["acceptance_rate"] = (double) res.n_accept / n_iter,
                      _["invalid_rate"] = (double) res.n_invalid / n_iter,
                      _["moves"] = moves);
}

// Step 1 for one DAG: the non-Gaussian nodes (1-based) and the p-values.
// [[Rcpp::export(rng = false)]]
List lewis_nongaussian_cpp(NumericMatrix X, IntegerMatrix dag, double alpha_sw) {
  check_data_matrix(X);
  int n = X.nrow(), p = X.ncol();
  if (!(alpha_sw > 0 && alpha_sw < 1)) stop("alpha_sw must lie strictly between 0 and 1.");
  std::vector<int> A = check_dag(dag, p);
  std::vector<double> pv;
  std::vector<int> nG = nongaussian_core(X.begin(), n, p, A, alpha_sw, &pv);
  IntegerVector out(nG.size());
  for (size_t k = 0; k < nG.size(); ++k) out[k] = nG[k] + 1;
  return List::create(_["nongaussian"] = out,
                      _["p_value"] = NumericVector(pv.begin(), pv.end()));
}

// Step 2 for one DAG and one non-Gaussian set (1-based).
// [[Rcpp::export(rng = false)]]
IntegerMatrix lewis_rpdag_cpp(IntegerMatrix dag, IntegerVector nG) {
  int p = dag.nrow();
  std::vector<int> A = check_dag(dag, p);
  std::vector<int> g(nG.size());
  for (R_xlen_t k = 0; k < nG.size(); ++k) {
    if (nG[k] == NA_INTEGER || nG[k] < 1 || nG[k] > p)
      stop("The non-Gaussian nodes must be indices between 1 and %d.", p);
    g[k] = nG[k] - 1;
  }
  return vec_to_mat(rpdag_core(A, p, g), p);
}

// Steps 1 and 2 for every distinct DAG of a chain. The DAGs are given as
// 1-based edge positions. The RPDAGs are returned in the same form.
// [[Rcpp::export(rng = false)]]
List lewis_project_cpp(NumericMatrix X, List dags, double alpha_sw) {
  check_data_matrix(X);
  int n = X.nrow(), p = X.ncol();
  if (!(alpha_sw > 0 && alpha_sw < 1)) stop("alpha_sw must lie strictly between 0 and 1.");
  int U = dags.size();
  List u_nG(U), u_D(U);
  for (int u = 0; u < U; ++u) {
    // Rebuild the DAG.
    IntegerVector pos = dags[u];
    std::vector<int> A = graph_from_positions(pos, p);
    if (!is_acyclic_core(A, p)) stop("Internal error: a sampled graph has a cycle.");
    // Step 1.
    std::vector<int> nG = nongaussian_core(X.begin(), n, p, A, alpha_sw, NULL);
    IntegerVector nGr(nG.size());
    for (size_t k = 0; k < nG.size(); ++k) nGr[k] = nG[k] + 1;
    u_nG[u] = nGr;
    // Step 2.
    u_D[u] = positions_to_r(rpdag_core(A, p, nG));
    if (u % 100 == 99) checkUserInterrupt();
  }
  return List::create(_["nongaussian"] = u_nG, _["rpdag"] = u_D);
}

// SHD between two partially directed graphs.
// [[Rcpp::export(rng = false)]]
int lewis_shd_cpp(IntegerMatrix A, IntegerMatrix B) {
  int p = A.nrow();
  std::vector<int> a = check_adjacency(A, p, "first graph");
  std::vector<int> b = check_adjacency(B, p, "second graph");
  return shd_core(a, b, p);
}

// SHD from one graph to each graph of a list, all given as edge positions.
// [[Rcpp::export(rng = false)]]
IntegerVector lewis_shd_list_cpp(IntegerVector ref, List graphs, int p) {
  if (p < 2) stop("Internal error: p must be at least 2.");
  std::vector<int> R = graph_from_positions(ref, p);
  int U = graphs.size();
  IntegerVector out(U);
  for (int u = 0; u < U; ++u) {
    IntegerVector pos = graphs[u];
    out[u] = shd_core(R, graph_from_positions(pos, p), p);
  }
  return out;
}

// Weighted frequency of each directed edge and each undirected edge over a
// list of graphs given as edge positions.
// [[Rcpp::export(rng = false)]]
List lewis_edge_freq_cpp(List graphs, NumericVector w, int p) {
  if (p < 2) stop("Internal error: p must be at least 2.");
  int U = graphs.size();
  if (w.size() != U) stop("Internal error: one weight per graph is needed.");
  NumericMatrix dir(p, p), und(p, p);
  for (int u = 0; u < U; ++u) {
    IntegerVector pos = graphs[u];
    std::vector<int> D = graph_from_positions(pos, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < p; ++i) {
        if (D[i + p * j] == 1 && D[j + p * i] == 0) dir(i, j) += w[u];
        if (D[i + p * j] == 1 && D[j + p * i] == 1) und(i, j) += w[u];
      }
  }
  return List::create(_["directed"] = dir, _["undirected"] = und);
}

// Step 3: the SHD medoid and its non-Gaussian set.
//   D_idx    index (1-based) of the RPDAG of each draw among D_unique
//   D_unique distinct RPDAGs as edge positions, in order of first visit
//   nG_idx   index (1-based) of the non-Gaussian set of each draw
//   nG_size  size of each distinct non-Gaussian set
// [[Rcpp::export(rng = false)]]
List lewis_medoid_cpp(IntegerVector D_idx, List D_unique, IntegerVector nG_idx,
                      IntegerVector nG_size, int p) {
  int M = D_idx.size(), U = D_unique.size();
  if (U == 0 || M == 0) stop("Internal error: no draws to summarize.");
  if (nG_idx.size() != M) stop("Internal error: one non-Gaussian set per draw is needed.");
  if (p < 2) stop("Internal error: p must be at least 2.");
  // Check the indices.
  for (int m = 0; m < M; ++m) {
    if (D_idx[m] == NA_INTEGER || D_idx[m] < 1 || D_idx[m] > U)
      stop("Internal error: graph index out of range.");
    if (nG_idx[m] == NA_INTEGER || nG_idx[m] < 1 || nG_idx[m] > nG_size.size())
      stop("Internal error: set index out of range.");
  }
  // Rebuild the distinct RPDAGs.
  std::vector< std::vector<int> > Du(U);
  for (int u = 0; u < U; ++u) {
    IntegerVector pos = D_unique[u];
    Du[u] = graph_from_positions(pos, p);
  }
  std::vector<int> di(D_idx.begin(), D_idx.end()), gi(nG_idx.begin(), nG_idx.end());
  std::vector<int> gs(nG_size.begin(), nG_size.end());
  int best = 0, best_g = 0;
  std::vector<double> mean_shd, edge_cnt;
  medoid_core(di, Du, p, gi, gs, best, best_g, mean_shd, edge_cnt);
  return List::create(_["best_D"] = best,
                      _["mean_shd"] = NumericVector(mean_shd.begin(), mean_shd.end()),
                      _["edge_counts"] = NumericVector(edge_cnt.begin(), edge_cnt.end()),
                      _["best_nG"] = best_g);
}
