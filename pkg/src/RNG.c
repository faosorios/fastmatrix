/* $ID: RNG.c, last updated 2023-07-22, F.Osorio */

#include "fastmatrix.h"

/* static functions */
static double chi_rand(double);
static void mnorm_rand(double *, int, int);
static void unif_ball_rand(double *, int, int);
static void unif_sphere_rand(double *, int, int);
/* ..end declarations */

void
rng_mnorm(double *y, int *nobs, int *nvar, double *center, double *cov)
{ /* multivariate normal random generation */
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  int info = 0, job = 1, n = *nobs, p = *nvar;

  GetRNGstate();
  FM_chol_decomp(cov, p, p, job, &info);
  if (info)
    error("cholesky factorization in mnorm_rand gave code %d", info);
  mnorm_rand(y, n, p);
  BLAS3_trmm(1.0, cov, p, p, n, side, uplo, trans, diag, y, p);
  for (int i = 0; i < n; i++) {
    BLAS1_axpy(1.0, center, 1, y, 1, p);
    y += p;
  }
  PutRNGstate();
}

void
mnorm_rand(double *y, int n, int p)
{ /* independent standard normal variates */
  int np = n * p;

  const int m = np % 8;

  for (int i = 0; i < m; i++)
    y[i] = norm_rand();

  for (int i = m; i + 7 < np; i += 8) {
    y[i] = norm_rand();
    y[i + 1] = norm_rand();
    y[i + 2] = norm_rand();
    y[i + 3] = norm_rand();
    y[i + 4] = norm_rand();
    y[i + 5] = norm_rand();
    y[i + 6] = norm_rand();
    y[i + 7] = norm_rand();
  }
}

void
rng_sphere(double *y, int *nobs, int *nvar)
{ /* random vector uniformly distributed on the unitary sphere */
  int n = *nobs, p = *nvar;

  GetRNGstate();
  unif_sphere_rand(y, n, p);
  PutRNGstate();
}

void
unif_sphere_rand(double *y, int n, int p)
{ /* RNG uniformly distributed on the unitary sphere */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    FM_normalize(y, 1, p);
    y += p;
  } 
}

void
rng_ball(double *y, int *nobs, int *nvar)
{ /* random vector uniformly distributed in the unitary ball */
  int n = *nobs, p = *nvar;

  GetRNGstate();
  unif_ball_rand(y, n, p);
  PutRNGstate();
}

void
unif_ball_rand(double *y, int n, int p)
{ /* RNG uniformly distributed in the unitary ball */
  double r, u, s;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    s = BLAS1_norm_two(y, 1, p);
    u = unif_rand();
    r = R_pow(u, 1. / p) / s;
    BLAS1_scale(r, y, 1, p);
    y += p;
  } 
}

void rng_chi(int *n, double *x, double *df, int *ndf)
{ /* random variates from the chi distribution */
  int nobs = *n, nd = *ndf;

  GetRNGstate();
  for (int i = 0; i < nobs; i++)
    x[i] = chi_rand(df[i % nd]);
  PutRNGstate();
}

double
chi_rand(double df)
{ /* random number generation from the chi family of distributions with
   * degrees of freedom parameter df >= 1, using the ratio of uniforms
   * method as described in:
   * Monahan, J.F. (1987). ACM Transactions on Mathematical Software 13, 168-172,
   * Monahan, J.F. (1988). ACM Transactions on Mathematical Software 14, 111. */

  /* Constants : */
  const static double S_EXP_M05    = 0.60653065971263342426; /* exp(-1/2)      */
  const static double S_2_EXP_025  = 2.56805083337548278877; /* 2 * exp(1/4)   */
  const static double S_4_EXP_M135 = 1.03696104258356602834; /* 4 * exp(-1.35) */

  double a, b, c, eta, h, r, u, v, z;

  /* Setup: */
  eta = sqrt(df - 1.0);
  a = S_EXP_M05 * (M_SQRT1_2 + eta) / (0.5 + eta);
  b = -S_EXP_M05 * (1.0 - 0.25 / (SQR(eta) + 1.0));
  c = MAX(-eta, b);

  repeat {
    u = unif_rand();
    v = a + (c - a) * unif_rand();
    z = v / u;

    /* Immediate rejection */
    if (z < -eta)
      continue;

    /* Quick acceptance */
    r = 2.5 - SQR(z);
    if (z < 0.0)
      r += z * SQR(z) / (3.0 * (z + eta));
    if (u < r / S_2_EXP_025)
      return z + eta;

    /* Quick rejection */
    if (SQR(z) > S_4_EXP_M135 / u + 1.4)
      continue;

    /* Regular test */
    h = SQR(eta) * log(1.0 + z / eta) - 0.5 * SQR(z) - z * eta;
    if (2 * log(u) < h)
      return z + eta;
  }
}
