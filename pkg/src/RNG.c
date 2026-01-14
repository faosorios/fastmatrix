/* $ID: RNG.c, last updated 2025-12-18, F.Osorio */

#include "fastmatrix.h"

/* static functions */
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
