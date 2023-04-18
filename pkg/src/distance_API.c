/* ID: distance_API.c, last updated 2020-10-17, F.Osorio */

#include "fastmatrix.h"

double
FM_pythag(double a, double b)
{ /* finds sqrt(a^2 + b^2) without overflow or destructive underflow
   * Moler, Morrison (1983). IBM Journal of Research and Development 27, 577–581.
   * doi: 10.1147/rd.276.0577 */
  double p, r, s, t, u;

  p = fmax2(fabs(a), fabs(b));
  if (p == 0.0)
    return p;

  r = fmin2(fabs(a), fabs(b) / p);
  r = r * r;

  repeat {
    t = 4.0 + r;
    if (t == 4.0)
      break;
    s = r / t;
    u = 1.0 + 2.0 * s;
    p *= u;
    r *= (s / u) * (s / u);
  }

  return p;
}

double
FM_mahalanobis(double *x, int p, double *center, double *Root)
{ /* Mahalanobis distance (just for a single observation), the argument 'Root'
   * corresponds to the lower triangular factor of the covariance matrix */
  int info = 0;
  double ans, *z;

  /* initialization */
  z = (double *) Calloc(p, double);

  /* computation of Mahalanobis distance */
  Memcpy(z, x, p);
  BLAS1_axpy(-1.0, center, 1, z, 1, p);
  FM_forwardsolve(Root, p, p, z, p, 1, &info);
  if (info) {
    Free(z);
    error("Covariance matrix is possibly singular");
  }
  ans = FM_norm_sqr(z, 1, p);

  Free(z);
  return ans;
}

void
FM_mahal_distances(double *x, int n, int p, double *center, double *cov, int inverted, double *distances)
{ /* returns the Mahalanobis distances (code re-use) */
  char *uplo, *notrans = "N", *diag = "N";
  int info = 0, job = inverted;
  double *z;

  /* attemp to decompose the 'cov' matrix or its inverse */
  FM_chol_decomp(cov, p, p, job, &info);
  if (info)
    error("Covariance matrix is possibly not positive-definite");

  /* invert triangular factor in-place (if requested) */
  if (!inverted) {
    FM_invert_triangular(cov, p, p, job, &info);
    if (info)
      error("Covariance matrix is possibly singular");
  }

  z = (double *) Calloc(p, double);
  uplo = (job) ? "U" : "L";

  /* computation of Mahalanobis distance */
  for (int i = 0; i < n; i++) {
    BLAS1_copy(z, 1, x + i, n, p);
    BLAS1_axpy(-1.0, center, 1, z, 1, p);
    BLAS2_trmv(cov, p, p, uplo, notrans, diag, z, 1);
    distances[i] = FM_norm_sqr(z, 1, p);
  }

  Free(z);
}
