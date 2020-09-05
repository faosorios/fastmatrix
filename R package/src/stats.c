/* ID: stats.c, last updated 2020-09-03, F.Osorio */

#include "fastmatrix.h"

void
center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter)
{ /* compute center and Scatter estimates using an online algorithm
   * based on AS 41: Applied Statistics 20, 1971, pp. 206-209 */
  double accum = 0.0, factor, wts, *diff;

  /* initialization */
  diff = (double *) Calloc(p, double);
  Memcpy(center, x, p); /* copying 1st observation */

  /* updating stage */
  for (int i = 1; i < n; i++) {
    wts = weights[i];
    accum += wts;
    factor = wts / accum;
    Memcpy(diff, x + i * p, p);
    BLAS1_axpy(-1.0, x, 1, diff, 1, p);
    BLAS1_axpy(1.0, center, 1, diff, 1, p);
    factor = wts - factor * wts;
    BLAS2_ger(factor, Scatter, p, p, p, diff, 1, diff, 1);
  }
  /* scaling Scatter */
  scale_mat(Scatter, p, 1.0 / n, Scatter, p, p, p);

  Free(diff);
}

void
MSSD(double *x, int n, int p, double *center, double *Scatter)
{ /* compute center and Scatter estimates using the Mean Square Successive Method (MSSD) */
  int accum = 1, nobs = n;
  double *curr, *diff, *prev;

  /* initialization */
  curr = (double *) Calloc(p, double);
  diff = (double *) Calloc(p, double);
  prev = (double *) Calloc(p, double);
  Memcpy(center, x, p); /* copying 1st observation */
  Memcpy(prev, x, p);   /* copying again */

  /* updating stage */
  for (int i = 1; i < nobs; i++) {
    accum++;
    Memcpy(curr, x + i * p, p); /* current observation */
    Memcpy(diff, curr, p);
    BLAS1_axpy(-1.0, prev, 1, curr, 1, p); /* successive difference */
    BLAS2_ger(0.5 / (nobs - 1.0), Scatter, p, p, p, curr, 1, curr, 1);
    Memcpy(prev, diff, p);
    BLAS1_axpy(-1.0, center, 1, diff, 1, p);
    BLAS1_axpy(1.0 / accum, diff, 1, center, 1, p);
  }

  Free(curr); Free(diff); Free(prev);
}
