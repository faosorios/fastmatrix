/* ID: mardia.c, last updated 2025-10-22, F.Osorio */

#include "fastmatrix.h"

void
skewness_and_kurtosis(double *x, int *n, int *p, double *center, double *Scatter, double *stats, int *task)
{ /* wrapper to 'FM_skewness_and_kurtosis' */
  int do_skewness = *task;

  FM_skewness_and_kurtosis(x, *n, *p, center, Scatter, stats, do_skewness);
}

void
mardia_stat(double *x, int *nobs, int *vars, double *center, double *Scatter, double *coef, double *stats)
{ /* computes Mardia test for multivariate normality based on multivariate 
   * skewness and kurtosis coefficients
   * Biometrika 57, 519-530, 1970. doi: 10.1093/biomet/57.3.519
   * Sankhya B 36, 115-128, 1974.  url: www.jstor.org/stable/25051892 */
  int n = *nobs, p = *vars, task = 1;
  double b1, b2, ex, z1, z2;

  /* compute multivariate skewness and kurtosis */
  FM_skewness_and_kurtosis(x, n, p, center, Scatter, coef, task);
  b1 = coef[0];
  b2 = coef[1];

  /* Mardia's test statistics */
  ex = (double) p * (p + 2.); /* expectation of b2 */
  z1 = (double) n * b1 / 6.;
  z2 = (b2 - ex) / sqrt(8 * ex / n);

  /* saving results */
  stats[0] = z1;
  stats[1] = z2;
}
