/* ID: cov_weighted.c, last updated 2021-03-03, F.Osorio */

#include "fastmatrix.h"

void
cov_weighted(double *x, int *nobs, int *vars, double *weights, double *mean, double *cov)
{ /* wrapper for 'FM_mean_and_Scatter' */
  int n = *nobs, p = *vars;

  FM_center_and_Scatter(x, n, p, weights, mean, cov);
}

void
cov_MSSD(double *x, int *nobs, int *vars, double *mean, double *cov)
{ /* wrapper for 'FM_cov_MSSD' */
  int n = *nobs, p = *vars;

  FM_cov_MSSD(x, n, p, mean, cov);
}

void
central_moments(double *x, int *nobs, double *mean, double *var, double *m3, double *m4)
{ /* wrapper for 'FM_moments' */
  int n = *nobs;

  FM_moments(x, n, mean, var, m3, m4);
}
