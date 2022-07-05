/* ID: stats.c, last updated 2022-06-01, F.Osorio */

#include "fastmatrix.h"

void
central_moments(double *x, int *nobs, double *mean, double *var, double *m3, double *m4)
{ /* wrapper for 'FM_moments' */
  int n = *nobs;

  FM_moments(x, n, mean, var, m3, m4);
}
