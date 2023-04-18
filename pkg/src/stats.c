/* ID: stats.c, last updated 2023-04-17, F.Osorio */

#include "fastmatrix.h"

void
central_moments(double *x, int *nobs, double *mean, double *m2, double *m3, double *m4)
{ /* wrapper for 'FM_moments' */
  int n = *nobs;

  FM_moments(x, n, mean, m2, m3, m4);
}
