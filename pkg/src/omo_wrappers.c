/* $ID: omo_wrappers.c, last updated 2022-07-01, F.Osorio */

#include "fastmatrix.h"

/* other matrix operations wrappers */

double
FM_blinf(double *a, int lda, int n, int p, double *x, double *y, int *info)
{ /* returns the value of t(x) %*% a %*% y */
  return F77_CALL(blinf)(a, &lda, &n, &p, x, y, info);
}

double
FM_quadf(double *a, int lda, int n, double *x, int *info)
{ /* returns the value of t(x) %*% a %*% x */
  return F77_CALL(quadf)(a, &lda, &n, x, info);
}

void
FM_murrv(double *y, double *a, int lda, int n, int p, double *x, int *info)
{ /* computes y <- a %*% x */
  F77_CALL(murrv)(y, a, &lda, &n, &p, x, info);
}
