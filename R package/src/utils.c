/* $ID: utils.c, last updated 2020-09-03, F.Osorio */

#include "fastmatrix.h"

/* operations on vectors */

double
norm_sqr(double *x, int inc, int n)
{ /* sum(x * x) */
  double length;

  length = BLAS1_norm_two(x, inc, n);
  return R_pow_di(length, 2);
}

void
normalize_vec(double *x, int inc, int n)
{ /* x <- x / sqrt(sum(x * x)) */
  double div = 1.0, length;

  length = BLAS1_norm_two(x, inc, n);
  div /= length;
  BLAS1_scale(div, x, inc, n);
}

/* matrix operations */

void
mat2vech(double *x, int *ldx, int *n, double *y)
{ /* y <- vech(x) */
  int p = *n, k = 0;

  for (int j = 0; j < p; j++) {
    for (int i = j; i < p; i++) {
      y[k] = x[i + j * *ldx];
      k++;
    }
  }
}
