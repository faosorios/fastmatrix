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
FM_cov2cor(double *cov, int p)
{ /* scales a 'covariance' matrix into the corresponding correlation matrix */
  double *s;

  s = (double *) Calloc(p, double);

  for (int i = 0; i < p; i++)
    s[i] = cov[i * (p + 1)];

  for (int i = 0; i < p; i++) {
    cov[i * (p + 1)] = 1.0;
    for (int j = i + 1; j < p; j++) {
      *(cov + i + j * p) = *(cov + i + j * p) / sqrt(s[i] * s[j]);
      *(cov + j + i * p) = *(cov + i + j * p);
    }
  }
  Free(s);
}

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
