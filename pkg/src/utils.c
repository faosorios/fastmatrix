/* $ID: utils.c, last updated 2020-09-03, F.Osorio */

#include "fastmatrix.h"

/* matrix operations */

void
FM_centering(double *x, int n, int p, double *center)
{ /* 'x' matrix is overwritten with its centered version */
  for (int i = 0; i < n; i++)
    BLAS1_axpy(-1.0, center, 1, x + i, n, p);
}

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
