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
