/* $ID: matrix.c, last updated 2020-08-08, F.Osorio */

#include "base.h"
#include "matrix.h"

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

void
hadamard_prod(double *x, double *y, int *n, double *prod)
{ /* prod <- x * y */
  for (int i = 0; i < *n; i++)
    *prod++ = *x++ * *y++;
}
