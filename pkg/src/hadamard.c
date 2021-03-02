/* $ID: hadamard.c, last updated 2021-02-14, F.Osorio */

#include "fastmatrix.h"

void
hadamard_prod(double *x, double *y, int *n, double *prod)
{ /* returns the element-wise product between vectors 'x' and 'y'
   * using unrolled loops */
  int m = *n % 8;

  if (*n <= 0) /* quick return if possible */
    return;

  for (int i = 0; i < m; i++)
    prod[i] = x[i] * y[i];

  for (int i = m; i + 7 < *n; i += 8) {
    prod[i] = x[i] * y[i];
    prod[i + 1] = x[i + 1] * y[i + 1];
    prod[i + 2] = x[i + 2] * y[i + 2];
    prod[i + 3] = x[i + 3] * y[i + 3];
    prod[i + 4] = x[i + 4] * y[i + 4];
    prod[i + 5] = x[i + 5] * y[i + 5];
    prod[i + 6] = x[i + 6] * y[i + 6];
    prod[i + 7] = x[i + 7] * y[i + 7];
  }
}
