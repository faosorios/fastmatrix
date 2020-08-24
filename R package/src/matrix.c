/* $ID: matrix.c, last updated 2020-08-08, F.Osorio */

#include "base.h"
#include "utils.h"
#include "matrix.h"

void
equilibrate(double *x, int *ldx, int *nrow, int *ncol, double *scales, double *condition, int *job)
{ /* columns equilibration of a rectangular matrix */
  int info = 0;

  F77_CALL(equilibrate_cols)(x, ldx, nrow, ncol, scales, condition, job, &info);
  if (info)
    error("equilibrate_cols gave code %d", info);
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
