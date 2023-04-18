/* $ID: bezier.c, last updated 2021-11-12, F.Osorio */

#include "fastmatrix.h"

void
bezier_smoother(double *x, double *y, int *n, double *grid, int *ngrid, double *xgrid, double *ygrid)
{ /* computation of Bezier smoother */
  double s, *work;

  if (*ngrid <= 0) /* quick return if possible */
    return;

  work = (double *) Calloc(2, double);
  for (int i = 0; i < *ngrid; i++) {
    s = grid[i];
    F77_CALL(decasteljau)(x, y, n, &s, work);
    xgrid[i] = work[0];
    ygrid[i] = work[1];
  }
  Free(work);
}
