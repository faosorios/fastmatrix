/* $ID: bezier.c, last updated 2024-09-03, F.Osorio */

#include "fastmatrix.h"

void
bezier_smoother(double *x, double *y, int *n, double *grid, int *ngrid, double *xgrid, double *ygrid)
{ /* computation of Bezier smoother */
  double s, *work;

  if (*ngrid <= 0) /* quick return if possible */
    return;

  work = (double *) R_Calloc(2, double);
  for (int i = 0; i < *ngrid; i++) {
    s = grid[i];
    F77_CALL(decasteljau)(x, y, n, &s, work);
    xgrid[i] = work[0];
    ygrid[i] = work[1];
  }
  R_Free(work);
}
