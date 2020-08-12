/* $ID: matrix.c, last updated 2020-08-08, F.Osorio */

#include "base.h"
#include "norms.h"

void
norm_one(double *x, int *inc, int *n, double *value)
{ /* absolute-value norm */
  *value = F77_CALL(dasum)(n, x, inc);
}

void
norm_two(double *x, int *inc, int *n, double *value)
{ /* Euclidean norm */
  *value = F77_CALL(dnrm2)(n, x, inc);
}

void
norm_inf(double *x, int *inc, int *n, double *value)
{ /* infinity norm */
  *value = F77_CALL(dnrminf)(n, x, inc);
}

void
norm_minkowski(double *x, int *inc, int *n, double *p, double *value)
{ /* Minkowski norm */
  *value = F77_CALL(minkowski)(n, x, inc, p);
}
