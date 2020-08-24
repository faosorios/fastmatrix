/* $ID: vector.c, last updated 2020-08-21, F.Osorio */

#include "base.h"
#include "vector.h"

void
normalize_vec(double *x, int inc, int n)
{ /* x <- x / sqrt(sum(x * x)) */
  double div = 1.0, length;

  length = F77_CALL(dnrm2)(&n, x, &inc);
  div /= length;
  F77_CALL(dscal)(&n, &div, x, &inc);
}
