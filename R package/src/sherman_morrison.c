/* $ID: sherman_morrison.c, last updated 2020-08-23, F.Osorio */

#include "base.h"
#include "sherman_morrison.h"

void
sherman_morrison(double *a, int *lda, int *n, double *b, double *d)
{ /* Sherman–Morrison formula */
  char *notrans = "N", *trans = "T";
  int inc = 1, p = *n;
  double alpha = -1.0, dot, one = 1.0, zero = 0.0, *u = NULL, *v = NULL;

  /* initializing */
  u = (double *) Calloc(p, double);
  v = (double *) Calloc(p, double);

  /* updating b and d */
  F77_CALL(dgemv)(notrans, &p, &p, &one, a, lda, b, &inc, &zero, u, &inc);
  dot = F77_CALL(ddot)(&p, d, &inc, u, &inc);
  F77_CALL(dgemv)(trans, &p, &p, &one, a, lda, d, &inc, &zero, v, &inc);

  /* applying rank-1 update */
  alpha /= 1.0 + dot;
  F77_CALL(dger)(&p, &p, &alpha, u, &inc, v, &inc, a, lda);

  Free(u); Free(v);
}
