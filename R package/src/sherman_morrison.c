/* $ID: sherman_morrison.c, last updated 2020-09-04, F.Osorio */

#include "fastmatrix.h"

void
sherman_morrison(double *a, int *lda, int *n, double *b, double *d)
{ /* Sherman–Morrison formula */
  char *notrans = "N", *trans = "T";
  int p = *n;
  double alpha = -1.0, dot, *u = NULL, *v = NULL;

  /* initializing */
  u = (double *) Calloc(p, double);
  v = (double *) Calloc(p, double);

  /* updating b and d */
  BLAS2_gemv(1.0, a, *lda, p, p, notrans, b, 1, 0.0, u, 1);
  dot = BLAS1_dot_product(d, 1, u, 1, p);
  BLAS2_gemv(1.0, a, *lda, p, p, trans, d, 1, 0.0, v, 1);

  /* applying rank-1 update */
  alpha /= 1.0 + dot;
  BLAS2_ger(alpha, a, *lda, p, p, u, 1, v, 1);

  Free(u); Free(v);
}
