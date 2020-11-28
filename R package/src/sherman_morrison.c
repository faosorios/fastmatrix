/* $ID: sherman_morrison.c, last updated 2020-11-18, F.Osorio */

#include "fastmatrix.h"

void
sherman_morrison(double *a, int *lda, int *n, double *b, double *d, int *inverted)
{ /* Sherman–Morrison formula (to be called by R) */
  FM_sherman_morrison(a, *lda, *n, b, d, *inverted);
}

void
FM_sherman_morrison(double *a, int lda, int n, double *b, double *d, int inverted)
{ /* Sherman–Morrison formula (symbol exported to API) */
  char *notrans = "N", *trans = "T";
  double alpha = -1.0, dot, *u = NULL, *v = NULL;
  int *pivot = NULL;

  /* initializing */
  u = (double *) Calloc(n, double);
  v = (double *) Calloc(n, double);
  pivot = (int *) Calloc(n, int);

  /* invert 'a' matrix in-place (if requested) */
  if (!inverted) {
    lu_dcmp(a, &lda, &n, &n, pivot);
    lu_inverse(a, &lda, &n, pivot);
  }

  /* updating b and d */
  BLAS2_gemv(1.0, a, lda, n, n, notrans, b, 1, 0.0, u, 1);
  dot = BLAS1_dot_product(d, 1, u, 1, n);
  BLAS2_gemv(1.0, a, lda, n, n, trans, d, 1, 0.0, v, 1);

  /* applying rank-1 update */
  alpha /= 1.0 + dot;
  BLAS2_ger(alpha, a, lda, n, n, u, 1, v, 1);

  Free(u); Free(v); Free(pivot);
}
