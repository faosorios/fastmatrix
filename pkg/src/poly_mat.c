/* $ID: poly_mat.c, last updated 2025-05-26, F.Osorio */

#include "fastmatrix.h"

void
matrix_polynomial(double *a, int *lda, int *n, double *coef, int *ncoef, double *b, int *ldb, int *info)
{ /* evaluates a real general matrix polynomial using a Horner's scheme */
  FM_polymat(a, *lda, *n, coef, *ncoef, b, *ldb, info);
}

void
FM_polymat(double *a, int lda, int n, double *coef, int ncoef, double *b, int ldb, int *info)
{ /* evaluates a real general matrix polynomial */
  int k = ncoef;
  double *s;

  /* test the input parameters */
  *info = 0;
  if (n < 0) {
    *info = -3;
  } else if (lda < MAX(1, n)) {
    *info = -2;
  } else if (ncoef < 0) {
    *info = -5;
  } else if (ldb < MAX(1, n)) {
    *info = -7;
  }
  if (*info != 0) return;

  /* quick return if possible */
  if (n == 0)
    return;

  /* quick return */
  if (ncoef == 1) {
    for (int j = 0; j < n; j++)
      b[j * (ldb + 1)] = coef[0];
    return;
  }

  /* initialization */
  s = (double *) R_Calloc(n * n, double);
  k--;
  FM_scale_mat(s, n, coef[k], a, lda, n, n);
  k--;
  for (int j = 0; j < n; j++)
    s[j * (n + 1)] += coef[k];

  /* Horner's scheme */
  while (k-- > 0) {
    FM_mult_mat(s, a, lda, n, n, s, n, n, n);
    for (int j = 0; j < n; j++)
      s[j * (n + 1)] += coef[k];
  }

  /* saving result */
  FM_copy_mat(b, ldb, s, n, n, n);

  R_Free(s);
}
