/* $ID: polymat, last updated 2022-07-04, F.Osorio */

#include "fastmatrix.h"

void
matrix_polynomial(double *a, int *lda, int *n, double *coef, int *ncoef, double *b, int *ldb, int *info)
{ /* evaluates a real general matrix polynomial using a Horner's scheme */
  FM_matrix_pol(a, *lda, *n, coef, *ncoef, b, *ldb, info);
}

void
FM_matrix_pol(double *a, int lda, int n, double *coef, int ncoef, double *b, int ldb, int *info)
{ /* evaluates a real general matrix polynomial */

  /* test the input parameters */
  *info = 0;
  if (n < 0) {
    *info = 3;
  } else if (lda < MAX(1, n)) {
    *info = 2;
  } else if (ncoef < 0) {
    *info = 5;
  } else if (ldb < MAX(1, n)) {
    *info = 7;
  }
  if (*info != 0) return;

  /* quick return if possible */
  if (n == 0)
    return;

  if (ncoef == 0) {
    for (int j = 0; j < n; j++)
      b[j * (n + 1)] = coef[0];
    return;
  }

  /* Horner's scheme */
  for (int j = 0; j < n; j++)
    b[j * (n + 1)] = coef[ncoef];
  while (ncoef-- > 0) {
    FM_mult_mat(b, a, lda, n, n, b, ldb, n, n);
    for (int j = 0; j < n; j++)
      b[j * (n + 1)] += coef[ncoef];
  }
}
