/* $ID: krylov, last updated 2024-09-03, F.Osorio */

#include "fastmatrix.h"

void
krylov_mat(double *a, int *lda, int *n, double *b, int *m, double *k, int *ldk, int *info)
{ /* construcs the Krylov matrix based on an n-by-n matrix a and an n-vector b */
  FM_krylov_mat(a, *lda, *n, b, *m, k, *ldk, info);
}

void
FM_krylov_mat(double *a, int lda, int n, double *b, int m, double *k, int ldk, int *info)
{ /* given an n-by-n matrix a and an n-vector b, this function constructs the Krylov
   * matrix k, where k = [b, Ab, ..., A^(m-1)b] */
  double *work;

  /* test the input parameters */
  *info = 0;
  if (n < 0) {
    *info = -3;
  } else if (lda < MAX(1, n)) {
    *info = -2;
  } else if (m < 0) {
    *info = -5;
  } else if (ldk < MAX(1, n)) {
    *info = -7;
  }
  if (*info != 0) return;

  /* quick return if possible */
  if ((n == 0) || (m == 0))
    return;

  work = (double *) R_Calloc(n, double);
  Memcpy(work, b, n);

  /* 1st column */
  Memcpy(k, work, n);
  k += ldk;

  /* next columns */
  for (int j = 1; j < m; j++) {
    FM_mult_mat(work, a, lda, n, n, work, n, n, 1);
    Memcpy(k, work, n);
    k += ldk;
  }

  R_Free(work);
}
