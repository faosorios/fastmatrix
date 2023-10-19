/* $ID: lu.c, last updated 2023-09-14, F.Osorio */

#include "fastmatrix.h"

void
lu_dcmp(double *a, int *lda, int *n, int *p, int *pivot, int *info)
{ /* wrapper to 'FM_lu_decomp' */
  FM_lu_decomp(a, *lda, *n, *p, pivot, info);
}

void
lu_inverse(double *a, int *lda, int *p, int *pivot, int *info)
{ /* wrapper to 'FM_lu_inverse' */
  FM_lu_inverse(a, *lda, *p, pivot, info);
}

void
lu_solve(double *a, int *lda, int *p, int *pivot, double *b, int *ldb, int *nrhs, int *info)
{ /* wrapper to 'FM_lu_solve' */
  FM_lu_solve(a, *lda, *p, pivot, b, *ldb, *nrhs, info);
}
