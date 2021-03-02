/* ID: cholesky.c, last updated 2021-02-17, F.Osorio */

#include "fastmatrix.h"

/* Cholesky decompositions */

void
chol_dcmp(double *a, int *lda, int *p, int *job, int *info)
{ /* wrapper to 'FM_chol_decomp' */
  FM_chol_decomp(a, *lda, *p, *job, info);
}
