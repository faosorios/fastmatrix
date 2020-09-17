/* $ID: duplication.c, last updated 2020-08-10, F.Osorio */

#include "fastmatrix.h"

void
symmetrizer_prod(double *a, int *lda, int *arow, int *acol, double *b, int *ldb)
{ /* computes: B <- N %*% A or B <- A %*% N */
  add_mat(a, *lda, 0.5, b, *ldb, *arow, *acol);
}
