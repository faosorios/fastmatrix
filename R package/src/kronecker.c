/* $ID: kronecker.c, last updated 2020-09-21, F.Osorio */

#include "fastmatrix.h"

void
kronecker_prod(double *a, int *arow, int *acol, double *b, int *brow, int *bcol, double *prod)
{ /* kronecker product between two matrices */
  int nrow = *arow * *brow, irow, jcol;
  double z;

  for (int i = 0; i < *arow; i++) {
    irow = i * *brow;
    for (int j = 0; j < *acol; j++) {
      jcol = j * *bcol;
      z = a[j * *arow + i];
      for (int k = 0; k < *brow; k++) {
        for (int l = 0; l < *bcol; l++)
          prod[(jcol + l) * nrow + irow + k] = z * b[l * *brow + k];
      }
    }
  }
}
