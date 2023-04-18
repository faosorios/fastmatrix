/* $ID: sweep_operator.c, last updated 2020-08-21, F.Osorio */

#include "fastmatrix.h"

void
sweep_operator(double *a, int *lda, int *p, int *which, int *r, int *job)
{ /* wrapper to Fortran 'sweepop' */
  int info = 0, k;

  for (int j = 0; j < *r; j++) {
    k = which[j];
    F77_CALL(sweepop)(a, lda, p, &k, job, &info);
    if (info)
      error("symmetric sweep operator gave code %d", info);
  }
}
