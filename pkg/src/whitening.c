/* $ID: whitening.c, last updated 2021-04-03, F.Osorio */

#include "fastmatrix.h"

void
whitening_chol(double *y, int *nrow, int *ncol, double *Scatter)
{ /* Cholesky whitening transformation  */
  char *side = "R", *uplo = "L", *trans = "T", *diag = "N";
  int info = 0, job = 0, n = *nrow, p = *ncol;

  FM_chol_decomp(Scatter, p, p, job, &info);
  if (info)
    error("DPOTRF in cholesky decomposition gave code %d", info);

  FM_invert_triangular(Scatter, p, p, job, &info);
  if (info)
    error("DTRTRI in matrix inversion gave code %d", info);

  BLAS3_trmm(1.0, Scatter, p, n, p, side, uplo, trans, diag, y, n);
}
