/* ID: cholesky.c, last updated 2022-02-10, F.Osorio */

#include "fastmatrix.h"

/* Cholesky decompositions */

void
chol_dcmp(double *a, int *lda, int *p, int *job, int *info)
{ /* wrapper to 'FM_chol_decomp' */
  FM_chol_decomp(a, *lda, *p, *job, info);
}

void
chol_update(double *r, int *ldr, int *p, double *x)
{ /* update the Cholesky decomposition */
  int n = *p;
  double *c, *s;

  /* cosines and sines of transforming rotations */
  c = (double *) Calloc(n, double);
  s = (double *) Calloc(n, double);

  /* update upper triangular matrix */
  for (int j = 0; j < n; j++) {
    double xj = x[j];

    /* apply the previous rotations */
    for (int i = 0; i < j; i++) {
      double aux = *(r + i + j * n) * c[i] + xj * s[i];
      xj = xj * c[i] - *(r + i + j * n) * s[i];
      *(r + i + j * n) = aux;
    }

    /* compute next rotation */
    BLAS1_rotg((r + j + j * n), &xj, (c + j), (s + j));
  }

  Free(c); Free(s);
}
