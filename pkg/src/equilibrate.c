/* $ID: equilibrate, last updated 2021-11-10, F.Osorio */

#include "fastmatrix.h"

void
equilibrate_pd(double *a, int *lda, int *p, double *scales, double *cond, double *largest, int *info)
{ /* equilibrate a symmetric positive definite matrix 'a' and reduce its
   * condition number (with respect to the two-norm) */
  F77_CALL(dpoequ)(p, a, lda, scales, cond, largest, info);
}

void
equilibrate_sym(double *a, int *lda, int *p, double *scales, double *cond, double *largest, int *info)
{ /* computes row and column scalings intended to equilibrate a symmetric
   * matrix 'a' (with respect to the Euclidean norm) and reduce its condition
   * number */
  char *task = "U";
  double *work;

  work = (double *) Calloc(2 * *p, double);
  F77_CALL(dsyequb)(task, p, a, lda, scales, cond, largest, work, info FCONE);
  Free(work);
}

void
equilibrating_sym(double *a, int *n, double *scales)
{ /* scales a symmetric matrix in order to reduce its condition number */
  int p = *n;

  for (int i = 0; i < p; i++) {
    a[i * (p + 1)] *= SQR(scales[i]);
    for (int j = i + 1; j < p; j++) {
      *(a + i + j * p) = *(a + i + j * p) * scales[i] * scales[j];
      *(a + j + i * p) = *(a + i + j * p);
    }
  }
}
