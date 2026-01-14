/* $ID: rank1_update.c, last updated 2026-01-13, F.Osorio */

#include "fastmatrix.h"

void
rank1_update(double *mat, int *ldmat, int *n, int *job, double *a, double *alpha, double *u)
{ /* Rank-1 update (to be called by R) */
  int p = *n, task = *job;

  if (task) {
    FM_copy_mat(mat, *ldmat, a, p, p, p);
    if (*alpha == 0.0)
      return;
    BLAS2_ger(*alpha, mat, *ldmat, p, p, u, 1, u, 1);
  }
  else 
    FM_diag_plus_rank1(mat, *ldmat, p, a, *alpha, u);
}

void 
FM_diag_plus_rank1(double *mat, int ldmat, int n, double *a, double alpha, double *u)
{ /* rank 1 update: mat <- diag(a) + alpha * u %*% t(u) (symbol exported to API) */

  if (n <= 0) /*quick return */
    return;

  if (alpha == 0.0) {
    for (int i = 0; i < n; i++)
      mat[i * (ldmat + 1)] = a[i];
    return;
  }

  for (int i = 0; i < n; i++) {
    mat[i * (ldmat + 1)] = a[i] + alpha * SQR(u[i]);
    for (int j = i + 1; j < n; j++) {
      *(mat + i + j * n) = alpha * u[i] * u[j];
      *(mat + j + i * n) = *(mat + i + j * n);
    }
  }
}
