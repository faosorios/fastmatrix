/* $ID: norms.c, last updated 2024-09-03, F.Osorio */

#include "fastmatrix.h"

void
norm_one(double *x, int *inc, int *n, double *value)
{ /* absolute-value norm */
  *value = BLAS1_sum_abs(x, *inc, *n);
}

void
norm_two(double *x, int *inc, int *n, double *value)
{ /* Euclidean norm */
  *value = BLAS1_norm_two(x, *inc, *n);
}

void
norm_inf(double *x, int *inc, int *n, double *value)
{ /* infinity norm */
  int idx = 0;

  idx = BLAS1_index_max(x, *inc, *n);
  idx--; /* index correction */
  *value = fabs(x[idx]);
}

void
norm_minkowski(double *x, int *inc, int *n, double *p, double *value)
{ /* Minkowski p-norm */
  *value = F77_CALL(minkowski)(n, x, inc, p);
}

void
matrix_norm(double *a, int *lda, int *nrow, int *ncol, int *job, double *value)
{ /* wrapper to LAPACK DLANGE */
  char *task = NULL;
  double *work = NULL;

  switch (*job) {
    case 0:
      task = "I";
      work = (double *) R_Calloc(*nrow, double);
      break;
    case 1:
      task = "1";
      break;
    case 2:
      task = "F";
      break;
    case 3:
      task = "M";
      break;
  }

  *value = F77_CALL(dlange)(task, nrow, ncol, a, lda, work FCONE);
  if (*job == 0)
    R_Free(work);
}
