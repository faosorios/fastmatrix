/* $ID: matrix.c, last updated 2020-08-12, F.Osorio */

#include "base.h"
#include "norms.h"

void
norm_one(double *x, int *inc, int *n, double *value)
{ /* absolute-value norm */
  *value = F77_CALL(dasum)(n, x, inc);
}

void
norm_two(double *x, int *inc, int *n, double *value)
{ /* Euclidean norm */
  *value = F77_CALL(dnrm2)(n, x, inc);
}

void
norm_inf(double *x, int *inc, int *n, double *value)
{ /* infinity norm */
  *value = F77_CALL(dnrminf)(n, x, inc);
}

void
norm_minkowski(double *x, int *inc, int *n, double *p, double *value)
{ /* Minkowski norm */
  *value = F77_CALL(minkowski)(n, x, inc, p);
}

void
matrix_norm(double *a, int *lda, int *nrow, int *ncol, int *job, double *value)
{ /* wrapper to LAPACK DLANGE */
  char *task;
  double *work = NULL;

  switch (*job) {
    case 0:
      task = "I";
      work = (double *) Calloc(*nrow, double);
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

  *value = F77_CALL(dlange)(task, nrow, ncol, a, lda, work);
  if (*job == 0)
    Free(work);
}
