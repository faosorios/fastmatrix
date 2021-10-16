/* $ID: lu.c, last updated 10-14-2021, F.Osorio */

#include "fastmatrix.h"

void
lu_dcmp(double *a, int *lda, int *n, int *p, int *pivot)
{ /* LU factorization of a real square matrix,
   * matrix 'a' is overwritten with the result */
  int info = 0;

  F77_CALL(dgetrf)(n, p, a, lda, pivot, &info);
  if (info)
    error("lu_dcmp gave code %d", info);
}

void
lu_inverse(double *a, int *lda, int *p, int *pivot)
{ /* computes the inverse of a matrix using the LU factorization */
  int info = 0, lwork = *p;
  double *work;

  work = (double *) Calloc(lwork, double);
  F77_CALL(dgetri)(p, a, lda, pivot, work, &lwork, &info);
  Free(work);
  if (info)
    error("lu_inverse gave code %d", info);
}

void
lu_solve(double *a, int *lda, int *p, int *pivot, double *b, int *ldb, int *nrhs)
{ /* solves a system of linear equations */
  char *notrans = "N";
  int info = 0;

  F77_CALL(dgetrs)(notrans, p, nrhs, a, lda, pivot, b, ldb, &info FCONE);
  if (info)
    error("lu_solve gave code %d", info);
}
