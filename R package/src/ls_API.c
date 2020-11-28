/* $ID: ls_API.c, last updated 2020-11-20, F.Osorio */

#include "fastmatrix.h"

void
FM_lsfit(double *x, int ldx, int nrow, int ncol, double *y, int ldy, int nrhs, double *coef, int *info)
{ /* solve (overdeterminated) least squares problems, multivariate responses are allowed */
  char *notrans = "N";
  int errcode, lwork;
  double opt, *work;

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dgels)(notrans, &nrow, &ncol, &nrhs, x, &ldx, y, &ldy, &opt, &lwork, &errcode);
  if (errcode != 0)
    error("DGELS in ordinary least squares gave error code %d", errcode);

  /* calling DGELS with optimal size of working array */
  lwork = (int) opt;
  work = (double *) Calloc(lwork, double);
  F77_CALL(dgels)(notrans, &nrow, &ncol, &nrhs, x, &ldx, y, &ldy, work, &lwork, info);
  FM_copy_mat(coef, ncol, y, ldy, ncol, nrhs);
  Free(work);
}

void
FM_gls_GQR(double *x, int ldx, int nrow, int ncol, double *y, double *cov, double *coef, int *info)
{ /* solve generalized least squares problems using a generalized QR factorization */
  int errcode, lwork, job = 0;
  double opt, dummy, *work, *b, *u;

  /* cholesky factorization of 'cov' */
  b = (double *) Calloc(nrow * nrow, double);
  FM_cpy_lower(cov, nrow, nrow, b, nrow);
  FM_chol_decomp(b, nrow, nrow, job, &errcode);
  if (errcode != 0)
    error("cholesky decomposition in generalized least squares gave error code %d", errcode);

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dggglm)(&nrow, &ncol, &nrow, x, &ldx, b, &nrow, y, coef, &dummy, &opt, &lwork, &errcode);

  if (errcode != 0)
    error("DGGGLM in generalized least squares gave error code %d", errcode);

  /* calling DGGGLM with optimal size of working array */
  lwork = (int) opt;
  work = (double *) Calloc(lwork, double);
  u = (double *) Calloc(nrow, double);
  F77_CALL(dggglm)(&nrow, &ncol, &nrow, x, &ldx, b, &nrow, y, coef, u, work, &lwork, &errcode);
  Free(b); Free(u); Free(work);
}
