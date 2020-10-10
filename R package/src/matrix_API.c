/* ID: matrix_API.c, last updated 2020-10-08, F.Osorio */

#include "fastmatrix.h"

/* basic matrix manipulations */

void
FM_add_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol)
{ /* y <- y + alpha * x */
  for (int j = 0; j < ncol; j++) {
    BLAS1_axpy(alpha, x, 1, y, 1, nrow);
    y += ldy; x += ldx;
  }
}

void
FM_copy_mat(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
{ /* y <- x[,] */
  for (int j = 0; j < ncol; j++) {
    Memcpy(y, x, nrow);
    y += ldy; x += ldx;
  }
}

void
FM_copy_trans(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
{ /* y <- t(x) */
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++)
      *(y + j + i * ldy) = *(x + i + j * ldx);
  }
}

void
FM_scale_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol)
{ /* y <- alpha * x[,] */
  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++)
      y[i] = alpha * x[i];
    y += ldy; x += ldx;
  }
}

void
FM_setzero(double *y, int ldy, int nrow, int ncol)
{ /* y[,] <- 0, sets all elements of y to 0 */
  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++)
      y[i] = 0.0;
    y += ldy;
  }
}

void
FM_GAXPY(double *y, double alpha, double *a, int lda, int nrow, int ncol, double *x, double beta, int job)
{ /* y <- alpha * a %*% x    + beta * y (job = 0), or
   * y <- alpha * t(a) %*% x + beta * y (job = 1) */
  char *trans;
  int inc = 1;

  trans = (job) ? "T" : "N";
  F77_CALL(dgemv)(trans, &nrow, &ncol, &alpha, a, &lda, x, &inc, &beta, y, &inc);
}

void
FM_mult_triangular(double *y, double *a, int lda, int n, double *x, int job)
{ /* y <- lower.tri(a) %*% x (job = 0), or
   * y <- upper.tri(a) %*% x (job = 1), where 'x' is a vector */
  char *uplo, *trans = "N", *diag = "N";
  int inc = 1;

  uplo = (job) ? "U" : "L";
  Memcpy(y, x, n);
  F77_CALL(dtrmv)(uplo, trans, diag, &n, a, &lda, y, &inc);
}

void
FM_mult_mat(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols)
{ /* matrix multiplication of two conformable matrices. z <- x %*% y */
  char *notransx = "N", *notransy = "N";
  double one = 1.0, zero = 0.0, *tmp = NULL;

  /* use tmp so z can be either x or y */
  tmp = (double *) Calloc(xrows * ycols, double);
  F77_CALL(dgemm)(notransx, notransy, &xrows, &ycols, &xcols, &one, x, &ldx, y, &ldy, &zero, tmp, &xrows);
  Memcpy(z, tmp, xrows * ycols);
  Free(tmp);
}

void
FM_crossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols)
{ /* cross product of two given matrices. z <- t(x) %*% y */
  char *transx = "T", *notransy = "N";
  double one = 1.0, zero = 0.0, *tmp = NULL;

  /* use tmp so z can be either x or y */
  tmp = (double *) Calloc(xcols * ycols, double);
  F77_CALL(dgemm)(transx, notransy, &xcols, &ycols, &xrows, &one, x, &ldx, y, &ldy, &zero, tmp, &xcols);
  Memcpy(z, tmp, xcols * ycols);
  Free(tmp);
}

void
FM_tcrossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols)
{ /* outer product of two given matrices. z <- x %*% t(y) */
  char *notransx = "N", *transy = "T";
  double one = 1.0, zero = 0.0, *tmp = NULL;

  /* use tmp so z can be either x or y */
  tmp = (double *) Calloc(xrows * yrows, double);
  F77_CALL(dgemm)(notransx, transy, &xrows, &yrows, &xcols, &one, x, &ldx, y, &ldy, &zero, tmp, &xrows);
  Memcpy(z, tmp, xrows * yrows);
  Free(tmp);
}

void
FM_rank1_update(double *a, int lda, int nrow, int ncol, double alpha, double *x, double *y)
{ /* rank 1 update: a <- alpha * x %*% t(y) + a */
  BLAS2_ger(alpha, a, lda, nrow, ncol, x, 1, y, 1);
}

double
FM_logAbsDet(double *a, int lda, int n)
{ /* log(abs(det(a))), where 'a' is a triangular matrix */
  double accum = 0.0;

  for (int i = 0; i < n; i++)
    accum += log(fabs(a[i * (lda + 1)]));
  return accum;
}

double
FM_trace(double *a, int lda, int n)
{ /* trace of a square matrix */
  double accum = 0.0;

  for (int i = 0; i < n; i++)
    accum += a[i * (lda + 1)];
  return accum;
}

/* routines for matrix decompositions */

void
FM_chol_decomp(double *a, int lda, int p, int job, int *info)
{ /* cholesky factorization of a real symmetric positive definite matrix a.
   * the factorization has the form:
   * a <- l %*% t(l), if job = 0, or
   * a <- t(u) %*% u, if job = 1,
   * where u is an upper triangular matrix and l is lower triangular */
  char *uplo;

  uplo = (job) ? "U" : "L";
  F77_CALL(dpotrf)(uplo, &p, a, &lda, info);
}

void
FM_svd_decomp(double *mat, int ldmat, int nrow, int ncol, double *u, int ldu, double *d, double *v, int ldv, int job, int *info)
{ /* return the SVD decomposition of a rectangular matrix */
  double *work, *upper;

  work  = (double *) Calloc(nrow, double);
  upper = (double *) Calloc(ncol, double);
  F77_CALL(dsvdc)(mat, &ldmat, &nrow, &ncol, d, upper, u, &ldu, v, &ldv, work, &job, info);
  Free(work); Free(upper);
}

void
FM_QR_decomp(double *mat, int ldmat, int nrow, int ncol, double *qraux, int *info)
{ /* return the QR decomposition of a rectangular matrix */
  int errcode = 0, lwork;
  double opt, *work;

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dgeqrf)(&nrow, &ncol, mat, &ldmat, qraux, &opt, &lwork, &errcode);
  if (errcode != 0)
    error("DGEQRF in QR decomposition gave error code %d", errcode);

  /* calling DGEQRF with optimal size of working array */
  lwork = (int) opt;
  work = (double *) Calloc(lwork, double);
  F77_CALL(dgeqrf)(&nrow, &ncol, mat, &ldmat, qraux, work, &lwork, info);
  Free(work);
}

void
FM_QR_qy(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info)
{ /* ymat <- qr.qy(qr, ymat) */
  char *side = "L", *notrans = "N";
  int errcode = 0, lwork, nrflc, nrhs;
  double opt, *work;

  nrhs  = ycol;
  nrflc = ncol; /* length of 'qraux' */

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dormqr)(side, notrans, &yrow, &nrhs, &nrflc, qr, &ldq, qraux, ymat, &ldy, &opt, &lwork, &errcode);
  if (errcode != 0)
    error("DORMQR in QR_qy gave error code %d", info);

  /* calling DORMQR with optimal size of working array */
  lwork = (int) opt;
  work = (double *) Calloc(lwork, double);
  F77_CALL(dormqr)(side, notrans, &yrow, &nrhs, &nrflc, qr, &ldq, qraux, ymat, &ldy, work, &lwork, info);
  Free(work);
}

void
FM_QR_qty(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info)
{ /* ymat <- qr.qty(qr, ymat) */
  char *side = "L", *trans = "T";
  int errcode = 0, lwork, nrflc, nrhs;
  double opt, *work;

  nrhs  = ycol;
  nrflc = ncol; /* length of 'qraux' */

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dormqr)(side, trans, &yrow, &nrhs, &nrflc, qr, &ldq, qraux, ymat, &ldy, &opt, &lwork, &errcode);
  if (errcode != 0)
    error("DORMQR in QR_qty gave error code %d", info);

  /* calling DORMQR with optimal size of working array */
  lwork = (int) opt;
  work = (double *) Calloc(lwork, double);
  F77_CALL(dormqr)(side, trans, &yrow, &nrhs, &nrflc, qr, &ldq, qraux, ymat, &ldy, work, &lwork, info);
  Free(work);
}

void
FM_QR_fitted(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int job, double *fitted)
{ /* fitted <- qr.fitted(qr, ymat) */
  int errcode = 0;

  FM_setzero(fitted, ldy, yrow, ycol); /* is not assumed that 'fitted' is zeroed */
  if (job) { /* ymat <- qr.qty(qr, ymat)*/
    FM_QR_qty(qr, ldq, nrow, ncol, qraux, ymat, ldy, yrow, ycol, &errcode);
    if (errcode != 0)
      error("DORMQR in QR_fitted gave error code %d", errcode);
  }
  FM_copy_mat(fitted, ldy, ymat, ldy, ncol, ycol); /* copying 1st 'ncol' rows */
  FM_QR_qy(qr, ldq, nrow, ncol, qraux, fitted, ldy, yrow, ycol, &errcode);
  if (errcode != 0)
    error("DORMQR in QR_fitted gave error code %d", errcode);
}

void
FM_QR_store_R(double *qr, int ldq, int nrow, int ncol, double *Dest, int ldDest)
{ /* copy the R part into Dest */
  int rows;

  for (int j = 0; j < ncol; j++) {
    rows = MIN(j + 1, nrow);
    Memcpy(Dest + j * ldDest, qr + j * ldq, rows);
  }
}

/* matrix inversion and linear solvers */

void
FM_invert_mat(double *a, int lda, int n, int *info)
{ /* performs matrix inversion */
  char *notrans = "N";
  int errcode = 0, lwork;
  double opt, *b, *dummy = NULL, *work;

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dgels)(notrans, &n, &n, &n, a, &lda, dummy, &n, &opt, &lwork, &errcode);
  if (errcode != 0)
    error("DGELS in invert_mat gave error code %d", errcode);

  /* calling DGELS with optimal size of working array */
  lwork = (int) opt;
  work  = (double *) Calloc(lwork, double);
  b     = (double *) Calloc(n * n, double);
  for (int j = 0; j < n; j++)
    b[j * (n + 1)] = 1.0;
  F77_CALL(dgels)(notrans, &n, &n, &n, a, &lda, b, &n, work, &lwork, info);
  Memcpy(a, b, n * n);
  Free(b); Free(work);
}

void
FM_invert_triangular(double *a, int lda, int n, int job, int *info)
{ /* computes the inverse of an upper (job = 1) or lower (job = 0)
   * triangular matrix in place */
  char *diag = "N", *uplo;

  uplo = (job) ? "U" : "L";
  F77_CALL(dtrtri)(uplo, diag, &n, a, &lda, info);
}

void
FM_chol_inverse(double *a, int lda, int p, int job, int *info)
{ /* inverse matrix from cholesky (or QR) decomposition
   * the lower triangular factor (job = 0), or the upper triangular
   * factor (job = 1) is used */
  char *uplo;

  uplo = (job) ? "U" : "L";
  F77_CALL(dpotri)(uplo, &p, a, &lda, info);
}

void
FM_backsolve(double *r, int ldr, int n, double *b, int ldb, int nrhs, int *info)
{ /* backsolve solve triangular systems of the form r %*% x = b, where r
   * is an upper triangular matrix and b is a matrix containing the right-hand
   * sides to equations */
  char *diag = "N", *uplo = "U", *notrans = "N";

  F77_CALL(dtrtrs)(uplo, notrans, diag, &n, &nrhs, r, &ldr, b, &ldb, info);
}

void
FM_forwardsolve(double *l, int ldl, int n, double *b, int ldb, int nrhs, int *info)
{ /* forwardsolve solve triangular systems of the form l %*% x = b, where l
   * is a lower triangular matrix and b is a matrix containing the right-hand
   * sides to equations */
  char *diag = "N", *uplo = "L", *notrans = "N";

  F77_CALL(dtrtrs)(uplo, notrans, diag, &n, &nrhs, l, &ldl, b, &ldb, info);
}

/* DEBUG routine */

void
FM_print_mat(double *x, int ldx, int nrow, int ncol, char *msg)
{ /* print matrix and message (used for printf debugging) */
  Rprintf( "%s\n", msg);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++)
      Rprintf( " %10.5g", x[i + j * ldx ]);
    Rprintf( "\n" );
  }
  Rprintf( "\n" );
}
