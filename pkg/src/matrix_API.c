/* ID: matrix_API.c, last updated 2026-01-12, F.Osorio */

#include "fastmatrix.h"

/* operations on vectors */

double
FM_norm_sqr(double *x, int inc, int n)
{ /* sum(x * x) */
  double length;

  length = BLAS1_norm_two(x, inc, n);
  return R_pow_di(length, 2);
}

void
FM_normalize(double *x, int inc, int n)
{ /* x <- x / sqrt(sum(x * x)) */
  double div = 1.0, length;

  length = BLAS1_norm_two(x, inc, n);
  div /= length;
  BLAS1_scale(div, x, inc, n);
}

double
FM_vecsum(double *x, int inc, int n)
{ /* sum(x) with increments 'inc' */
  int ix = 0;
  double accum = 0.0;

  for (int i = 0; i < n; i++) {
    accum += x[ix];
    ix += inc;
  }
  return accum;
}

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
  F77_CALL(dgemv)(trans, &nrow, &ncol, &alpha, a, &lda, x, &inc, &beta, y, &inc FCONE);
}

void
FM_mult_triangular(double *y, double *a, int lda, int n, double *x, int job)
{ /* y <- lower.tri(a) %*% x (job = 0), or
   * y <- upper.tri(a) %*% x (job = 1), where 'x' is a vector */
  char *uplo, *trans = "N", *diag = "N";
  int inc = 1;

  uplo = (job) ? "U" : "L";
  Memcpy(y, x, n);
  F77_CALL(dtrmv)(uplo, trans, diag, &n, a, &lda, y, &inc FCONE FCONE FCONE);
}

void
FM_mult_mat_vec(double *y, double *a, int lda, int n, int p, double *x)
{ /* performs matrix-vector multiplication, y <- a %*% x */
  for (int j = 0; j < p; j++)
    BLAS1_axpy(x[j], a + j * lda, 1, y, 1, n);
}

void
FM_mult_mat(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols)
{ /* matrix multiplication of two conformable matrices. z <- x %*% y */
  char *notrans = "N";
  double one = 1.0, zero = 0.0, *tmp = NULL;

  /* use tmp so z can be either x or y */
  tmp = (double *) R_Calloc(xrows * ycols, double);
  F77_CALL(dgemm)(notrans, notrans, &xrows, &ycols, &xcols, &one, x, &ldx, y, &ldy, &zero, tmp, &xrows FCONE FCONE);
  Memcpy(z, tmp, xrows * ycols);
  R_Free(tmp);
}

void
FM_crossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols)
{ /* cross product of two given matrices. z <- t(x) %*% y */
  char *trans = "T", *notrans = "N";
  double one = 1.0, zero = 0.0, *tmp = NULL;

  /* use tmp so z can be either x or y */
  tmp = (double *) R_Calloc(xcols * ycols, double);
  F77_CALL(dgemm)(trans, notrans, &xcols, &ycols, &xrows, &one, x, &ldx, y, &ldy, &zero, tmp, &xcols FCONE FCONE);
  Memcpy(z, tmp, xcols * ycols);
  R_Free(tmp);
}

void
FM_tcrossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols)
{ /* outer product of two given matrices. z <- x %*% t(y) */
  char *notrans = "N", *trans = "T";
  double one = 1.0, zero = 0.0, *tmp = NULL;

  /* use tmp so z can be either x or y */
  tmp = (double *) R_Calloc(xrows * yrows, double);
  F77_CALL(dgemm)(notrans, trans, &xrows, &yrows, &xcols, &one, x, &ldx, y, &ldy, &zero, tmp, &xrows FCONE FCONE);
  Memcpy(z, tmp, xrows * yrows);
  R_Free(tmp);
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

/* operations on triangular matrices */

void
FM_cpy_upper(double *x, int ldx, int p, double *y, int ldy)
{ /* upper.tri(y) <- upper.tri(x) */
  for (int j = 0; j < p; j++) {
    for (int i = j; i < p; i++)
      y[j + i * ldy] = x[j + i * ldx];
  }
}

void
FM_cpy_lower(double *x, int ldx, int p, double *y, int ldy)
{ /* lower.tri(y) <- lower.tri(x) */
  for (int j = 0; j < p; j++) {
    for (int i = j; i < p; i++)
      y[i + j * ldy] = x[i + j * ldx];
  }
}

void
FM_cpy_upper2lower(double *x, int ldx, int p, double *y)
{ /* lower.tri(y) <- upper.tri(x), only the strictly upper triangular
   * part is copied to the strictlty lower triangular part of 'y',
   * 'x' and 'y' must be square matrices, 'y' can be 'x' */
  for (int j = 0; j < p; j++) {
    for (int i = j; i < p; i++)
      y[i + j * ldx] = x[j + i * ldx];
  }
}

void
FM_cpy_lower2upper(double *x, int ldx, int p, double *y)
{ /* upper.tri(y) <- lower.tri(x), only the strictly lower triangular
   * part is copied to the strictlty upper triangular part of 'y',
   * 'x' and 'y' must be square matrices, 'y' can be 'x' */
  for (int j = 0; j < p; j++) {
    for (int i = j; i < p; i++)
      y[j + i * ldx] = x[i + j * ldx];
  }
}

double
FM_sum_upper_tri(double *x, int ldx, int p, int job)
{ /* sum(upper.tri(x)) */
  double accum = 0.0;

  if (job) {
    for (int j = 0; j < p; j++) {
      for (int i = j; i < p; i++)
        accum += x[j + i * ldx];
    }
  } else { /* strictly upper part */
    for (int j = 0; j < p; j++) {
      for (int i = j + 1; i < p; i++)
        accum += x[j + i * ldx];
    }
  }

  return accum;
}

double
FM_sum_lower_tri(double *x, int ldx, int p, int job)
{ /* sum(lower.tri(x)) */
  double accum = 0.0;

  if (job) {
    for (int j = 0; j < p; j++) {
      for (int i = j; i < p; i++)
        accum += x[i + j * ldx];
    }
  } else { /* strictly lower part */
    for (int j = 0; j < p; j++) {
      for (int i = j + 1; i < p; i++)
        accum += x[i + j * ldx];
    }
  }

  return accum;
}

/* other matrix operations */

double 
FM_bilinear_form(double *a, int lda, int n, int p, double *x, double *y)
{ /* this function computes the bilinear form, t(x) %*% A %*% y */
  char *notrans = "N";
  double value = 0.0, *z;

  /* quick return if possible */
  if (n <= 0) 
    return 0.0;
  if (p <= 0)
    return 0.0;
  if (lda < (n > 1 ? n : 1)) 
    return 0.0;

  /* start operations */
  z = (double *) R_Calloc(n, double);
  BLAS2_gemv(1.0, a, lda, n, p, notrans, y, 1, 0.0, z, 1);
  value = BLAS1_dot_product(x, 1, z, 1, n);
  R_Free(z);

  return value;
}

double 
FM_quadratic_form(double *a, int lda, int n, double *x)
{ /* this function computes the quadratic form, t(x) %*% A %*% x */
  char *notrans = "N";
  double value = 0.0, *z;

  /* quick return if possible */
  if (n <= 0) 
    return 0.0;
  if (lda < (n > 1 ? n : 1)) 
    return 0.0;

  /* start operations */
  z = (double *) R_Calloc(n, double);
  BLAS2_gemv(1.0, a, lda, n, n, notrans, x, 1, 0.0, z, 1);
  value = BLAS1_dot_product(x, 1, z, 1, n);
  R_Free(z);

  return value;
}

void 
FM_murrv(double *y, double *a, int lda, int n, int p, double *x, int *info)
{ /* multiplies a real rectangular matrix by a vector, y = a %*% x */
  char *notrans = "N";

  *info = 0;
  /* quick return if possible */
  if ((n == 0) || (p == 0)) 
    return;

  /* test the input parameters */
  if (n < 0) {
    *info = -4;
  } else if (p < 0) {
    *info = -5;
  } else if (lda < (n > 1 ? n : 1)) {
    *info = -3;
  }
  if (*info != 0)
    return;

  /* y <- 1.0 * a %*% x + 0.0 * y */
  BLAS2_gemv(1.0, a, lda, n, p, notrans, x, 1, 0.0, y, 1);
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
  F77_CALL(dpotrf)(uplo, &p, a, &lda, info FCONE);
}

void
FM_lu_decomp(double *a, int lda, int n, int p, int *pivot, int *info)
{ /* LU factorization of a real square matrix,
   * matrix 'a' is overwritten with the result */

  F77_CALL(dgetrf)(&n, &p, a, &lda, pivot, info);
}

void
FM_svd_decomp(double *mat, int ldmat, int nrow, int ncol, double *u, int ldu, double *d, double *v, int ldv, int job, int *info)
{ /* return the SVD decomposition of a rectangular matrix */
  double *work, *upper;

  work  = (double *) R_Calloc(nrow, double);
  upper = (double *) R_Calloc(ncol, double);
  F77_CALL(dsvdc)(mat, &ldmat, &nrow, &ncol, d, upper, u, &ldu, v, &ldv, work, &job, info);
  R_Free(work); R_Free(upper);
}

void
FM_schur_decomp(double *a, int lda, int n, int task, double *re, double *im, double *v, int ldv, int *info)
{ /* return the Schur decomposition of a nonsymmetric matrix */
  int errcode = 0, bwork = 0, lwork;
  double opt, *work;

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(schur_wrapper)(a, &lda, &n, &task, re, im, v, &ldv, &opt, &lwork, &bwork, &errcode);
  if (errcode != 0)
    error("DGEES in Schur decomposition gave error code %d", errcode);

  /* calling DGEES with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(schur_wrapper)(a, &lda, &n, &task, re, im, v, &ldv, work, &lwork, &bwork, &errcode);
  R_Free(work);
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
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dgeqrf)(&nrow, &ncol, mat, &ldmat, qraux, work, &lwork, info);
  R_Free(work);
}

void
FM_QL_decomp(double *mat, int ldmat, int nrow, int ncol, double *qlaux, int *info)
{ /* return the QL decomposition of a rectangular matrix */
  int errcode = 0, lwork;
  double opt, *work;

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dgeqlf)(&nrow, &ncol, mat, &ldmat, qlaux, &opt, &lwork, &errcode);
  if (errcode != 0)
    error("DGEQLF in QL decomposition gave error code %d", errcode);

  /* calling DGEQLF with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dgeqlf)(&nrow, &ncol, mat, &ldmat, qlaux, work, &lwork, info);
  R_Free(work);
}

void
FM_LQ_decomp(double *mat, int ldmat, int nrow, int ncol, double *lqaux, int *info)
{ /* return the LQ decomposition of a rectangular matrix */
  int errcode = 0, lwork;
  double opt, *work;

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dgelqf)(&nrow, &ncol, mat, &ldmat, lqaux, &opt, &lwork, &errcode);
  if (errcode != 0)
    error("DGELQF in LQ decomposition gave error code %d", errcode);

  /* calling DGEQRF with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dgelqf)(&nrow, &ncol, mat, &ldmat, lqaux, work, &lwork, info);
  R_Free(work);
}

/* routines for QR, QL and LQ operations */

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
  F77_CALL(dormqr)(side, notrans, &yrow, &nrhs, &nrflc, qr, &ldq, qraux, ymat, &ldy, &opt, &lwork, &errcode FCONE FCONE);
  if (errcode != 0)
    error("DORMQR in QR_qy gave error code %d", errcode);

  /* calling DORMQR with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dormqr)(side, notrans, &yrow, &nrhs, &nrflc, qr, &ldq, qraux, ymat, &ldy, work, &lwork, info FCONE FCONE);
  R_Free(work);
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
  F77_CALL(dormqr)(side, trans, &yrow, &nrhs, &nrflc, qr, &ldq, qraux, ymat, &ldy, &opt, &lwork, &errcode FCONE FCONE);
  if (errcode != 0)
    error("DORMQR in QR_qty gave error code %d", errcode);

  /* calling DORMQR with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dormqr)(side, trans, &yrow, &nrhs, &nrflc, qr, &ldq, qraux, ymat, &ldy, work, &lwork, info FCONE FCONE);
  R_Free(work);
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
FM_QR_store_R(double *qr, int ldq, int ncol, double *Dest, int ldDest)
{ /* copy the R part into Dest, it is assumed that nrow >= ncol */
  FM_cpy_upper(qr, ldq, ncol, Dest, ldDest);
}

void
FM_QL_qy(double *ql, int ldq, int nrow, int ncol, double *qlaux, double *ymat, int ldy, int yrow, int ycol, int *info)
{ /* ymat <- ql.qy(ql, ymat) */
  char *side = "L", *notrans = "N";
  int errcode = 0, lwork, nrflc, nrhs;
  double opt, *work;

  nrhs  = ycol;
  nrflc = ncol; /* length of 'qlaux' */

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dormql)(side, notrans, &yrow, &nrhs, &nrflc, ql, &ldq, qlaux, ymat, &ldy, &opt, &lwork, &errcode FCONE FCONE);
  if (errcode != 0)
    error("DORMQL in QL_qy gave error code %d", errcode);

  /* calling DORMQL with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dormql)(side, notrans, &yrow, &nrhs, &nrflc, ql, &ldq, qlaux, ymat, &ldy, work, &lwork, info FCONE FCONE);
  R_Free(work);
}

void
FM_QL_qty(double *ql, int ldq, int nrow, int ncol, double *qlaux, double *ymat, int ldy, int yrow, int ycol, int *info)
{ /* ymat <- ql.qty(ql, ymat) */
  char *side = "L", *trans = "T";
  int errcode = 0, lwork, nrflc, nrhs;
  double opt, *work;

  nrhs  = ycol;
  nrflc = ncol; /* length of 'qlaux' */

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dormql)(side, trans, &yrow, &nrhs, &nrflc, ql, &ldq, qlaux, ymat, &ldy, &opt, &lwork, &errcode FCONE FCONE);
  if (errcode != 0)
    error("DORMQL in QL_qty gave error code %d", errcode);

  /* calling DORMQL with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dormql)(side, trans, &yrow, &nrhs, &nrflc, ql, &ldq, qlaux, ymat, &ldy, work, &lwork, info FCONE FCONE);
  R_Free(work);
}

void
FM_LQ_yq(double *lq, int ldl, int nrow, int ncol, double *lqaux, double *ymat, int ldy, int yrow, int ycol, int *info)
{ /* ymat <- lq.yq(lq, ymat) */
  char *side = "R", *notrans = "T";
  int errcode = 0, lwork, nrflc, nrhs;
  double opt, *work;

  nrhs  = ycol;
  nrflc = MIN(nrow, ncol); /* length of 'lqaux' */

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dormlq)(side, notrans, &yrow, &nrhs, &nrflc, lq, &ldl, lqaux, ymat, &ldy, &opt, &lwork, &errcode FCONE FCONE);
  if (errcode != 0)
    error("DORMLQ in LQ_yq gave error code %d", errcode);

  /* calling DORMLQ with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dormlq)(side, notrans, &yrow, &nrhs, &nrflc, lq, &ldl, lqaux, ymat, &ldy, work, &lwork, info FCONE FCONE);
  R_Free(work);
}

void
FM_LQ_yqt(double *lq, int ldl, int nrow, int ncol, double *lqaux, double *ymat, int ldy, int yrow, int ycol, int *info)
{ /* ymat <- lq.yqt(lq, ymat) */
  char *side = "R", *trans = "T";
  int errcode = 0, lwork, nrflc, nrhs;
  double opt, *work;

  nrhs  = ycol;
  nrflc = MIN(nrow, ncol); /* length of 'lqaux' */

  /* ask for optimal size of work array */
  lwork = -1;
  F77_CALL(dormlq)(side, trans, &yrow, &nrhs, &nrflc, lq, &ldl, lqaux, ymat, &ldy, &opt, &lwork, &errcode FCONE FCONE);
  if (errcode != 0)
    error("DORMLQ in LQ_yqt gave error code %d", errcode);

  /* calling DORMLQ with optimal size of working array */
  lwork = (int) opt;
  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dormlq)(side, trans, &yrow, &nrhs, &nrflc, lq, &ldl, lqaux, ymat, &ldy, work, &lwork, info FCONE FCONE);
  R_Free(work);
}

void
FM_LQ_store_L(double *lq, int ldl, int nrow, double *Dest, int ldDest)
{ /* copy the L part into Dest, it is assumed that nrow <= ncol */
  FM_cpy_lower(lq, ldl, nrow, Dest, ldDest);
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
  F77_CALL(dgels)(notrans, &n, &n, &n, a, &lda, dummy, &n, &opt, &lwork, &errcode FCONE);
  if (errcode != 0)
    error("DGELS in invert_mat gave error code %d", errcode);

  /* calling DGELS with optimal size of working array */
  lwork = (int) opt;
  work  = (double *) R_Calloc(lwork, double);
  b     = (double *) R_Calloc(n * n, double);
  for (int j = 0; j < n; j++)
    b[j * (n + 1)] = 1.0;
  F77_CALL(dgels)(notrans, &n, &n, &n, a, &lda, b, &n, work, &lwork, info FCONE);
  Memcpy(a, b, n * n);
  R_Free(b); R_Free(work);
}

void
FM_invert_triangular(double *a, int lda, int n, int job, int *info)
{ /* computes the inverse of an upper (job = 1) or lower (job = 0)
   * triangular matrix in place */
  char *diag = "N", *uplo;

  uplo = (job) ? "U" : "L";
  F77_CALL(dtrtri)(uplo, diag, &n, a, &lda, info FCONE FCONE);
}

void
FM_chol_inverse(double *a, int lda, int p, int job, int *info)
{ /* inverse matrix from cholesky (or QR) decomposition
   * the lower triangular factor (job = 0), or the upper triangular
   * factor (job = 1) is used */
  char *uplo;

  uplo = (job) ? "U" : "L";
  F77_CALL(dpotri)(uplo, &p, a, &lda, info FCONE);
  /* copying triangular (lower/upper) part */
  if (job)
    FM_cpy_upper2lower(a, lda, p, a);
  else
    FM_cpy_lower2upper(a, lda, p, a);
}

void
FM_lu_inverse(double *a, int lda, int p, int *pivot, int *info)
{ /* computes the inverse of a matrix using the LU factorization */
  int lwork = p;
  double *work;

  work = (double *) R_Calloc(lwork, double);
  F77_CALL(dgetri)(&p, a, &lda, pivot, work, &lwork, info);
  R_Free(work);
}

void
FM_backsolve(double *r, int ldr, int n, double *b, int ldb, int nrhs, int *info)
{ /* backsolve solve triangular systems of the form r %*% x = b, where r
   * is an upper triangular matrix and b is a matrix containing the right-hand
   * sides to equations */
  char *diag = "N", *uplo = "U", *notrans = "N";

  F77_CALL(dtrtrs)(uplo, notrans, diag, &n, &nrhs, r, &ldr, b, &ldb, info FCONE FCONE FCONE);
}

void
FM_forwardsolve(double *l, int ldl, int n, double *b, int ldb, int nrhs, int *info)
{ /* forwardsolve solve triangular systems of the form l %*% x = b, where l
   * is a lower triangular matrix and b is a matrix containing the right-hand
   * sides to equations */
  char *diag = "N", *uplo = "L", *notrans = "N";

  F77_CALL(dtrtrs)(uplo, notrans, diag, &n, &nrhs, l, &ldl, b, &ldb, info FCONE FCONE FCONE);
}

void
FM_lu_solve(double *a, int lda, int p, int *pivot, double *b, int ldb, int nrhs, int *info)
{ /* solves a system of linear equations */
  char *notrans = "N";

  F77_CALL(dgetrs)(notrans, &p, &nrhs, a, &lda, pivot, b, &ldb, info FCONE);
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
