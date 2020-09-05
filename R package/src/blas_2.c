/* $ID: blas_2.c, last updated 2020-09-01, F.Osorio */

#include "fastmatrix.h"

/* BLAS level 2 wrappers */

void
BLAS2_gemv(double alpha, double *a, int lda, int nrow, int ncol, char *trans, double *x,
  int incx, double beta, double *y, int incy)
{ /* y <- alpha * a %*% x + beta * y, or
     y <- alpha * t(a) %*% x + beta * y */
  F77_CALL(dgemv)(trans, &nrow, &ncol, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void
BLAS2_symv(double alpha, double *a, int lda, int n, char *uplo, double *x, int incx,
  double beta, double *y, int incy)
{ /* y <- alpha * a %*% x + beta * y, with 'a' symmetric matrix */
  F77_CALL(dsymv)(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void
BLAS2_trmv(double *a, int lda, int n, char *uplo, char *trans, char *diag, double *x, int inc)
{ /* x <- a %*% x, or x <- t(a) %*% x, with 'a' upper or lower triangular matrix */
  F77_CALL(dtrmv)(uplo, trans, diag, &n, a, &lda, x, &inc);
}

void
BLAS2_trsv(double *a, int lda, int n, char *uplo, char *trans, char *diag, double *x, int inc)
{ /* solve triangular systems:
     solve(a, x), or solve(t(a), x),
     with 'a' upper or lower triangular matrix */
  F77_CALL(dtrsv)(uplo, trans, diag, &n, a, &lda, x, &inc);
}

void
BLAS2_ger(double alpha, double *a, int lda, int nrow, int ncol, double *x, int incx,
  double *y, int incy)
{ /* a <- alpha * x %*% t(y) + a */
  F77_CALL(dger)(&nrow, &ncol, &alpha, x, &incx, y, &incy, a, &lda);
}

void
BLAS2_syr(double alpha, double *a, int lda, int n, char *uplo, double *x, int inc)
{ /* a <- alpha * x %*% t(x) + a, with 'a' symmetric matrix */
  F77_CALL(dsyr)(uplo, &n, &alpha, x, &inc, a, &lda);
}

void
BLAS2_syr2(double alpha, double *a, int lda, int n, char *uplo, double *x, int incx,
  double *y, int incy)
{ /* a <- alpha * x %*% t(y) + alpha * y %*% t(x) + a, with 'a' symmetric matrix */
  F77_CALL(dsyr2)(uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
