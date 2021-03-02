/* $ID: blas_3.c, last updated 2020-09-10, F.Osorio */

#include "fastmatrix.h"

/* BLAS level 3 wrappers */

void
BLAS3_gemm(double alpha, double *a, int lda, double *b, int ldb, int m, int n, int k,
  char *transa, char *transb, double beta, double *y, int ldy)
{ /* y <- alpha * op(a) %*% op(b) + beta * y,
   * with op(x) = x, or op(x) = t(x) */
  F77_CALL(dgemm)(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, y, &ldy);
}

void
BLAS3_symm(double alpha, double *a, int lda, double *b, int ldb, int nrow, int ncol,
  char *side, char *uplo, double beta, double *y, int ldy)
{ /* y <- alpha * a %*% b + beta * y, or
     y <- alpha * b %*% a + beta * y, with 'a' symmetric matrix */
  F77_CALL(dsymm)(side, uplo, &nrow, &ncol, &alpha, a, &lda, b, &ldb, &beta, y, &ldy);
}

void
BLAS3_syrk(double alpha, double *a, int lda, int n, int k, char *uplo, char *trans,
  double beta, double *y, int ldy)
{ /* y <- alpha * a %*% t(a) + beta * y, or
     y <- alpha * t(a) %*% a + beta * y, with 'y' symmetric matrix */
  F77_CALL(dsyrk)(uplo, trans, &n, &k, &alpha, a, &lda, &beta, y, &ldy);
}

void
BLAS3_trmm(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo,
  char *trans, char *diag, double *y, int ldy)
{ /* y <- alpha * op(a) %*% y, or y <- alpha * y %*% op(a),
   * with op(x) = x, or op(x) = t(x), with 'a' upper or lower triangular matrix */
  F77_CALL(dtrmm)(side, uplo, trans, diag, &nrow, &ncol, &alpha, a, &lda, y, &ldy);
}

void
BLAS3_trsm(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo,
  char *trans, char *diag, double *y, int ldy)
{ /* solve triangular systems:
     solve(op(a), alpha * y), or alpha * y %*% solve(t(a)),
     with 'a' upper or lower triangular matrix, solution is overwritten on 'y' */
  F77_CALL(dtrsm)(side, uplo, trans, diag, &nrow, &ncol, &alpha, a, &lda, y, &ldy);
}
