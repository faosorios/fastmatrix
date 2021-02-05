/* $ID: fastmatrix_API.h, last updated 2020-09-05, F.Osorio */

#ifndef FASTMATRIX_API_H
#define FASTMATRIX_API_H

#include "fastmatrix.h"
#include <R.h>
#include <Rconfig.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <R_ext/Rdynload.h>

/* Rdynload defines typedef void * (*DL_FUNC) (), which is just right
 * for almost all the routines that return void. We add two more */
typedef int (*INT_FUNC)();           /* pointer to a function which returns an int */
typedef double (*DBL_FUNC)();        /* pointer to a function which returns a double */
typedef double (*f)(double, void *); /* pointer to a function, used by 'FM_brent' */

/* BLAS-1: external API */

void BLAS1_axpy(double alpha, double *x, int incx, double *y, int incy, int n) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS1_axpy");
  fun(alpha, x, incx, y, incy, n);
}

void BLAS1_copy(double *y, int incy, double *x, int incx, int n) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS1_copy");
  fun(y, incy, x, incx, n);
}

double BLAS1_dot_product(double *x, int incx, double *y, int incy, int n) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "BLAS1_dot_product");
    if (fun == NULL) Rf_error("cannot find function 'BLAS1_dot_product'");
  }
  return(fun(x, incx, y, incy, n));
}

int BLAS1_index_max(double *x, int inc, int n) {
  static INT_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (INT_FUNC) R_GetCCallable("fastmatrix", "BLAS1_index_max");
    if (fun == NULL) Rf_error("cannot find function 'BLAS1_index_max'");
  }
  return(fun(x, inc, n));
}

double BLAS1_norm_two(double *x, int inc, int n) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "BLAS1_norm_two");
    if (fun == NULL) Rf_error("cannot find function 'BLAS1_norm_two'");
  }
  return(fun(x, inc, n));
}

void BLAS1_scale(double alpha, double *x, int inc, int n) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS1_scale");
  fun(alpha, x, inc, n);
}

double BLAS1_sum_abs(double *x, int inc, int n) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "BLAS1_sum_abs");
    if (fun == NULL) Rf_error("cannot find function 'BLAS1_sum_abs'");
  }
  return(fun(x, inc, n));
}

void BLAS1_swap(double *x, int incx, double *y, int incy, int n) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS1_swap");
  fun(x, incx, y, incy, n);
}

/* BLAS-2: external API */

void BLAS2_gemv(double alpha, double *a, int lda, int nrow, int ncol, char *trans,
  double *x, int incx, double beta, double *y, int incy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS2_gemv");
  fun(alpha, a, lda, nrow, ncol, trans, x, incx, beta, y, incy);
}

void BLAS2_symv(double alpha, double *a, int lda, int n, char *uplo, double *x,
  int incx, double beta, double *y, int incy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS2_symv");
  fun(alpha, a, lda, n, uplo, x, incx, beta, y, incy);
}

void BLAS2_trmv(double *a, int lda, int n, char *uplo, char *trans, char *diag,
  double *x, int inc) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS2_trmv");
  fun(a, lda, n, uplo, trans, diag, x, inc);
}

void BLAS2_trsv(double *a, int lda, int n, char *uplo, char *trans, char *diag,
  double *x, int inc) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS2_trsv");
  fun(a, lda, n, uplo, trans, diag, x, inc);
}

void BLAS2_ger(double alpha, double *a, int lda, int nrow, int ncol, double *x,
  int incx, double *y, int incy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS2_ger");
  fun(alpha, a, lda, nrow, ncol, x, incx, y, incy);
}

void BLAS2_syr(double alpha, double *a, int lda, int n, char *uplo, double *x, int inc) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS2_syr");
  fun(alpha, a, lda, n, uplo, x, inc);
}

void BLAS2_syr2(double alpha, double *a, int lda, int n, char *uplo, double *x,
  int incx, double *y, int incy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS2_syr2");
  fun(alpha, a, lda, n, uplo, x, incx, y, incy);
}

/* BLAS-3: external API */

void BLAS3_gemm(double alpha, double *a, int lda, double *b, int ldb, int m, int n, int k, char *transa, char *transb, double beta, double *y, int ldy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS3_gemm");
  fun(alpha, a, lda, b, ldb, m, n, k, transa, transb, beta, y, ldy);
}

void BLAS3_symm(double alpha, double *a, int lda, double *b, int ldb, int nrow, int ncol, char *side, char *uplo, double beta, double *y, int ldy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS3_symm");
  fun(alpha, a, lda, b, ldb, nrow, ncol, side, uplo, beta, y, ldy);
}

void BLAS3_syrk(double alpha, double *a, int lda, int n, int k, char *uplo, char *trans, double beta, double *y, int ldy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS3_syrk");
  fun(alpha, a, lda, n, k, uplo, trans, beta, y, ldy);
}

void BLAS3_trmm(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo, char *trans, char *diag, double *y, int ldy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS3_trmm");
  fun(alpha, a, lda, nrow, ncol, side, uplo, trans, diag, y, ldy);
}

void BLAS3_trsm(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo, char *trans, char *diag, double *y, int ldy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "BLAS3_trsm");
  fun(alpha, a, lda, nrow, ncol, side, uplo, trans, diag, y, ldy);
}

/*  operations on vectors: external API */

double FM_norm_sqr(double *x, int inc, int n) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_norm_sqr");
    if (fun == NULL) Rf_error("cannot find function 'FM_norm_sqr'");
  }
  return(fun(x, inc, n));
}

void FM_normalize(double *x, int inc, int n) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_normalize");
  fun(x, inc, n);
}

double FM_vecsum(double *x, int inc, int n) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_vecsum");
    if (fun == NULL) Rf_error("cannot find function 'FM_vecsum'");
  }
  return(fun(x, inc, n));
}

/* basic matrix manipulations: external API */

void FM_add_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_add_mat");
  fun(y, ldy, alpha, x, ldx, nrow, ncol);
}

void FM_copy_mat(double *y, int ldy, double *x, int ldx, int nrow, int ncol) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_copy_mat");
  fun(y, ldy, x, ldx, nrow, ncol);
}

void FM_copy_trans(double *y, int ldy, double *x, int ldx, int nrow, int ncol) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_copy_trans");
  fun(y, ldy, x, ldx, nrow, ncol);
}

void FM_crossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_crossprod");
  fun(z, x, ldx, xrows, xcols, y, ldy, yrows, ycols);
}

void FM_GAXPY(double *y, double alpha, double *a, int lda, int nrow, int ncol, double *x, double beta, int job) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_GAXPY");
  fun(y, alpha, a, lda, nrow, ncol, x, beta, job);
}

double FM_logAbsDet(double *a, int lda, int n) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_logAbsDet");
    if (fun == NULL) Rf_error("cannot find function 'FM_logAbsDet'");
  }
  return(fun(a, lda, n));
}

void FM_mult_mat(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_mult_mat");
  fun(z, x, ldx, xrows, xcols, y, ldy, yrows, ycols);
}

void FM_mult_triangular(double *y, double *a, int lda, int n, double *x, int job) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_mult_triangular");
  fun(y, a, lda, n, x, job);
}

void FM_rank1_update(double *a, int lda, int nrow, int ncol, double alpha, double *x, double *y) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_rank1_update");
  fun(a, lda, nrow, ncol, alpha, x, y);
}

void FM_scale_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_scale_mat");
  fun(y, ldy, alpha, x, ldx, nrow, ncol);
}

void FM_setzero(double *y, int ldy, int nrow, int ncol) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_setzero");
  fun(y, ldy, nrow, ncol);
}

double FM_trace(double *a, int lda, int n) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_trace");
    if (fun == NULL) Rf_error("cannot find function 'FM_trace'");
  }
  return(fun(a, lda, n));
}

void FM_tcrossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_tcrossprod");
  fun(z, x, ldx, xrows, xcols, y, ldy, yrows, ycols);
}

/* operations on triangular matrices: external API */

void FM_cpy_lower(double *x, int ldx, int p, double *y, int ldy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_cpy_lower");
  fun(x, ldx, p, y, ldy);
}

void FM_cpy_upper(double *x, int ldx, int p, double *y, int ldy) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_cpy_upper");
  fun(x, ldx, p, y, ldy);
}

void FM_cpy_lower2upper(double *x, int ldx, int p, double *y) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_cpy_lower2upper");
  fun(x, ldx, p, y);
}

void FM_cpy_upper2lower(double *x, int ldx, int p, double *y) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_cpy_upper2lower");
  fun(x, ldx, p, y);
}

double FM_sum_lower_tri(double *x, int ldx, int p, int job) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_sum_lower_tri");
    if (fun == NULL) Rf_error("cannot find function 'FM_sum_lower_tri'");
  }
  return(fun(x, ldx, p, job));
}

double FM_sum_upper_tri(double *x, int ldx, int p, int job) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_sum_upper_tri");
    if (fun == NULL) Rf_error("cannot find function 'FM_sum_upper_tri'");
  }
  return(fun(x, ldx, p, job));
}

/* matrix factorizations: external API */

void FM_chol_decomp(double *a, int lda, int p, int job, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_chol_decomp");
  fun(a, lda, p, job, info);
}

void FM_QR_decomp(double *mat, int ldmat, int nrow, int ncol, double *qraux, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_QR_decomp");
  fun(mat, ldmat, nrow, ncol, qraux, info);
}

void FM_QL_decomp(double *mat, int ldmat, int nrow, int ncol, double *qlaux, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_QL_decomp");
  fun(mat, ldmat, nrow, ncol, qlaux, info);
}

void FM_LQ_decomp(double *mat, int ldmat, int nrow, int ncol, double *lqaux, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_LQ_decomp");
  fun(mat, ldmat, nrow, ncol, lqaux, info);
}

void FM_svd_decomp(double *mat, int ldmat, int nrow, int ncol, double *u, int ldu, double *d, double *v, int ldv, int job, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_svd_decomp");
  fun(mat, ldmat, nrow, ncol, u, ldu, d, v, ldv, job, info);
}

/* QR operations: external API */

void FM_QR_qy(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_QR_qy");
  fun(qr, ldq, nrow, ncol, qraux, ymat, ldy, yrow, ycol, info);
}

void FM_QR_qty(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_QR_qty");
  fun(qr, ldq, nrow, ncol, qraux, ymat, ldy, yrow, ycol, info);
}

void FM_QR_fitted(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int job, double *fitted) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_QR_fitted");
  fun(qr, ldq, nrow, ncol, qraux, ymat, ldy, yrow, ycol, job, fitted);
}

void FM_QR_store_R(double *qr, int ldq, int nrow, int ncol, double *Dest, int ldDest) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_QR_store_R");
  fun(qr, ldq, nrow, ncol, Dest, ldDest);
}

void FM_QL_qy(double *ql, int ldq, int nrow, int ncol, double *qlaux, double *ymat, int ldy, int yrow, int ycol, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_QL_qy");
  fun(ql, ldq, nrow, ncol, qlaux, ymat, ldy, yrow, ycol, info);
}

void FM_QL_qty(double *ql, int ldq, int nrow, int ncol, double *qlaux, double *ymat, int ldy, int yrow, int ycol, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_QL_qty");
  fun(ql, ldq, nrow, ncol, qlaux, ymat, ldy, yrow, ycol, info);
}

void FM_LQ_yq(double *lq, int ldl, int nrow, int ncol, double *lqaux, double *ymat, int ldy, int yrow, int ycol, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_LQ_yq");
  fun(lq, ldl, nrow, ncol, lqaux, ymat, ldy, yrow, ycol, info);
}

void FM_LQ_yqt(double *lq, int ldl, int nrow, int ncol, double *lqaux, double *ymat, int ldy, int yrow, int ycol, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_LQ_yqt");
  fun(lq, ldl, nrow, ncol, lqaux, ymat, ldy, yrow, ycol, info);
}

void FM_LQ_store_L(double *lq, int ldq, int nrow, double *Dest, int ldDest) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_LQ_store_L");
  fun(lq, ldq, nrow, Dest, ldDest);
}

/* matrix inversion and linear solvers: external API */

void FM_backsolve(double *r, int ldr, int n, double *b, int ldb, int nrhs, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_backsolve");
  fun(r, ldr, n, b, ldb, nrhs, info);
}

void FM_forwardsolve(double *l, int ldl, int n, double *b, int ldb, int nrhs, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_forwardsolve");
  fun(l, ldl, n, b, ldb, nrhs, info);
}

void FM_chol_inverse(double *a, int lda, int p, int job, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_chol_inverse");
  fun(a, lda, p, job, info);
}

void FM_invert_mat(double *a, int lda, int n, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_invert_mat");
  fun(a, lda, n, info);
}

void FM_invert_triangular(double *a, int lda, int n, int job, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_invert_triangular");
  fun(a, lda, n, job, info);
}

/* least square procedures: external API */

void FM_lsfit(double *x, int ldx, int nrow, int ncol, double *y, int ldy, int nrhs, double *coef, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_lsfit");
  fun(x, ldx, nrow, ncol, y, ldy, nrhs, coef, info);
}

void FM_gls_GQR(double *x, int ldx, int nrow, int ncol, double *y, double *cov, double *coef, int *info) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_gls_GQR");
  fun(x, ldx, nrow, ncol, y, cov, coef, info);
}

/* Distances: external API */

double FM_pythag(double a, double b) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_pythag");
    if (fun == NULL) Rf_error("cannot find function 'FM_pythag'");
  }
  return(fun(a, b));
}

double FM_mahalanobis(double *x, int p, double *center, double *Root) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_mahalanobis");
    if (fun == NULL) Rf_error("cannot find function 'FM_mahalanobis'");
  }
  return(fun(x, p, center, Root));
}

void FM_mahal_distances(double *x, int n, int p, double *center, double *cov, int inverted, double *distances) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_mahal_distances");
  fun(x, n, p, center, cov, inverted, distances);
}

void FM_WH_chisq(double *distances, int n, int p, double *z) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_WH_chisq");
  fun(distances, n, p, z);
}

void FM_WH_F(double *distances, int n, int p, double eta, double *z) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_WH_F");
  fun(distances, n, p, eta, z);
}

/* Descriptive statistics: external API */

void FM_mean_and_var(double *x, int nobs, double *mean, double *var) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_mean_and_var");
  fun(x, nobs, mean, var);
}

void FM_online_covariance(double *x, double *y, int nobs, double *xbar, double *ybar, double *xvar, double *yvar, double *cov) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_online_covariance");
  fun(x, y, nobs, xbar, ybar, xvar, yvar, cov);
}

void FM_online_center(double *x, int n, int p, double *weights, double *center) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_online_center");
  fun(x, n, p, weights, center);
}

void FM_center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_center_and_Scatter");
  fun(x, n, p, weights, center, Scatter);
}

void FM_skewness_and_kurtosis(double *x, int n, int p, double *center, double *Scatter, double *stats, int do_skewness) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_skewness_and_kurtosis");
  fun(x, n, p, center, Scatter, stats, do_skewness);
}

void FM_cov_MSSD(double *x, int n, int p, double *center, double *Scatter) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_cov_MSSD");
  fun(x, n, p, center, Scatter);
}

double FM_find_quantile(double *a, int n, int k) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_find_quantile");
    if (fun == NULL) Rf_error("cannot find function 'FM_find_quantile'");
  }
  return(fun(a, n, k));
}

/* Brent's method for unidimensional optimization: external API */

double FM_brent(double ax, double bx, double (*f)(double, void *), void *info, double tolerance) {
  static DBL_FUNC fun = NULL;
  if (fun == NULL) {
    fun = (DBL_FUNC) R_GetCCallable("fastmatrix", "FM_brent");
    if (fun == NULL) Rf_error("cannot find function 'FM_brent'");
  }
  return(fun(ax, bx, f, info, tolerance));
}

/* Misc: external API */

void FM_centering(double *x, int n, int p, double *center) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_centering");
  fun(x, n, p, center);
}

void FM_cov2cor(double *cov, int p) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_cov2cor");
  fun(cov, p);
}

void FM_sherman_morrison(double *a, int lda, int n, double *b, double *d, int inverted) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_sherman_morrison");
  fun(a, lda, n, b, d, inverted);
}

/* 'DEBUG' routine: external API */

void FM_print_mat(double *x, int ldx, int nrow, int ncol, char *msg) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_print_mat");
  fun(x, ldx, nrow, ncol, msg);
}

#endif /* FASTMATRIX_API_H */
