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

/* Rdynload defines "typedef void * (*DL_FUNC) (), which is just right
 * for almost all the routines that return void. We add two more */
typedef int (*INT_FUNC)();      /* pointer to a function which returns an int */
typedef double (*DBL_FUNC)();   /* pointer to a function which returns a double */

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

void FM_center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_center_and_Scatter");
  fun(x, n, p, weights, center, Scatter);
}

void FM_MSSD(double *x, int n, int p, double *center, double *Scatter) {
  static void (*fun)() = NULL;
  if (fun == NULL)
    fun = (void (*)) R_GetCCallable("fastmatrix", "FM_MSSD");
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

#endif /* FASTMATRIX_API_H */
