/* ID: fastmatrix.h, last updated 2020-09-10, F.Osorio */

#ifndef FASTMATRIX_H
#define FASTMATRIX_H

/* This contains the prototype calls for all the .C functions that
 * are called by another C function */

/* BLAS-1 wrappers */
void BLAS1_axpy(double alpha, double *x, int incx, double *y, int incy, int n);
void BLAS1_copy(double *y, int incy, double *x, int incx, int n);
double BLAS1_dot_product(double *x, int incx, double *y, int incy, int n);
int BLAS1_index_max(double *x, int inc, int n);
double BLAS1_norm_two(double *x, int inc, int n);
void BLAS1_scale(double alpha, double *x, int inc, int n);
double BLAS1_sum_abs(double *x, int inc, int n);
void BLAS1_swap(double *x, int incx, double *y, int incy, int n);

/* BLAS-2 wrappers */
void BLAS2_gemv(double alpha, double *a, int lda, int nrow, int ncol, char *trans, double *x, int incx, double beta, double *y, int incy);
void BLAS2_symv(double alpha, double *a, int lda, int n, char *uplo, double *x, int incx, double beta, double *y, int incy);
void BLAS2_trmv(double *a, int lda, int n, char *uplo, char *trans, char *diag, double *x, int inc);
void BLAS2_trsv(double *a, int lda, int n, char *uplo, char *trans, char *diag, double *x, int inc);
void BLAS2_ger(double alpha, double *a, int lda, int nrow, int ncol, double *x, int incx, double *y, int incy);
void BLAS2_syr(double alpha, double *a, int lda, int n, char *uplo, double *x, int inc);
void BLAS2_syr2(double alpha, double *a, int lda, int n, char *uplo, double *x, int incx, double *y, int incy);

/* BLAS-3 wrappers */
void BLAS3_gemm(double alpha, double *a, int lda, double *b, int ldb, int m, int n, int k, char *transa, char *transb, double beta, double *y, int ldy);
void BLAS3_symm(double alpha, double *a, int lda, double *b, int ldb, int nrow, int ncol, char *side, char *uplo, double beta, double *y, int ldy);
void BLAS3_syrk(double alpha, double *a, int lda, int n, int k, char *uplo, char *trans, double beta, double *y, int ldy);
void BLAS3_trmm(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo, char *trans, char *diag, double *y, int ldy);
void BLAS3_trsm(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo, char *trans, char *diag, double *y, int ldy);

/* descriptive statistics */
void FM_mean_and_var(double *x, int nobs, double *mean, double *var);
void FM_online_covariance(double *x, double *y, int nobs, double *xbar, double *ybar, double *xvar, double *yvar, double *cov);
void FM_center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter);
void FM_MSSD(double *x, int n, int p, double *center, double *Scatter);
double FM_find_quantile(double *a, int n, int k);

#endif /* FASTMATRIX_H */
