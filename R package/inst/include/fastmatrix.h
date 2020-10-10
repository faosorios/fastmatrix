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

/* basic matrix manipulations */
void FM_add_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol);
void FM_copy_mat(double *y, int ldy, double *x, int ldx, int nrow, int ncol);
void FM_copy_trans(double *y, int ldy, double *x, int ldx, int nrow, int ncol);
void FM_crossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols);
void FM_GAXPY(double *y, double alpha, double *a, int lda, int nrow, int ncol, double *x, double beta, int job);
double FM_logAbsDet(double *a, int lda, int n);
void FM_mult_mat(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols);
void FM_mult_triangular(double *y, double *a, int lda, int n, double *x, int job);
void FM_rank1_update(double *a, int lda, int nrow, int ncol, double alpha, double *x, double *y);
void FM_scale_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol);
void FM_setzero(double *y, int ldy, int nrow, int ncol);
double FM_trace(double *a, int lda, int n);
void FM_tcrossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols);

/* matrix factorizations */
void FM_chol_decomp(double *a, int lda, int p, int job, int *info);
void FM_QR_decomp(double *mat, int ldmat, int nrow, int ncol, double *qraux, int *info);
void FM_svd_decomp(double *mat, int ldmat, int nrow, int ncol, double *u, int ldu, double *d, double *v, int ldv, int job, int *info);

/* QR operations */
void FM_QR_qy(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info);
void FM_QR_qty(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info);
void FM_QR_fitted(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int job, double *fitted);
void FM_QR_store_R(double *qr, int ldq, int nrow, int ncol, double *Dest, int ldDest);

/* matrix inversion and linear solvers */
void FM_backsolve(double *r, int ldr, int n, double *b, int ldb, int nrhs, int *info);
void FM_forwardsolve(double *l, int ldl, int n, double *b, int ldb, int nrhs, int *info);
void FM_chol_inverse(double *a, int lda, int p, int job, int *info);
void FM_invert_mat(double *a, int lda, int n, int *info);
void FM_invert_triangular(double *a, int lda, int n, int job, int *info);

/* descriptive statistics */
void FM_mean_and_var(double *x, int nobs, double *mean, double *var);
void FM_online_covariance(double *x, double *y, int nobs, double *xbar, double *ybar, double *xvar, double *yvar, double *cov);
void FM_center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter);
void FM_cov_MSSD(double *x, int n, int p, double *center, double *Scatter);
double FM_find_quantile(double *a, int n, int k);

/* misc */
void FM_cov2cor(double *cov, int p);

/* 'DEBUG' routine */
void FM_print_mat(double *x, int ldx, int nrow, int ncol, char *msg);

#endif /* FASTMATRIX_H */
