/* ID: fastmatrix.h, last updated 2022-02-10, F.Osorio */

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
void BLAS1_rot(double *x, int incx, double *y, int incy, int n, double c, double s);
void BLAS1_rotg(double *a, double *b, double *c, double *s);
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

/*  operations on vectors */
double FM_norm_sqr(double *x, int inc, int n);
void FM_normalize(double *x, int inc, int n);
double FM_vecsum(double *x, int inc, int n);

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

/* operations on triangular matrices */
void FM_cpy_lower(double *x, int ldx, int p, double *y, int ldy);
void FM_cpy_upper(double *x, int ldx, int p, double *y, int ldy);
void FM_cpy_lower2upper(double *x, int ldx, int p, double *y);
void FM_cpy_upper2lower(double *x, int ldx, int p, double *y);
double FM_sum_lower_tri(double *x, int ldx, int p, int job);
double FM_sum_upper_tri(double *x, int ldx, int p, int job);

/* matrix factorizations */
void FM_chol_decomp(double *a, int lda, int p, int job, int *info);
void FM_QR_decomp(double *mat, int ldmat, int nrow, int ncol, double *qraux, int *info);
void FM_QL_decomp(double *mat, int ldmat, int nrow, int ncol, double *qlaux, int *info);
void FM_LQ_decomp(double *mat, int ldmat, int nrow, int ncol, double *lqaux, int *info);
void FM_svd_decomp(double *mat, int ldmat, int nrow, int ncol, double *u, int ldu, double *d, double *v, int ldv, int job, int *info);

/* QR, QL and LQ operations */
void FM_QR_qy(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info);
void FM_QR_qty(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int *info);
void FM_QR_fitted(double *qr, int ldq, int nrow, int ncol, double *qraux, double *ymat, int ldy, int yrow, int ycol, int job, double *fitted);
void FM_QR_store_R(double *qr, int ldq, int nrow, int ncol, double *Dest, int ldDest);
void FM_QL_qy(double *ql, int ldq, int nrow, int ncol, double *qlaux, double *ymat, int ldy, int yrow, int ycol, int *info);
void FM_QL_qty(double *ql, int ldq, int nrow, int ncol, double *qlaux, double *ymat, int ldy, int yrow, int ycol, int *info);
void FM_LQ_yq(double *lq, int ldl, int nrow, int ncol, double *lqaux, double *ymat, int ldy, int yrow, int ycol, int *info);
void FM_LQ_yqt(double *lq, int ldl, int nrow, int ncol, double *lqaux, double *ymat, int ldy, int yrow, int ycol, int *info);
void FM_LQ_store_L(double *lq, int ldl, int nrow, double *Dest, int ldDest);

/* matrix inversion and linear solvers */
void FM_backsolve(double *r, int ldr, int n, double *b, int ldb, int nrhs, int *info);
void FM_forwardsolve(double *l, int ldl, int n, double *b, int ldb, int nrhs, int *info);
void FM_chol_inverse(double *a, int lda, int p, int job, int *info);
void FM_invert_mat(double *a, int lda, int n, int *info);
void FM_invert_triangular(double *a, int lda, int n, int job, int *info);

/* least square procedures */
void FM_lsfit(double *x, int ldx, int nrow, int ncol, double *y, int ldy, int nrhs, double *coef, int *info);
void FM_gls_GQR(double *x, int ldx, int nrow, int ncol, double *y, double *cov, double *coef, int *info);

/* distances */
double FM_pythag(double a, double b);
double FM_mahalanobis(double *x, int p, double *center, double *Root);
void FM_mahal_distances(double *x, int n, int p, double *center, double *cov, int inverted, double *distances);
void FM_WH_chisq(double *distances, int n, int p, double *z);
void FM_WH_F(double *distances, int n, int p, double eta, double *z);

/* products */
void FM_two_product_FMA(double a, double b, double *x, double *y);
void FM_compensated_product(double *x, int nobs, double *prod);

/* descriptive statistics */
void FM_center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter);
void FM_cov_MSSD(double *x, int n, int p, double *center, double *Scatter);
double FM_find_quantile(double *a, int n, int k);
void FM_geometric_mean(double *x, int nobs, double *mean);
void FM_mean_and_var(double *x, int nobs, double *mean, double *var);
void FM_moments(double *x, int nobs, double *mean, double *s2, double *s3, double *s4);
void FM_online_center(double *x, int n, int p, double *weights, double *center);
void FM_online_covariance(double *x, double *y, int nobs, double *xbar, double *ybar, double *xvar, double *yvar, double *cov);
void FM_skewness_and_kurtosis(double *x, int n, int p, double *center, double *Scatter, double *stats, int do_skewness);

/* misc */
void FM_centering(double *x, int n, int p, double *center);
void FM_cov2cor(double *cov, int p);
void FM_sherman_morrison(double *a, int lda, int n, double *b, double *d, int inverted);

/* 'DEBUG' routine */
void FM_print_mat(double *x, int ldx, int nrow, int ncol, char *msg);

#endif /* FASTMATRIX_H */
