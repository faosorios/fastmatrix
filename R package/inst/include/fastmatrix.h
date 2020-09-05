/* ID: fastmatrix.h, last updated 2020-09-05, F.Osorio */

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

/* center and Scatter estimation */
void center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter);
void MSSD(double *x, int n, int p, double *center, double *Scatter);

#endif /* FASTMATRIX_H */
