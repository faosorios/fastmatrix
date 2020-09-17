/* ID: fastmatrix.h, last updated 2020-09-04, F.Osorio */

#ifndef FASTMATRIX_H
#define FASTMATRIX_H

#include <R.h>
#include <Rconfig.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

/* some definitions */
#define DNULLP   (double *) 0
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SQR(x)   R_pow_di(x, 2)
#define repeat   for(;;)

/* operations on arrays */
void F77_NAME(arraymult)(double *, int *, int *, int *, double *, int *, int *, int *, double *, int *, int *, double *, int *, int *);
void F77_NAME(bracketprod)(double *, int *, int *, int *, double *, int *, int *, int *, double *, int *, int *);

/* BLAS-1 wrappers */
void BLAS1_axpy(double, double *, int, double *, int, int);
void BLAS1_copy(double *, int, double *, int, int);
double BLAS1_dot_product(double *, int, double *, int, int);
int BLAS1_index_max(double *, int, int);
double BLAS1_norm_two(double *, int, int);
void BLAS1_scale(double, double *, int, int);
double BLAS1_sum_abs(double *, int, int);
void BLAS1_swap(double *, int, double *, int, int);

/* BLAS-2 wrappers */
void BLAS2_gemv(double, double *, int, int, int, char *, double *, int, double, double *, int);
void BLAS2_symv(double, double *, int, int, char *, double *, int, double, double *, int);
void BLAS2_trmv(double *, int, int, char *, char *, char *, double *, int);
void BLAS2_trsv(double *, int, int, char *, char *, char *, double *, int);
void BLAS2_ger(double, double *, int, int, int, double *, int, double *, int);
void BLAS2_syr(double, double *, int, int, char *, double *, int);
void BLAS2_syr2(double, double *, int, int, char *, double *, int, double *, int);

/* BLAS-3 wrappers */
void BLAS3_gemm(double, double *, int, double *, int, int, int, int, char *, char *, double, double *, int);
void BLAS3_symm(double, double *, int, double *, int, int, int, char *, char *, double, double *, int);
void BLAS3_syrk(double, double *, int, int, int , char *, char *, double, double *, int);
void BLAS3_trmm(double, double *, int, int, int, char *, char *, char *, char *, double *, int);
void BLAS3_trsm(double, double *, int, int, int, char *, char *, char *, char *, double *, int);

/* routines for operations on commutation matrices */
void F77_NAME(comm_rows)(int *, int *, int *);
void F77_NAME(commutation_mat)(int *, int *, int *, int *, int *, int *);
void F77_NAME(comm_left_mult)(int *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
void F77_NAME(comm_right_mult)(int *, int *, int *, double *, int *, int *, int *, double *, int *, int *);

/* routines for operations on duplication matrices */
void dupl_cols(int *, int *);
void duplication_mat(int *, int *, int *, int *);
void dupl_left_mult(double *, int *, int *, int *, int *, int *, double *, int *);
void dupl_left_trans(double *, int *, int *, int *, int *, int *, int *, double *, int *);
void dupl_right_mult(double *, int *, int *, int *, int *, int *, int *, double *, int *);
void dupl_right_trans(double *, int *, int *, int *, int *, int *, double *, int *);

/* routines for operations on symmetrizer matrices */
void F77_NAME(symmetrizer_mat)(double *, int *, int *, int *, int *, double *, int *, int *);
void symmetrizer_prod(double *, int *, int *, int *, double *, int *);

/* basic matrix manipulations */
void add_mat(double *, int, double, double *, int, int, int);
void copy_mat(double *, int, double *, int, int, int);
void copy_trans(double *, int, double *, int, int, int);
void scale_mat(double *, int, double, double *, int, int, int);
void setzero(double *, int, int, int);

/* vector norms */
void norm_one(double *, int *, int *, double *);
void norm_two(double *, int *, int *, double *);
void norm_inf(double *, int *, int *, double *);
void norm_minkowski(double *, int *, int *, double *, double *);
void matrix_norm(double *, int *, int *, int *, int *, double *);
double F77_NAME(minkowski)(int *, double *, int *, double *);

/* power method */
void power_method(double *, int *, int *, int *, double *, double *, int *, double *, int *);

/* Sherman-Morrison formula */
void sherman_morrison(double *, int *, int *, double *, double *);

/* matrix factorizations */
void lu_decomp(double *, int *, int *, int *, int *);
void lu_inverse(double *, int *, int *, int *);
void lu_solve(double *, int *, int *, int *, double *, int *, int *);

/* descriptive statistics */
void FM_mean_and_var(double *, int, double *, double *);
void FM_online_covariance(double *, double *, int, double *, double *, double *, double *, double *);
void FM_center_and_Scatter(double *, int, int, double *, double *, double *);
void FM_MSSD(double *, int, int, double *, double *);
double FM_find_quantile(double *, int, int);

/* sweep operator for symmetric matrices */
void sweep_operator(double *, int *, int *, int *, int *, int *);
void F77_NAME(sweepop)(double *, int *, int *, int *, int *, int *);

/* utils on vectors */
double norm_sqr(double *, int, int);
void normalize_vec(double *, int, int);

/* utils on matrices */
void equilibrate_mat(double *, int *, int *, int *, double *, double *, int *);
void F77_NAME(equilibrate_cols)(double *, int *, int *, int *, double *, double *, int *, int *);
void F77_NAME(hadamard_prod)(double *, double *, int *, double *);
void F77_NAME(inner_frobenius)(double *, int *, double *, int *, int *, int *, double *);
void mat2vech(double *, int *, int *, double *);
void F77_NAME(pivot_mat)(double *, int *, int *, int *);

#endif /* FASTMATRIX_H */
