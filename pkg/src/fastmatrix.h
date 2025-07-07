/* ID: fastmatrix.h, last updated 2025-05-26, F.Osorio */

#ifndef FASTMATRIX_H
#define FASTMATRIX_H

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
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
#ifndef FCONE
# define FCONE
#endif
#define CUBE(x)         R_pow_di(x, 3)
#define DNULLP          (double *) 0
#define EPS_CONV        1.0e-2
#define EPS_TINY        3.6668528625010358e-11
#define FOURTH(x)       R_pow_di(x, 4)
#define GOLDEN          0.3819660112501051
#define IZERO(x)        (((x) == 0) ? 1 : 0)
#define MAX(a,b)        (((a)>(b)) ? (a) : (b))
#define MIN(a,b)        (((a)<(b)) ? (a) : (b))
#define OFFSET(n, inc)  (((inc) > 0) ? 0 : ((n) - 1) * (-(inc)))
#define repeat          for(;;)
#define SGN(x)          (((x) >= 0) ? 1.0 : -1.0)
#define SQR(x)          R_pow_di(x, 2)

/* operations on arrays */
void F77_NAME(arraymult)(double *, int *, int *, int *, double *, int *, int *, int *, double *, int *, int *, double *, int *, int *);
void F77_NAME(bracketprod)(double *, int *, int *, int *, double *, int *, int *, int *, double *, int *, int *);

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

/* other matrix operations */
double F77_NAME(blinf)(double *, int *, int *, int *, double *, double *);
double F77_NAME(quadf)(double *, int *, int *, double *);
double F77_NAME(murrv)(double *, double *, int *, int *, int *, double *, int *);

/* constructors of special matrices */
void F77_NAME(frank_mat)(double *, int *, int *, int *);
void F77_NAME(helmert_mat)(double *, int *, int *, int *);
void krylov_mat(double *, int *, int *, double *, int *, double *, int *, int *);

/* vector norms */
void norm_one(double *, int *, int *, double *);
void norm_two(double *, int *, int *, double *);
void norm_inf(double *, int *, int *, double *);
void norm_minkowski(double *, int *, int *, double *, double *);
void matrix_norm(double *, int *, int *, int *, int *, double *);
double F77_NAME(minkowski)(int *, double *, int *, double *);

/* equilibration of matrices */
void F77_NAME(equilibrate_cols)(double *, int *, int *, int *, double *, double *, int *, int *);
void equilibrate_pd(double *, int *, int *, double *, double *, double *, int *);
void equilibrate_sym(double *, int *, int *, double *, double *, double *, int *);
void equilibrating_sym(double *, int *, double *);

/* kronecker product */
void kronecker_prod(double *, int *, int *, double *, int *, int *, double *);

/* power method */
void power_method(double *, int *, int *, double *, double *, int *, double *, int *);
void inverse_power(double *, int *, int *, double *, double *, int *, double *, int *);

/* Sherman-Morrison formula */
void sherman_morrison(double *, int *, int *, double *, double *, int *);

/* LU factorization */
void lu_dcmp(double *, int *, int *, int *, int *, int *);
void lu_inverse(double *, int *, int *, int *, int *);
void lu_solve(double *, int *, int *, int *, double *, int *, int *, int *);

/* iterative methods to solve linear systems */
void cg_solver(double *, int *, int *, double *, double *, int *, double *, int *, int *);
void jacobi_solver(double *, int *, int *, double *, double *, int *, double *, int *, int *);
void seidel_solver(double *, int *, int *, double *, double *, int *, double *, int *, int *);

/* matrix decompositions and operations */
void chol_dcmp(double *, int *, int *, int *, int *);
void chol_update(double *, int *, int *, double *);
void F77_NAME(mchol_dcmp)(double *, int *, int *, double *, double *, int *);
void F77_NAME(ldl_dcmp)(double *, int *, int *, double *, int *);
void svd_dcmp(double *, int *, int *, int *, double *, int *, double *, double *, int *, int *, int *);

/* OLS methods */
void OLS_cg(double *, int *, int *, int *, double *, double *, double *, int *, int *);
void OLS_qr(double *, int *, int *, int *, double *, double *, double *, double *, double *, double *);

/* ridge regression */
void OLS_ridge(double *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, double *, int *);

/* descriptive statistics */
void central_moments(double *, int *, double *, double *, double *, double *);
void cov_weighted(double *, int *, int *, double *, double *, double *);
void cov_MSSD(double *, int *, int *, double *, double *);
void cov4th(double *, int *, int *, double *, double *);
void F77_NAME(median_center)(double *, int *, int *, int *, double *, int *, int *);
void geometric_mean(double *, int *, double *);
void mahal_distances(double *, int *, int *, double *, double *, int *, double *);
void skewness_and_kurtosis(double *, int *, int *, double *, double *, double *, int *);

/* Wilson-Hilferty transformation */
void wilson_hilferty_chisq(double *, int *, int *, double *);
void wilson_hilferty_gamma(double *, int *, double *, double *, double *);

/* tests for normality based on the standardized third and fourth moments */
void doornik_hansen(double *, int *, double *, double *, double *);
void jarque_bera(double *, int *, double *, double *, double *);
void robust_JB(double *, double *, int *, double *, double *, double *);
void urzua_ALM(double *, int *, double *, double *, double *);

/* random number generation */
void rng_mnorm(double *, int *, int *, double *, double *);
void rng_ball(double *, int *, int *);
void rng_sphere(double *, int *, int *);

/* correlation structures */
void cor_AR1(double *, int *, double *);
void cor_CS(double *, int *, double *);

/* kernel U-statistics for asymptotic variance estimation of CCC */
void F77_NAME(rhoc_ustat)(double *, double *, int *, double *, double *, double *);

/* sweep operator for symmetric matrices */
void sweep_operator(double *, int *, int *, int *, int *, int *);
void F77_NAME(sweepop)(double *, int *, int *, int *, int *, int *);

/* Brent's method for unidimensional optimization */
double brent(double, double, double (*f)(double, void *), void *, double);

/* Bezier curve smoother */
void bezier_smoother(double *, double *, int *, double *, int *, double *, double *);
void F77_NAME(decasteljau)(double *, double *, int *, double *, double *);

/* misc */
void F77_NAME(circulant_mat)(double *, int *, double *, int *, int *);
void hadamard_prod(double *, double *, int *, double *);
void F77_NAME(hankel_mat)(double *, double *, int *, double *, int *, int *);
void F77_NAME(inner_frobenius)(double *, int *, double *, int *, int *, int *, double *);
void mat2vech(double *, int *, int *, double *);
void matrix_polynomial(double *, int *, int *, double *, int *, double *, int *, int *);
void F77_NAME(pivot_mat)(double *, int *, int *, int *);
void Psi2Q(double *, double *, int *);
void sqrt_mat_DB(double *, int *, int *, int *, int *, double *, int *);
void whitening_chol(double *, int *, int *, double *);

/* ========================================================================== *
 * additional symbols from LAPACK
 * ========================================================================== */

/* DSYEQUB - compute row and column scalings intended to equilibrate
 * a symmetric matrix A (with respect to the two-norm) and reduce
 * its condition number */
extern void
F77_NAME(dsyequb)(const char *uplo, const int *n, const double *a, const int *lda,
    double *s, double *scond, double *amax, double *work, int *info FCLEN);

/* ========================================================================== *
 * symbols callable from other packages
 * ========================================================================== */

/* BLAS-1 wrappers */
void BLAS1_axpy(double, double *, int, double *, int, int);
void BLAS1_copy(double *, int, double *, int, int);
double BLAS1_dot_product(double *, int, double *, int, int);
int BLAS1_index_max(double *, int, int);
double BLAS1_norm_two(double *, int, int);
void BLAS1_rot(double *, int, double *, int, int, double, double);
void BLAS1_rotg(double *, double *, double *, double *);
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

/* other matrix operations */
double OMO_blinf(double *, int, int, int, double *, double *);
double OMO_quadf(double *, int, int, double *);
void OMO_murrv(double *, double *, int, int, int, double *, int *);

/* operations on vectors */
double FM_norm_sqr(double *, int, int);
void FM_normalize(double *, int, int);
double FM_vecsum(double *, int, int);

/* basic matrix manipulations */
void FM_add_mat(double *, int, double, double *, int, int, int);
void FM_copy_mat(double *, int, double *, int, int, int);
void FM_copy_trans(double *, int, double *, int, int, int);
void FM_crossprod(double *, double *, int, int, int, double *, int, int, int);
void FM_GAXPY(double *, double, double *, int, int, int, double *, double, int);
double FM_logAbsDet(double *, int, int);
void FM_mult_mat(double *, double *, int, int, int, double *, int, int, int);
void FM_mult_mat_vec(double *, double *, int, int, int, double *);
void FM_mult_triangular(double *, double *, int, int, double *, int);
void FM_rank1_update(double *, int, int, int, double, double *, double *);
void FM_scale_mat(double *, int, double, double *, int, int, int);
void FM_setzero(double *, int, int, int);
double FM_trace(double *, int, int);
void FM_tcrossprod(double *, double *, int, int, int, double *, int, int, int);

/* Householder reflexion */
void FM_householder(double *, int, double *);

/* operations on triangular matrices */
void FM_cpy_lower(double *, int, int, double *, int);
void FM_cpy_upper(double *, int, int, double *, int);
void FM_cpy_lower2upper(double *, int, int, double *);
void FM_cpy_upper2lower(double *, int, int, double *);
double FM_sum_lower_tri(double *, int, int, int);
double FM_sum_upper_tri(double *, int, int, int);

/* matrix factorizations */
void FM_chol_decomp(double *, int, int, int, int *);
void FM_lu_decomp(double *, int, int, int, int *, int *);
void FM_QR_decomp(double *, int, int, int, double *, int *);
void FM_QL_decomp(double *, int, int, int, double *, int *);
void FM_LQ_decomp(double *, int, int, int, double *, int *);
void FM_svd_decomp(double *, int, int, int, double *, int, double *, double *, int, int, int *);

/* QR, QL and LQ operations */
void FM_QR_qy(double *, int, int, int, double *, double *, int, int, int, int *);
void FM_QR_qty(double *, int, int, int, double *, double *, int, int, int, int *);
void FM_QR_fitted(double *, int, int, int, double *, double *, int, int, int, int, double *);
void FM_QR_store_R(double *, int, int, double *, int);
void FM_QL_qy(double *, int, int, int, double *, double *, int, int, int, int *);
void FM_QL_qty(double *, int, int, int, double *, double *, int, int, int, int *);
void FM_LQ_yq(double *, int, int, int, double *, double *, int, int, int, int *);
void FM_LQ_yqt(double *, int, int, int, double *, double *, int, int, int, int *);
void FM_LQ_store_L(double *, int, int, double *, int);

/* matrix inversion and linear solvers */
void FM_backsolve(double *, int, int, double *, int, int, int *);
void FM_chol_inverse(double *, int, int, int, int *);
void FM_forwardsolve(double *, int, int, double *, int, int, int *);
void FM_invert_mat(double *, int, int, int *);
void FM_invert_triangular(double *, int, int, int, int *);
void FM_lu_inverse(double *, int, int, int *, int *);
void FM_lu_solve(double *, int, int, int *, double *, int, int, int *);

/* least square procedures */
void FM_lsfit(double *, int, int, int, double *, int, int, double *, int *);
void FM_gls_GQR(double *, int, int, int, double *, double *, double *, int *);

/* distances */
double FM_pythag(double, double);
double FM_mahalanobis(double *, int, double *, double *);
void FM_mahal_distances(double *, int, int, double *, double *, int, double *);

/* Wilson-Hilferty transformation */
void FM_WH_chisq(double *, int, int, double *);
void FM_WH_gamma(double *, int, double, double, double *);
void FM_WH_F(double *, int, int, double, double *);
void FM_WH_Laplace(double *, int, int, double *);

/* products */
void FM_two_product_FMA(double, double, double *, double *);
void FM_compensated_product(double *, int, double *);

/* descriptive statistics */
void FM_center_and_Scatter(double *, int, int, double *, double *, double *);
void FM_cov4th(double *, int, int, double *, double *);
void FM_cov_MSSD(double *, int, int, double *, double *);
double FM_find_quantile(double *, int, int);
void FM_geometric_mean(double *, int, double *);
void FM_mean_and_var(double *, int, double *, double *);
void FM_mediancenter(double *, int, int, double *, int *);
void FM_moments(double *, int, double *, double *, double *, double *);
void FM_online_covariance(double *, double *, int, double *, double *, double *, double *, double *);
void FM_online_center(double *, int, int, double *, double *);
void FM_online_mad(double *, int, double *);
void FM_skewness_and_kurtosis(double *, int, int, double *, double *, double *, int);

/* correlation structures */
void FM_cor_AR1(double *, int, double);
void FM_cor_CS(double *, int, double);

/* misc */
void FM_centering(double *, int, int, double *);
void FM_cov2cor(double *, int);
void FM_krylov_mat(double *, int, int, double *, int, double *, int, int *);
void FM_polymat(double *, int, int, double *, int, double *, int, int *);
void FM_sherman_morrison(double *, int, int, double *, double *, int);

/* 'DEBUG' routine */
void FM_print_mat(double *, int, int, int, char *);

#endif /* FASTMATRIX_H */
