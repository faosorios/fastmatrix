/* $ID: R_init_fastmatrix.c, last updated 2025-05-25, F.Osorio */

#include "fastmatrix.h"
#include <R_ext/Rdynload.h>

/* borrowed from Matrix */
#define CALLDEF(name, nargs)  {#name, (DL_FUNC) &name, nargs}
#define EXTDEF(name, nargs)   {#name, (DL_FUNC) &name, nargs}
#define F77DEF(name, nargs)   {#name, (DL_FUNC) &F77_NAME(name), nargs}

static const R_CMethodDef CEntries[]  = {
  CALLDEF(bezier_smoother,        7),
  CALLDEF(central_moments,        6),
  CALLDEF(cg_solver,              9),
  CALLDEF(chol_dcmp,              5),
  CALLDEF(chol_update,            4),
  CALLDEF(cor_AR1,                3),
  CALLDEF(cor_CS,                 3),
  CALLDEF(cov_MSSD,               5),
  CALLDEF(cov_weighted,           6),
  CALLDEF(cov4th,                 5),
  CALLDEF(doornik_hansen,         5),
  CALLDEF(dupl_cols,              2),
  CALLDEF(dupl_left_mult,         8),
  CALLDEF(dupl_left_trans,        9),
  CALLDEF(dupl_right_mult,        9),
  CALLDEF(dupl_right_trans,       8),
  CALLDEF(duplication_mat,        4),
  CALLDEF(equilibrate_pd,         7),
  CALLDEF(equilibrate_sym,        7),
  CALLDEF(equilibrating_sym,      3),
  CALLDEF(geometric_mean,         3),
  CALLDEF(hadamard_prod,          4),
  CALLDEF(jacobi_solver,          9),
  CALLDEF(jarque_bera,            5),
  CALLDEF(kronecker_prod,         7),
  CALLDEF(krylov_mat,             8),
  CALLDEF(lu_dcmp,                6),
  CALLDEF(lu_inverse,             5),
  CALLDEF(lu_solve,               8),
  CALLDEF(mahal_distances,        7),
  CALLDEF(mat2vech,               4),
  CALLDEF(matrix_norm,            6),
  CALLDEF(matrix_polynomial,      8),
  CALLDEF(norm_one,               4),
  CALLDEF(norm_two,               4),
  CALLDEF(norm_inf,               4),
  CALLDEF(norm_minkowski,         5),
  CALLDEF(OLS_cg,                 9),
  CALLDEF(OLS_qr,                10),
  CALLDEF(OLS_ridge,             21),
  CALLDEF(power_method,           8),
  CALLDEF(Psi2Q,                  3),
  CALLDEF(robust_JB,              6),
  CALLDEF(rng_ball,               3),
  CALLDEF(rng_mnorm,              5),
  CALLDEF(rng_sphere,             3),
  CALLDEF(seidel_solver,          9),
  CALLDEF(sherman_morrison,       6),
  CALLDEF(sqrt_mat_DB,            7),
  CALLDEF(skewness_and_kurtosis,  7),
  CALLDEF(svd_dcmp,              11),
  CALLDEF(sweep_operator,         6),
  CALLDEF(symmetrizer_prod,       6),
  CALLDEF(urzua_ALM,              5),
  CALLDEF(wilson_hilferty_chisq,  4),
  CALLDEF(wilson_hilferty_gamma,  5),
  CALLDEF(whitening_chol,         4),
  {NULL, NULL, 0}
};

static const R_FortranMethodDef F77Entries[] = {
  F77DEF(arraymult,              14),
  F77DEF(blinf,                   6),
  F77DEF(bracketprod,            11),
  F77DEF(circulant_mat,           5),
  F77DEF(comm_rows,               3),
  F77DEF(comm_left_mult,         10),
  F77DEF(comm_right_mult,        10),
  F77DEF(commutation_mat,         6),
  F77DEF(equilibrate_cols,        8),
  F77DEF(frank_mat,               4),
  F77DEF(hankel_mat,              6),
  F77DEF(helmert_mat,             4),
  F77DEF(inner_frobenius,         7),
  F77DEF(ldl_dcmp,                5),
  F77DEF(mchol_dcmp,              6),
  F77DEF(median_center,           7),
  F77DEF(murrv,                   7),
  F77DEF(pivot_mat,               4),
  F77DEF(quadf,                   4),
  F77DEF(rhoc_ustat,              6),
  F77DEF(symmetrizer_mat,         8),
  {NULL, NULL, 0}
};

void R_init_fastmatrix(DllInfo *info) {
  /* Register internal routines. We have no .Call or .External calls */
  R_registerRoutines(info,
                     CEntries,    /* .C         R_CMethodDef */
                     NULL,        /* .Call      R_CallMethodDef */
                     F77Entries,  /* .Fortran   R_FortranMethodDef */
                     NULL);       /* .External  R_ExternalMethodDef */
  R_useDynamicSymbols(info, FALSE);

  /* Register functions callable from other packages' C code */
  #define FM_REGDEF(name) R_RegisterCCallable("fastmatrix", #name, (DL_FUNC) name)

  /* BLAS-1 wrappers */
  FM_REGDEF(BLAS1_axpy);
  FM_REGDEF(BLAS1_copy);
  FM_REGDEF(BLAS1_dot_product);
  FM_REGDEF(BLAS1_index_max);
  FM_REGDEF(BLAS1_norm_two);
  FM_REGDEF(BLAS1_rot);
  FM_REGDEF(BLAS1_rotg);
  FM_REGDEF(BLAS1_scale);
  FM_REGDEF(BLAS1_sum_abs);
  FM_REGDEF(BLAS1_swap);
  /* BLAS-2 wrappers */
  FM_REGDEF(BLAS2_gemv);
  FM_REGDEF(BLAS2_symv);
  FM_REGDEF(BLAS2_trmv);
  FM_REGDEF(BLAS2_trsv);
  FM_REGDEF(BLAS2_ger);
  FM_REGDEF(BLAS2_syr);
  FM_REGDEF(BLAS2_syr2);
  /* BLAS-3 wrappers */
  FM_REGDEF(BLAS3_gemm);
  FM_REGDEF(BLAS3_symm);
  FM_REGDEF(BLAS3_syrk);
  FM_REGDEF(BLAS3_trmm);
  FM_REGDEF(BLAS3_trsm);
  /* OMO wrappers */
  FM_REGDEF(OMO_blinf);
  FM_REGDEF(OMO_quadf);
  FM_REGDEF(OMO_murrv);
  /* operations on vectors */
  FM_REGDEF(FM_norm_sqr);
  FM_REGDEF(FM_normalize);
  FM_REGDEF(FM_vecsum);
  /* basic matrix manipulations */
  FM_REGDEF(FM_add_mat);
  FM_REGDEF(FM_copy_mat);
  FM_REGDEF(FM_copy_trans);
  FM_REGDEF(FM_crossprod);
  FM_REGDEF(FM_GAXPY);
  FM_REGDEF(FM_logAbsDet);
  FM_REGDEF(FM_mult_mat);
  FM_REGDEF(FM_mult_mat_vec);
  FM_REGDEF(FM_mult_triangular);
  FM_REGDEF(FM_rank1_update);
  FM_REGDEF(FM_scale_mat);
  FM_REGDEF(FM_setzero);
  FM_REGDEF(FM_trace);
  FM_REGDEF(FM_tcrossprod);
  /* operations on triangular matrices */
  FM_REGDEF(FM_cpy_lower);
  FM_REGDEF(FM_cpy_upper);
  FM_REGDEF(FM_cpy_lower2upper);
  FM_REGDEF(FM_cpy_upper2lower);
  FM_REGDEF(FM_sum_lower_tri);
  FM_REGDEF(FM_sum_upper_tri);
  /* matrix factorizations */
  FM_REGDEF(FM_chol_decomp);
  FM_REGDEF(FM_lu_decomp);
  FM_REGDEF(FM_QR_decomp);
  FM_REGDEF(FM_QL_decomp);
  FM_REGDEF(FM_LQ_decomp);
  FM_REGDEF(FM_svd_decomp);
  /* QR, QL and LQ operations */
  FM_REGDEF(FM_QR_qy);
  FM_REGDEF(FM_QR_qty);
  FM_REGDEF(FM_QR_fitted);
  FM_REGDEF(FM_QR_store_R);
  FM_REGDEF(FM_QL_qy);
  FM_REGDEF(FM_QL_qty);
  FM_REGDEF(FM_LQ_yq);
  FM_REGDEF(FM_LQ_yqt);
  FM_REGDEF(FM_LQ_store_L);
  /* matrix inversion and linear solvers */
  FM_REGDEF(FM_backsolve);
  FM_REGDEF(FM_forwardsolve);
  FM_REGDEF(FM_chol_inverse);
  FM_REGDEF(FM_invert_mat);
  FM_REGDEF(FM_invert_triangular);
  FM_REGDEF(FM_lu_inverse);
  FM_REGDEF(FM_lu_solve);
  /* least square procedures */
  FM_REGDEF(FM_gls_GQR);
  FM_REGDEF(FM_lsfit);
  /* distances */
  FM_REGDEF(FM_pythag);
  FM_REGDEF(FM_mahalanobis);
  /* Wilson-Hilferty transformation */
  FM_REGDEF(FM_WH_chisq);
  FM_REGDEF(FM_WH_gamma);
  FM_REGDEF(FM_WH_F);
  FM_REGDEF(FM_WH_Laplace);
  /* products */
  FM_REGDEF(FM_compensated_product);
  FM_REGDEF(FM_two_product_FMA);
  /* descriptive statistics */
  FM_REGDEF(FM_center_and_Scatter);
  FM_REGDEF(FM_cov_MSSD);
  FM_REGDEF(FM_cov4th);
  FM_REGDEF(FM_find_quantile);
  FM_REGDEF(FM_geometric_mean);
  FM_REGDEF(FM_mean_and_var);
  FM_REGDEF(FM_mediancenter);
  FM_REGDEF(FM_moments);
  FM_REGDEF(FM_online_center);
  FM_REGDEF(FM_online_covariance);
  FM_REGDEF(FM_skewness_and_kurtosis);
  /* misc */
  FM_REGDEF(FM_centering);
  FM_REGDEF(FM_cov2cor);
  FM_REGDEF(FM_krylov_mat);
  FM_REGDEF(FM_sherman_morrison);
  /* 'DEBUG' routine */
  FM_REGDEF(FM_print_mat);

  #undef FM_REGDEF
}
