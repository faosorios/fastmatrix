/* $ID: init.c, last updated 2020-09-05, F.Osorio */

#include "fastmatrix.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, nargs)  {#name, (DL_FUNC) &name, nargs}
#define F77DEF(name, nargs)   {#name, (DL_FUNC) &F77_NAME(name), nargs}

static const R_CMethodDef CEntries[]  = {
  CALLDEF(dupl_cols,        2),
  CALLDEF(dupl_left_mult,   8),
  CALLDEF(dupl_left_trans,  9),
  CALLDEF(dupl_right_mult,  9),
  CALLDEF(dupl_right_trans, 8),
  CALLDEF(duplication_mat,  4),
  CALLDEF(lu_decomp,        5),
  CALLDEF(lu_inverse,       4),
  CALLDEF(lu_solve,         7),
  CALLDEF(mat2vech,         4),
  CALLDEF(matrix_norm,      6),
  CALLDEF(norm_one,         4),
  CALLDEF(norm_two,         4),
  CALLDEF(norm_inf,         4),
  CALLDEF(norm_minkowski,   5),
  CALLDEF(power_method,     9),
  CALLDEF(sherman_morrison, 5),
  CALLDEF(sweep_operator,   6),
  CALLDEF(symmetrizer_prod, 6),
  {NULL, NULL, 0}
};

static const R_FortranMethodDef F77Entries[] = {
  F77DEF(arraymult,         14),
  F77DEF(bracketprod,       11),
  F77DEF(comm_rows,          3),
  F77DEF(comm_left_mult,    10),
  F77DEF(comm_right_mult,   10),
  F77DEF(commutation_mat,    6),
  F77DEF(equilibrate_cols,   8),
  F77DEF(hadamard_prod,      4),
  F77DEF(inner_frobenius,    7),
  F77DEF(pivot_mat,          4),
  F77DEF(symmetrizer_mat,    8),
  {NULL, NULL, 0}
};

void R_init_fastmatrix(DllInfo *dll)
{
  /* Register the internal routines. We have no .Call or .External calls */
  R_registerRoutines(dll, CEntries, NULL, F77Entries, NULL);
  R_useDynamicSymbols(dll, FALSE);

  /* BLAS-1 wrappers callable from other packages */
  R_RegisterCCallable("fastmatrix", "BLAS1_axpy",         (DL_FUNC) &BLAS1_axpy);
  R_RegisterCCallable("fastmatrix", "BLAS1_copy",         (DL_FUNC) &BLAS1_copy);
  R_RegisterCCallable("fastmatrix", "BLAS1_dot_product",  (DL_FUNC) &BLAS1_dot_product);
  R_RegisterCCallable("fastmatrix", "BLAS1_index_max",    (DL_FUNC) &BLAS1_index_max);
  R_RegisterCCallable("fastmatrix", "BLAS1_norm_two",     (DL_FUNC) &BLAS1_norm_two);
  R_RegisterCCallable("fastmatrix", "BLAS1_scale",        (DL_FUNC) &BLAS1_scale);
  R_RegisterCCallable("fastmatrix", "BLAS1_sum_abs",      (DL_FUNC) &BLAS1_sum_abs);
  R_RegisterCCallable("fastmatrix", "BLAS1_swap",         (DL_FUNC) &BLAS1_swap);

  /* BLAS-2 wrappers callable from other packages */
  R_RegisterCCallable("fastmatrix", "BLAS2_gemv",   (DL_FUNC) &BLAS2_gemv);
  R_RegisterCCallable("fastmatrix", "BLAS2_symv",   (DL_FUNC) &BLAS2_symv);
  R_RegisterCCallable("fastmatrix", "BLAS2_trmv",   (DL_FUNC) &BLAS2_trmv);
  R_RegisterCCallable("fastmatrix", "BLAS2_trsv",   (DL_FUNC) &BLAS2_trsv);
  R_RegisterCCallable("fastmatrix", "BLAS2_ger",    (DL_FUNC) &BLAS2_ger);
  R_RegisterCCallable("fastmatrix", "BLAS2_syr",    (DL_FUNC) &BLAS2_syr);
  R_RegisterCCallable("fastmatrix", "BLAS2_syr2",   (DL_FUNC) &BLAS2_syr2);

  /* descriptive statistics code callable from other packages */
  R_RegisterCCallable("fastmatrix", "FM_mean_and_var",       (DL_FUNC) &FM_mean_and_var);
  R_RegisterCCallable("fastmatrix", "FM_online_covariance",  (DL_FUNC) &FM_online_covariance);
  R_RegisterCCallable("fastmatrix", "FM_center_and_Scatter", (DL_FUNC) &FM_center_and_Scatter);
  R_RegisterCCallable("fastmatrix", "FM_MSSD",               (DL_FUNC) &FM_MSSD);
  R_RegisterCCallable("fastmatrix", "FM_find_quantile",      (DL_FUNC) &FM_find_quantile);
}
