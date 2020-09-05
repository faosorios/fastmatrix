/* $ID: init.c, last updated 2020-09-05, F.Osorio */

#include "fastmatrix.h"
#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[]  = {
  {"dupl_cols",         (DL_FUNC) &dupl_cols,                    2},
  {"dupl_left_mult",    (DL_FUNC) &dupl_left_mult,               8},
  {"dupl_left_trans",   (DL_FUNC) &dupl_left_trans,              9},
  {"dupl_right_mult",   (DL_FUNC) &dupl_right_mult,              9},
  {"dupl_right_trans",  (DL_FUNC) &dupl_right_trans,             8},
  {"duplication_mat",   (DL_FUNC) &duplication_mat,              4},
  {"mat2vech",          (DL_FUNC) &mat2vech,                     4},
  {"matrix_norm",       (DL_FUNC) &matrix_norm,                  6},
  {"norm_one",          (DL_FUNC) &norm_one,                     4},
  {"norm_two",          (DL_FUNC) &norm_two,                     4},
  {"norm_inf",          (DL_FUNC) &norm_inf,                     4},
  {"norm_minkowski",    (DL_FUNC) &norm_minkowski,               5},
  {"power_method",      (DL_FUNC) &power_method,                 9},
  {"sherman_morrison",  (DL_FUNC) &sherman_morrison,             5},
  {"sweep_operator",    (DL_FUNC) &sweep_operator,               6},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
  {"arraymult",         (DL_FUNC) &F77_NAME(arraymult),         14},
  {"bracketprod",       (DL_FUNC) &F77_NAME(bracketprod),       11},
  {"comm_rows",         (DL_FUNC) &F77_NAME(comm_rows),          3},
  {"comm_left_mult",    (DL_FUNC) &F77_NAME(comm_left_mult),    10},
  {"comm_right_mult",   (DL_FUNC) &F77_NAME(comm_right_mult),   10},
  {"commutation_mat",   (DL_FUNC) &F77_NAME(commutation_mat),    6},
  {"equilibrate_cols",  (DL_FUNC) &F77_NAME(equilibrate_cols),   8},
  {"hadamard_prod",     (DL_FUNC) &F77_NAME(hadamard_prod),      4},
  {"inner_frobenius",   (DL_FUNC) &F77_NAME(inner_frobenius),    7},
  {NULL, NULL, 0}
};

void R_init_fastmatrix(DllInfo *dll)
{
  /* Register the internal routines. We have no .Call or .External calls */
  R_registerRoutines(dll, CEntries, NULL, FortEntries, NULL);
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

  /* C code callable from other packages */
  R_RegisterCCallable("fastmatrix", "center_and_Scatter", (DL_FUNC) &center_and_Scatter);
  R_RegisterCCallable("fastmatrix", "MSSD",               (DL_FUNC) &MSSD);
}
