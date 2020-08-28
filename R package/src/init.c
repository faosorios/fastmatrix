/* $ID: init.c, last updated 2020-08-23, F.Osorio */

#include <R_ext/Rdynload.h>
#include "array.h"
#include "commutation.h"
#include "duplication.h"
#include "matrix.h"
#include "norms.h"
#include "power_method.h"
#include "sweep_operator.h"
#include "sherman_morrison.h"
#include "utils.h"

static const R_CMethodDef CEntries[]  = {
  {"dupl_cols",         (DL_FUNC) &dupl_cols,                    2},
  {"dupl_left_mult",    (DL_FUNC) &dupl_left_mult,               8},
  {"dupl_left_trans",   (DL_FUNC) &dupl_left_trans,              9},
  {"dupl_right_mult",   (DL_FUNC) &dupl_right_mult,              9},
  {"dupl_right_trans",  (DL_FUNC) &dupl_right_trans,             8},
  {"duplication_mat",   (DL_FUNC) &duplication_mat,              4},
  {"equilibrate",       (DL_FUNC) &equilibrate,                  7},
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
  {"hadamard",          (DL_FUNC) &F77_NAME(hadamard),           4},
  {"inner_frobenius",   (DL_FUNC) &F77_NAME(inner_frobenius),    7},
  {NULL, NULL, 0}
};

void R_init_heavy(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
