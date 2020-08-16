/* $ID: init.c, last updated 2020-08-10, F.Osorio */

#include <R_ext/Rdynload.h>
#include "array.h"
#include "duplication.h"
#include "matrix.h"
#include "norms.h"
#include "utils.h"

static const R_CMethodDef CEntries[]  = {
  {"dupl_cols",         (DL_FUNC) &dupl_cols,                    2},
  {"dupl_left_mult",    (DL_FUNC) &dupl_left_mult,               8},
  {"dupl_left_trans",   (DL_FUNC) &dupl_left_trans,              9},
  {"dupl_right_mult",   (DL_FUNC) &dupl_right_mult,              9},
  {"dupl_right_trans",  (DL_FUNC) &dupl_right_trans,             8},
  {"duplication_mat",   (DL_FUNC) &duplication_mat,              4},
  {"hadamard_prod",     (DL_FUNC) &hadamard_prod,                4},
  {"mat2vech",          (DL_FUNC) &mat2vech,                     4},
  {"matrix_norm",       (DL_FUNC) &matrix_norm,                  6},
  {"norm_one",          (DL_FUNC) &norm_one,                     4},
  {"norm_two",          (DL_FUNC) &norm_two,                     4},
  {"norm_inf",          (DL_FUNC) &norm_inf,                     4},
  {"norm_minkowski",    (DL_FUNC) &norm_minkowski,               5},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
  {"arraymult",         (DL_FUNC) &F77_NAME(arraymult),         14},
  {"bracketprod",       (DL_FUNC) &F77_NAME(bracketprod),       11},
  {"hadamard",          (DL_FUNC) &F77_NAME(hadamard),           4},
  {"inner_frobenius",   (DL_FUNC) &F77_NAME(inner_frobenius),    7},
  {"sweepop",           (DL_FUNC) &F77_NAME(sweepop),            7},
  {"vech",              (DL_FUNC) &F77_NAME(vech),               4},
  {NULL, NULL, 0}
};

void R_init_heavy(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
