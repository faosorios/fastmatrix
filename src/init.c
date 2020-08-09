#include <R_ext/Rdynload.h>
#include "duplication.h"
#include "matrix.h"
#include "specmat.h"
#include "utils.h"

static const R_CMethodDef CEntries[]  = {
  {"dupl_cols",         (DL_FUNC) &dupl_cols,           2},
  {"dupl_left_mult",    (DL_FUNC) &dupl_left_mult,      8},
  {"dupl_left_trans",   (DL_FUNC) &dupl_left_trans,     9},
  {"dupl_right_mult",   (DL_FUNC) &dupl_right_mult,     9},
  {"dupl_right_trans",  (DL_FUNC) &dupl_right_trans,    8},
  {"duplication_mat",   (DL_FUNC) &duplication_mat,     4},
  {"hadamard_prod",     (DL_FUNC) &hadamard_prod,       4},
  {"mat2vech",          (DL_FUNC) &mat2vech,            4},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
  {"hadamard",          (DL_FUNC) &F77_NAME(hadamard),  4},
  {"vech",              (DL_FUNC) &F77_NAME(vech),      4},
  {NULL, NULL, 0}
};

void R_init_heavy(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
