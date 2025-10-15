/* $ID: schur.c, last updated 2025-10-13, F.Osorio */

#include "fastmatrix.h"

void
schur_dcmp(double *a, int *lda, int *n, int *task, double *re, double *im, double *v, int *ldv, int *info)
{ /* wrapper to 'FM_schur_decomp' */
  FM_schur_decomp(a, *lda, *n, *task, re, im, v, *ldv, info);
}
