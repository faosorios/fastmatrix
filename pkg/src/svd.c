/* ID: svd.c, last updated 2020-10-02, F.Osorio */

#include "fastmatrix.h"

/* SVD decomposition */

void
svd_dcmp(double *a, int *lda, int *n, int *p, double *u, int *ldu, double *d, double *vt, int *ldvt, int *job, int *info)
{ /* wrapper to 'FM_svd_decomp' */
  FM_svd_decomp(a, *lda, *n, *p, u, *ldu, d, vt, *ldvt, *job, info);
}
