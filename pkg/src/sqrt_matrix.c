/* ID: sqrt_matrix.c, last updated 2025-05-26, F.Osorio */

#include "fastmatrix.h"

void
sqrt_mat_DB(double *a, int *lda, int *n, int *info, int *maxiter, double *tolerance, int *numIter)
{ /* Newton's method for solving x %*% x = a. where a is an n-by-n 
   * matrix using the Denman & Beavers iteration by
   * Denman & Beavers (1976). Appl. Math. Comput. 2, 63-94 */
  double conv, *b, *diff, *old, *x, *y;
  int errcode, iter = 0, job = 1, p = *n;

  /* test the input parameters */
  *info = 0;
  if (p < 0) {
    *info = -3;
  } else if (*lda < MAX(1, p)) {
    *info = -2;
  }
  if (*info != 0) return;

  /* quick return if possible */
  if (p == 0)
    return;

  /* initialization */
  b     = (double *) R_Calloc(p * p, double);
  diff  = (double *) R_Calloc(p * p, double);
  old   = (double *) R_Calloc(p * p, double);
  x     = (double *) R_Calloc(p * p, double);
  y     = (double *) R_Calloc(p * p, double);
  FM_copy_mat(x, p, a, *lda, p, p);
  for (int j = 0; j < p; j++)
    y[j * (p + 1)] = 1.0;

  /* main loop */
  repeat {
    FM_copy_mat(old, p, x, p, p, p);

    /* inversion of 'y' matrix */
    FM_copy_mat(b, p, y, p, p, p);
    errcode = 0;
    FM_invert_mat(b, p, p, &errcode);
    if (errcode != 0)
      error("DGELS in sqrt_mat_DB gave error code %d", errcode);

    /* inversion of 'x_old' matrix */
    FM_copy_mat(a, p, old, p, p, p);
    errcode = 0;
    FM_invert_mat(a, p, p, &errcode);
    if (errcode != 0)
      error("DGELS in sqrt_mat_DB gave error code %d", errcode);

    /* updating 'x' and 'y' matrices */
    for (int j = 0; j < p; j++) {
      for (int i = 0; i < p; i++) {
        x[i + j * p] = 0.5 * (old[i + j * p] + b[i + j * p]);
        y[i + j * p] = 0.5 * (y[i + j * p] + a[i + j * p]);
      }
    }

    iter++;

    /* eval convergence */
    for (int j = 0; j < p; j++) {
      for (int i = 0; i < p; i++) {
        diff[i + j * p] = x[i + j * p] - old[i + j * p];
      }
    }
    matrix_norm(diff, &p, &p, &p, &job, &conv);
    if (conv < *tolerance)
      break; /* successful completion */
    if (iter >= *maxiter)
      break; /* maximum number of iterations exceeded */
  }
  FM_copy_mat(a, *lda, x, p, p, p); /* a <- x */
  *numIter = iter;

  R_Free(b), R_Free(diff), R_Free(old), R_Free(x), R_Free(y); 
}
