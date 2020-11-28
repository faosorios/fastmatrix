/* $ID: power_method.c, last updated 2020-09-04, F.Osorio */

#include "fastmatrix.h"

void
power_method(double *a, int *lda, int *nrow, int *ncol, double *x, double *lambda,
  int *maxiter, double *tolerance, int *numIter)
{ /* power method to approximate dominant eigenvalue and eigenvector */
  char *notrans = "N";
  double conv, div, newLambda, *z = NULL;
  int iter = 0, n = *nrow, p = *ncol;

  z = (double *) Calloc(n, double);

  /* normalizing initial vector */
  FM_normalize(x, 1, n);

  /* initial estimate */
  BLAS2_gemv(1.0, a, *lda, n, p, notrans, x, 1, 0.0, z, 1);
  FM_normalize(z, 1, n);
  *lambda = BLAS1_dot_product(x, 1, z, 1, n);
  Memcpy(x, z, p); /* x <- z */

  /* main loop */
  repeat {
    BLAS2_gemv(1.0, a, *lda, n, p, notrans, x, 1, 0.0, z, 1);
    FM_normalize(z, 1, n);
    BLAS2_gemv(1.0, a, *lda, n, p, notrans, z, 1, 0.0, x, 1);
    newLambda = BLAS1_dot_product(x, 1, z, 1, n);

    iter++;

    /* eval convergence */
    conv = fabs(newLambda - *lambda);
    if (conv < *tolerance)
      break; /* successful completion */
    if (iter >= *maxiter)
      break; /* maximum number of iterations exceeded */

    *lambda = newLambda;
    Memcpy(x, z, p); /* x <- z */
  }
  *lambda = newLambda;
  div = 1.0 / *lambda;
  BLAS1_scale(div, x, 1, n);

  *numIter = iter;

  Free(z);
}
