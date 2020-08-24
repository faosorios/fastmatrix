/* $ID: power_method.c, last updated 2020-08-21, F.Osorio */

#include "base.h"
#include "power_method.h"
#include "vector.h"

void
power_method(double *a, int *lda, int *nrow, int *ncol, double *x, double *lambda,
  int *maxiter, double *tolerance, int *numIter)
{ /* power method to approximate dominant eigenvalue and eigenvector */
  char *notrans = "N";
  double conv, div, newLambda, one = 1.0, zero = 0.0, *z = NULL;
  int inc = 1, iter = 0, n = *nrow, p = *ncol;

  z = (double *) Calloc(n, double);

  /* normalizing initial vector */
  normalize_vec(x, inc, n);

  /* initial estimate */
  F77_CALL(dgemv)(notrans, &n, &p, &one, a, lda, x, &inc, &zero, z, &inc);
  normalize_vec(z, inc, n);
  *lambda = F77_CALL(ddot)(&n, x, &inc, z, &inc);
  Memcpy(x, z, p); /* x <- z */

  /* main loop */
  repeat {
    F77_CALL(dgemv)(notrans, &n, &p, &one, a, lda, x, &inc, &zero, z, &inc);
    normalize_vec(z, inc, n);
    F77_CALL(dgemv)(notrans, &n, &p, &one, a, lda, z, &inc, &zero, x, &inc);
    newLambda = F77_CALL(ddot)(&n, x, &inc, z, &inc);

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
  F77_CALL(dscal)(&n, &div, x, &inc);

  *numIter = iter;

  Free(z);
}
