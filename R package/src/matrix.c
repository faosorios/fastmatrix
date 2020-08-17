/* $ID: matrix.c, last updated 2020-08-08, F.Osorio */

#include "base.h"
#include "matrix.h"

void
mat2vech(double *x, int *ldx, int *n, double *y)
{ /* y <- vech(x) */
  int p = *n, k = 0;

  for (int j = 0; j < p; j++) {
    for (int i = j; i < p; i++) {
      y[k] = x[i + j * *ldx];
      k++;
    }
  }
}

void
hadamard_prod(double *x, double *y, int *n, double *prod)
{ /* prod <- x * y */
  for (int i = 0; i < *n; i++)
    *prod++ = *x++ * *y++;
}

void
power_method(double *a, int *lda, int *nrow, int *ncol, double *x, double *lambda,
  int *maxiter, double *tolerance, int *numIter)
{ /* power method to approximate dominant eigenvalue and eigenvector */
  char *notrans = "N";
  double one = 1.0, zero = 0.0, conv, div, norm, newLambda, *z = NULL;
  int inc = 1, iter = 0, n = *nrow, p = *ncol;

  z = (double *) Calloc(n, double);

  /* normalizing initial vector */
  norm = F77_CALL(dnrm2)(&n, x, &inc);
  div  = 1.0 / norm;
  F77_CALL(dscal)(&n, &div, x, &inc);

  /* initial estimate */
  F77_CALL(dgemv)(notrans, &n, &p, &one, a, lda, x, &inc, &zero, z, &inc);
  norm = F77_CALL(dnrm2)(&n, z, &inc);
  div  = 1.0 / norm;
  F77_CALL(dscal)(&n, &div, z, &inc);
  *lambda = F77_CALL(ddot)(&n, x, &inc, z, &inc);
  Memcpy(x, z, p); /* x <- z */

  /* main loop */
  repeat {
    F77_CALL(dgemv)(notrans, &n, &p, &one, a, lda, x, &inc, &zero, z, &inc);
    norm = F77_CALL(dnrm2)(&n, z, &inc);
    div  = 1.0 / norm;
    F77_CALL(dscal)(&n, &div, z, &inc);
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
