/* $ID: power_method.c, last updated 2023-02-25, F.Osorio */

#include "fastmatrix.h"

void
power_method(double *a, int *lda, int *p, double *x, double *lambda, int *maxiter, double *tolerance, int *numIter)
{ /* power method to approximate dominant eigenvalue and eigenvector */
  char *notrans = "N";
  double conv, newLambda, *u = NULL, *v = NULL;
  int iter = 0, n = *p;

  u = (double *) Calloc(n, double);
  v = (double *) Calloc(n, double);

  /* normalizing initial vector */
  Memcpy(u, x, n); /* u <- x */
  FM_normalize(u, 1, n);

  /* main loop */
  repeat {
    BLAS2_gemv(1.0, a, *lda, n, n, notrans, u, 1, 0.0, v, 1);
    FM_normalize(v, 1, n);
    newLambda = OMO_quadf(a, *lda, n, v);

    iter++;

    /* eval convergence */
    conv = fabs(newLambda - *lambda);
    if (conv < *tolerance)
      break; /* successful completion */
    if (iter >= *maxiter)
      break; /* maximum number of iterations exceeded */

    *lambda = newLambda;
    Memcpy(u, v, n); /* u <- v */
  }
  Memcpy(x, v, n); /* x <- v */
  *lambda = newLambda;
  *numIter = iter;

  Free(u); Free(v);
}

void
inverse_power(double *a, int *lda, int *p, double *x, double *lambda, int *maxiter, double *tolerance, int *numIter)
{ /* inverse power method to find the smallest eigenvalue and eigenvector */
  double conv, newLambda, *u = NULL, *v = NULL;
  int iter = 0, n = *p, one = 1, *pivot = NULL;

  u = (double *) Calloc(n, double);
  v = (double *) Calloc(n, double);
  pivot = (int *) Calloc(n, int);

  /* normalizing initial vector */
  Memcpy(u, x, n); /* u <- x */
  FM_normalize(u, 1, n);

  /* perform LU decomposition */
  lu_dcmp(a, lda, &n, p, pivot);

  /* main loop */
  repeat {
    Memcpy(v, u, n); /* v <- u */
    lu_solve(a, lda, &n, pivot, v, &n, &one);
    FM_normalize(v, 1, n);
    newLambda = OMO_quadf(a, *lda, n, v);

    iter++;

    /* eval convergence */
    conv = fabs(newLambda - *lambda);
    if (conv < *tolerance)
      break; /* successful completion */
    if (iter >= *maxiter)
      break; /* maximum number of iterations exceeded */

    *lambda = newLambda;
    Memcpy(u, v, n); /* u <- v */
  }
  Memcpy(x, v, n); /* x <- v */
  *lambda = newLambda;
  *numIter = iter;

  Free(u); Free(v); Free(pivot);
}

