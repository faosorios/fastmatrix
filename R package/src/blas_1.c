/* $ID: blas_1.c, last updated 2020-08-31, F.Osorio */

#include "fastmatrix.h"

/* BLAS level 1 wrappers */

void
BLAS1_axpy(double alpha, double *x, int incx, double *y, int incy, int n)
{ /* y <- alpha * x + y (AXPY operation) */
  F77_CALL(daxpy)(&n, &alpha, x, &incx, y, &incy);
}

void
BLAS1_copy(double *y, int incy, double *x, int incx, int n)
{ /* y <- x (alternative to Memcpy with increments not equal to 1) */
  F77_CALL(dcopy)(&n, x, &incx, y, &incy);
}

double
BLAS1_dot_product(double *x, int incx, double *y, int incy, int n)
{ /* sum(x * y) */
  return F77_CALL(ddot)(&n, x, &incx, y, &incy);
}

int
BLAS1_index_max(double *x, int inc, int n)
{ /* index of element having maximum absolute value */
  return F77_CALL(idamax)(&n, x, &inc);
}

double
BLAS1_norm_two(double *x, int inc, int n)
{ /* sqrt(sum(x * x)) */
  return F77_CALL(dnrm2)(&n, x, &inc);
}

void
BLAS1_scale(double alpha, double *x, int inc, int n)
{ /* x <- alpha * x (x is overwritten) */
  F77_CALL(dscal)(&n, &alpha, x, &inc);
}

double
BLAS1_sum_abs(double *x, int inc, int n)
{ /* sum(abs(x)) */
  return F77_CALL(dasum)(&n, x, &inc);
}

void
BLAS1_swap(double *x, int incx, double *y, int incy, int n)
{ /* x <-> y */
  F77_CALL(dswap)(&n, x, &incx, y, &incy);
}
