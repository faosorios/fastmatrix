/* $ID: blas_1.c, last updated 2021-02-14, F.Osorio */

#include "fastmatrix.h"

/* BLAS level 1 wrappers */

void
BLAS1_axpy(double alpha, double *x, int incx, double *y, int incy, int n)
{ /* y <- alpha * x + y (AXPY operation) */

  /* quick return if possible */
  if (n <= 0 || incx <= 0 || incy <= 0)
    return;
  if (alpha == 0.0)
    return;

  if (incx == 1 && incy == 1) {
    /* code for increments equal to 1 */
    const int m = n % 8;

    for (int i = 0; i < m; i++)
      y[i] += alpha * x[i];

    for (int i = m; i + 7 < n; i += 8) {
      y[i] += alpha * x[i];
      y[i + 1] += alpha * x[i + 1];
      y[i + 2] += alpha * x[i + 2];
      y[i + 3] += alpha * x[i + 3];
      y[i + 4] += alpha * x[i + 4];
      y[i + 5] += alpha * x[i + 5];
      y[i + 6] += alpha * x[i + 6];
      y[i + 7] += alpha * x[i + 7];
    }
  } else {
    /* code for increments not equal to 1 */
    int ix = OFFSET(n, incx);
    int iy = OFFSET(n, incy);

    for (int i = 0; i < n; i++) {
      y[iy] += alpha * x[ix];
      ix += incx;
      iy += incy;
    }
  }
}

void
BLAS1_copy(double *y, int incy, double *x, int incx, int n)
{ /* y <- x (alternative to Memcpy with increments not equal to 1) */

  /* quick return if possible */
  if (n <= 0 || incx <= 0 || incy <= 0)
    return;

  if (incx == 1 && incy == 1) {
    /* code for increments equal to 1 */
    const int m = n % 8;

    for (int i = 0; i < m; i++)
      y[i] = x[i];

    for (int i = m; i + 7 < n; i += 8) {
      y[i] = x[i];
      y[i + 1] = x[i + 1];
      y[i + 2] = x[i + 2];
      y[i + 3] = x[i + 3];
      y[i + 4] = x[i + 4];
      y[i + 5] = x[i + 5];
      y[i + 6] = x[i + 6];
      y[i + 7] = x[i + 7];
    }
  } else {
    /* code for increments not equal to 1 */
    int ix = OFFSET(n, incx);
    int iy = OFFSET(n, incy);

    for (int i = 0; i < n; i++) {
      y[iy] = x[ix];
      ix += incx;
      iy += incy;
    }
  }
}

double
BLAS1_dot_product(double *x, int incx, double *y, int incy, int n)
{ /* sum(x * y) */
  double accum = 0.0;

  /* quick return if possible */
  if (n <= 0 || incx <= 0 || incy <= 0)
    return 0.0;

  if (incx == 1 && incy == 1) {
    /* code for increments equal to 1 */
    const int m = n % 8;

    for (int i = 0; i < m; i++)
      accum += x[i] * y[i];

    for (int i = m; i + 7 < n; i += 8) {
      accum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
            + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4] + x[i + 5] * y[i + 5]
            + x[i + 6] * y[i + 6] + x[i + 7] * y[i + 7];
    }
  } else {
    /* code for increments not equal to 1 */
    int ix = OFFSET(n, incx);
    int iy = OFFSET(n, incy);

    for (int i = 0; i < n; i++) {
      accum += y[iy] * x[ix];
      ix += incx;
      iy += incy;
    }
  }

  return accum;
}

int
BLAS1_index_max(double *x, int inc, int n)
{ /* index of element having maximum absolute value */
  return F77_CALL(idamax)(&n, x, &inc);
}

double
BLAS1_norm_two(double *x, int inc, int n)
{ /* sqrt(sum(x * x)) */
  int ix = 0;
  double az, z, scale = 0.0, ssq = 1.0;

  /* quick return if possible */
  if (n <= 0 || inc <= 0)
    return 0.0;
  if (n == 1)
    return fabs(x[0]);

  for (int i = 0; i < n; i++) {
    /* this loop is equivalent to the call to the LAPACK auxiliary
     * routine: DLASSQ which updates a sum of squares represented
     * in scaled form */
    z = x[ix];

    if (z != 0.0) {
      az = fabs(z);

      if (scale < az) {
        ssq = 1.0 + ssq * (scale / az) * (scale / az);
        scale = az;
      } else
        ssq += (az / scale) * (az / scale);
    }

    ix += inc;
  }

  return scale * sqrt(ssq);
}

void
BLAS1_scale(double alpha, double *x, int inc, int n)
{ /* x <- alpha * x (x is overwritten) */

  /* quick return if possible */
  if (n <= 0 || inc <= 0)
    return;

  if (inc == 1) {
    /* code for increment equal to 1 */
    const int m = n % 8;

    for (int i = 0; i < m; i++)
      x[i] *= alpha;

    for (int i = m; i + 7 < n; i += 8) {
      x[i] *= alpha;
      x[i + 1] *= alpha;
      x[i + 2] *= alpha;
      x[i + 3] *= alpha;
      x[i + 4] *= alpha;
      x[i + 5] *= alpha;
      x[i + 6] *= alpha;
      x[i + 7] *= alpha;
    }
  } else {
    /* code for increment not equal to 1 */
    int ix = OFFSET(n, inc);

    for (int i = 0; i < n; i++) {
      x[ix] *= alpha;
      ix += inc;
    }
  }
}

double
BLAS1_sum_abs(double *x, int inc, int n)
{ /* sum(abs(x)) */
  double accum = 0.0;

  /* quick return if possible */
  if (n <= 0 || inc <= 0)
    return 0.0;
  if (n == 1)
    return fabs(x[0]);

  if (inc == 1) {
    /* code for increment equal to 1 */
    const int m = n % 8;

    for (int i = 0; i < m; i++)
      accum += fabs(x[i]);

    for (int i = m; i + 7 < n; i += 8) {
      accum += fabs(x[i]) + fabs(x[i + 1]) + fabs(x[i + 2]) + fabs(x[i + 3])
            + fabs(x[i + 4]) + fabs(x[i + 5]) + fabs(x[i + 6]) + fabs(x[i + 7]);
    }
  } else {
    /* code for increment not equal to 1 */
    int ix = OFFSET(n, inc);

    for (int i = 0; i < n; i++) {
      accum += fabs(x[ix]);
      ix += inc;
    }
  }

  return accum;
}

void
BLAS1_swap(double *x, int incx, double *y, int incy, int n)
{ /* interchanges two vectors: x <-> y */
  double aux;

  /* quick return if possible */
  if (n <= 0 || incx <= 0 || incy <= 0)
    return;

  if (incx == 1 && incy == 1) {
    /* code for increments equal to 1 */
    const int m = n % 8;

    for (int i = 0; i < m; i++) {
      aux = x[i]; x[i] = y[i]; y[i] = aux;
    }

    for (int i = m; i + 7 < n; i += 8) {
      aux = x[i]; x[i] = y[i]; y[i] = aux;
      aux = x[i + 1]; x[i + 1] = y[i + 1]; y[i + 1] = aux;
      aux = x[i + 2]; x[i + 2] = y[i + 2]; y[i + 2] = aux;
      aux = x[i + 3]; x[i + 3] = y[i + 3]; y[i + 3] = aux;
      aux = x[i + 4]; x[i + 4] = y[i + 4]; y[i + 4] = aux;
      aux = x[i + 5]; x[i + 5] = y[i + 5]; y[i + 5] = aux;
      aux = x[i + 6]; x[i + 6] = y[i + 6]; y[i + 6] = aux;
      aux = x[i + 7]; x[i + 7] = y[i + 7]; y[i + 7] = aux;
    }
  } else {
    /* code for increments not equal to 1 */
    int ix = OFFSET(n, incx);
    int iy = OFFSET(n, incy);

    for (int i = 0; i < n; i++) {
      aux = x[ix];
      x[ix] = y[iy];
      y[iy] = aux;
      ix += incx;
      iy += incy;
    }
  }
}
