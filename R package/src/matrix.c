/* ID: matrix.c, last updated 2020-09-03, F.Osorio */

#include "fastmatrix.h"

/* basic matrix manipulations */

void
add_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol)
{ /* y <- y + alpha * x */
  for (int j = 0; j < ncol; j++) {
    BLAS1_axpy(alpha, x, 1, y, 1, nrow);
    y += ldy; x += ldx;
  }
}

void
copy_mat(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
{ /* y <- x[,] */
  for (int j = 0; j < ncol; j++) {
    Memcpy(y, x, nrow);
    y += ldy; x += ldx;
  }
}

void
copy_trans(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
{ /* y <- t(x) */
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++)
      *(y + j + i * ldy) = *(x + i + j * ldx);
  }
}

void
scale_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol)
{ /* y <- alpha * x[,] */
  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++)
      y[i] = alpha * x[i];
    y += ldy; x += ldx;
  }
}

void
setzero(double *y, int ldy, int nrow, int ncol)
{ /* y[,] <- 0, sets all elements of y to 0 */
  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++)
      y[i] = 0.0;
    y += ldy;
  }
}
