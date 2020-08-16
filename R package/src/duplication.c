/* $ID: duplication.c, last updated 2020-08-10, F.Osorio */

#include "base.h"
#include "duplication.h"

void
dupl_cols(int *order, int *cols)
{ /* compact information to build a duplication matrix */
  int n = *order, k = 0;

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      k++;
      cols[i + j * n] = cols[j + i * n] = k;
    }
  }
}

void
duplication_mat(int *x, int *ldx, int *n, int *col)
{ /* creates the duplication matrix of order 'n' */
  int nrow = SQR(*n), pos;

  for (int i = 0; i < nrow; i++) {
    pos = col[i] - 1; /* index correction */
    x[i + pos * *ldx] = 1;
  }
}

void
dupl_left_mult(double *a, int *lda, int *arow, int *acol, int *col, int *n, double *b, int *ldb)
{ /* computes: B <- Dn %*% A */
  int nrow = SQR(*n), ncol = *n * (*n + 1) / 2, pos;

  if (*arow != ncol)
    return;

  for (int j = 0; j < *acol; j++) {
    for (int i = 0; i < nrow; i++) {
      pos = col[i] - 1; /* index correction */
      b[i + j * *ldb] = a[pos + j * *lda];
    }
  }
}

void
dupl_left_trans(double *a, int *lda, int *arow, int *acol, int *col, int *n, int *counts, double *b, int *ldb)
{ /* computes: B <- t(Dn) %*% A */
  int nrow = *n * (*n + 1) / 2, ncol = SQR(*n), pos1 = 0, pos2 = 0, k;

  if (*arow != ncol)
    return;

  for (int j = 0; j < *acol; j++) {
    k = 0;
    for (int i = 0; i < nrow; i++) {
      if (counts[i] < 2) {
        pos1 = col[k] - 1; /* index correction */
        b[i + j * *ldb] = a[pos1 + j * *lda];
        k++;
      } else {
        pos1 = col[k] - 1;     /* index correction */
        pos2 = col[k + 1] - 1; /* index correction */
        b[i + j * *ldb] = a[pos1 + j * *lda] + a[pos2 + j * *lda];
        k += 2;
      }
    }
  }
}

void
dupl_right_mult(double *a, int *lda, int *arow, int *acol, int *col, int *n, int *counts, double *b, int *ldb)
{ /* computes: B <- A %*% Dn */
  int nrow = SQR(*n), ncol = *n * (*n + 1) / 2, pos1 = 0, pos2 = 0, k;

  if (*acol != nrow)
    return;

  k = 0;
  for (int j = 0; j < ncol; j++) {
    if (counts[j] < 2) {
      pos1 = col[k] - 1; /* index correction */
      k++;
    } else {
      pos1 = col[k] - 1;     /* index correction */
      pos2 = col[k + 1] - 1; /* index correction */
      k += 2;
    }
    for (int i = 0; i < *arow; i++) {
      if (counts[j] < 2)
        b[i + j * *ldb] = a[i + pos1 * *lda];
      else
        b[i + j * *ldb] = a[i + pos1 * *lda] + a[i + pos2 * *lda];
    }
  }
}

void
dupl_right_trans(double *a, int *lda, int *arow, int *acol, int *col, int *n, double *b, int *ldb)
{ /* computes: B <- A %*% t(Dn) */
  int nrow = *n * (*n + 1) / 2, ncol = SQR(*n), pos;

  if (*acol != nrow)
    return;

  for (int j = 0; j < ncol; j++) {
    pos = col[j] - 1; /* index correction */
    for (int i = 0; i < *arow; i++) {
      b[i + j * *ldb] = a[i + pos * *lda];
    }
  }
}
