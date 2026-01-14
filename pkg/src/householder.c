/* $ID: householder.c, last updated 2026-01-12, F.Osorio */

#include "fastmatrix.h"

void
house_vec(double *x, int *n, double *u, double *tau)
{ /* Householder vector */
  int p = *n;

  *tau = FM_house(x, p, u);
}

void
house_mat(double *u, double *tau, int *n, double *mat)
{ /* Householder reflector */
  int p = *n;
  double beta = *tau;

  FM_householder_mat(u, beta, p, mat);
}

void 
house_prod_vec(double *x, int *n, double *u, double *tau)
{ /* applying a Householder reflection to a vector */
  int p = *n;
  double beta = *tau;

  FM_householder_vec(x, p, u, beta);
}

void
house_prod_mat(double *a, int *lda, int *nrow, int *ncol, int *job, double *u, double *tau)
{ /* applying a Householder reflection to a matrix */
  int n = *nrow, p = *ncol, task = *job;
  double beta = *tau;

  FM_householder_prod(a, *lda, n, p, task, u, beta);
}

double 
FM_house(double *x, int n, double *u)
{ /* creates the Householder vector */
  double beta, mu, v, s;

  if (n <= 1) /* quick return */
    return 0.0;

  s = FM_norm_sqr(x + 1, 1, n - 1);
  u[0] = 1.0;
  Memcpy(u + 1, x + 1, n - 1);

  v = x[0];
  if ((s == 0.0) && (v >= 0.0))
    return 0.0;
  else if ((s == 0.0) && (v < 0.0))
    return -2.0;
  
  mu = sqrt(SQR(v) + s);
  if (v <= 0.0)
    v = v - mu;
  else
    v = -s / (v + mu);
  u[0] = v;

  beta = 2.0 * SQR(v) / (s + SQR(v));
  BLAS1_scale(1.0 / v, u, 1, n);
  return beta;
}

void
FM_householder_mat(double *u, double tau, int n, double *mat)
{ /* creates the Householder matrix based on the vector u(1:n) */

  if (n <= 0) /* quick return */
    return;
  
  for (int i = 0; i < n; i++) {
    mat[i * (n + 1)] = 1.0 - tau * SQR(u[i]);
    for (int j = i + 1; j < n; j++) {
      *(mat + i + j * n) = -tau * u[i] * u[j];
      *(mat + j + i * n) = *(mat + i + j * n);
    }
  }
}

void
FM_householder_vec(double *x, int n, double *u, double tau)
{ /* applying a Householder reflection to a vector */
  double prod;

  if (n <= 0) /* quick return */
    return;

  prod = BLAS1_dot_product(u, 1, x, 1, n);
  tau *= -prod;
  BLAS1_axpy(tau, u, 1, x, 1, n);
}

void
FM_householder_prod(double *a, int lda, int n, int p, int job, double *u, double tau)
{ /* applying a Householder reflection to a matrix */
  char *trans = "T", *notrans = "N";
  double *v;

  if (n <= 0) /* quick return */
    return;

  if (job) { /* left: P %*% A */
    v = (double *) R_Calloc(p, double);
    BLAS2_gemv(1.0, a, lda, n, p, trans, u, 1, 0.0, v, 1);
  } else { /* right:  A %*% P */
    v = (double *) R_Calloc(n, double);
    BLAS2_gemv(1.0, a, lda, n, p, notrans, u, 1, 0.0, v, 1);
  }

  if (job) /* left: P %*% A */
    FM_rank1_update(a, lda, n, p, -tau, u, v);
  else /* right:  A %*% P */
    FM_rank1_update(a, lda, n, p, -tau, v, u);

  R_Free(v);
}
