/* $ID: ols.c, last updated 2020-10-08, F.Osorio */

#include "fastmatrix.h"

void
OLS_cg(double *x, int *ldx, int *nrow, int *ncol, double *y, double *coef, double *tol, int *maxiter, int *info)
{ /* ordinary least-squares fit using conjugate gradients */
  int iter = 0, n = *nrow, p = *ncol;
  double delta, gamma, lambda, u, v, *d, *g, *h, *work;

  /* initialization */
  d    = (double *) Calloc(p, double);
  g    = (double *) Calloc(p, double);
  h    = (double *) Calloc(p, double);
  work = (double *) Calloc(n, double);

  /* warming-up */
  FM_crossprod(g, x, *ldx, n, p, y, n, n, 1);
  BLAS1_scale(-1.0, g, 1, p);
  Memcpy(d, g, p);
  gamma = FM_norm_sqr(g, 1, p);

  /* iteration */
  while (gamma > *tol) {
    FM_mult_mat(work, x, *ldx, n, p, d, p, p, 1);
    FM_crossprod(h, x, *ldx, n, p, work, n, n, 1);
    u = BLAS1_dot_product(d, 1, h, 1, p);
    v = FM_norm_sqr(g, 1, p);
    lambda = -v / u;
    BLAS1_axpy(lambda, d, 1, coef, 1, p);
    BLAS1_axpy(lambda, h, 1, g, 1, p);
    delta = FM_norm_sqr(g, 1, p) / v;
    FM_norm_sqr(g, 1, p);
    BLAS1_scale(delta, d, 1, p);
    BLAS1_axpy(1.0, g, 1, d, 1, p);

    gamma = FM_norm_sqr(g, 1, p);
    iter++;
    if (iter > *maxiter)
      break; /* maximum of iterations exceeded */
  }
  *info = iter;

  Free(d); Free(g); Free(h); Free(work);
}

void
OLS_qr(double *x, int *ldx, int *nrow, int *ncol, double *y, double *qraux, double *coef, double *fitted, double *resid, double *RSS)
{ /* ordinary least-squares fit using the QR decomposition */
  int info = 0, job = 0, n = *nrow, p = *ncol, rdf;

  /* QR decomposition of model matrix */
  FM_QR_decomp(x, *ldx, n, p, qraux, &info);
  if (info)
    error("QR_decomp in OLS_qr gave error code %d", info);

  /* copying response */
  Memcpy(resid, y, n);

  /* solve the LS problem */
  FM_QR_qty(x, *ldx, n, p, qraux, y, n, n, 1, &info);
  if (info)
    error("QR_qty in OLS_qr gave error code %d", info);
  Memcpy(coef, y, p);
  FM_backsolve(x, *ldx, p, coef, p, 1, &info);
  if (info)
    error("DTRTRS in OLS_qr gave error code %d", info);

  /* compute fitted values */
  FM_QR_fitted(x, *ldx, n, p, qraux, y, n, n, 1, job, fitted);

  /* compute residuals */
  BLAS1_axpy(-1.0, fitted, 1, resid, 1, n);

  /* residual sum of squares */
  rdf  = n - p;
  *RSS = FM_norm_sqr(y + p, 1, rdf);
}
