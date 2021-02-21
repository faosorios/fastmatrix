/* $ID: ols.c, last updated 2020-10-08, F.Osorio */

#include "fastmatrix.h"

void
OLS_cg(double *x, int *ldx, int *nrow, int *ncol, double *y, double *coef, double *tol, int *maxiter, int *info)
{ /* ordinary least-squares fit using conjugate gradients */
  int iter = 0, n = *nrow, p = *ncol;
  double delta, gamma, lambda, u, v, *d, *g, *h, *work;
  double az, z, accum, scale, ssq;

  /* initialization */
  d    = (double *) Calloc(p, double);
  g    = (double *) Calloc(p, double);
  h    = (double *) Calloc(p, double);
  work = (double *) Calloc(n, double);

  /* warming-up */
  FM_crossprod(g, x, *ldx, n, p, y, n, n, 1);
  /* computing the gradient, search direction and convergence criterion (gamma) */
  scale = 0.0; ssq = 1.0;
  for (int i = 0; i < p; i++) {
    z = -g[i];
    d[i] = g[i] = z;

    if (z != 0.0) {
      az = fabs(z);

      if (scale < az) {
        ssq = 1.0 + ssq * (scale / az) * (scale / az);
        scale = az;
      } else
        ssq += (az / scale) * (az / scale);
    }
  }
  gamma = SQR(scale * sqrt(ssq));

  /* iteration */
  while (gamma > *tol) {
    FM_mult_mat(work, x, *ldx, n, p, d, p, p, 1);
    FM_crossprod(h, x, *ldx, n, p, work, n, n, 1);
    /* computing dot product and step-length (lambda) */
    accum = 0.0; scale = 0.0; ssq = 1.0;
    for (int i = 0; i < p; i++) { /* code re-use! */
      accum += d[i] * h[i];
      z = g[i];

      if (z != 0.0) {
        az = fabs(z);

        if (scale < az) {
          ssq = 1.0 + ssq * (scale / az) * (scale / az);
          scale = az;
        } else
          ssq += (az / scale) * (az / scale);
      }
    }
    u = accum;
    v = SQR(scale * sqrt(ssq));
    lambda = -v / u;
    /* updating coefficients and gradient */
    for (int i = 0; i < p; i++) {
      coef[i] += lambda * d[i];
      g[i] += lambda * h[i];
    }
    delta = FM_norm_sqr(g, 1, p) / v;
    /* updating search direction and convergence criterion (gamma) */
    scale = 0.0; ssq = 1.0;
    for (int i = 0; i < p; i++) {
      z = g[i];
      d[i] = z + delta * d[i];

      if (z != 0.0) {
        az = fabs(z);

        if (scale < az) {
          ssq = 1.0 + ssq * (scale / az) * (scale / az);
          scale = az;
        } else
          ssq += (az / scale) * (az / scale);
      }
    }
    gamma = SQR(scale * sqrt(ssq));

    iter++;
    if (iter > *maxiter)
      break; /* maximum number of iterations exceeded */
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
