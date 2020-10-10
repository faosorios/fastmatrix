/* $ID: ols.c, last updated 2020-10-08, F.Osorio */

#include "fastmatrix.h"

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
  *RSS = norm_sqr(y + p, 1, rdf);
}
