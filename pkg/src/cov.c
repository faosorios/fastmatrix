/* ID: cov.c, last updated 2022-06-01, F.Osorio */

#include "fastmatrix.h"

void
cov_weighted(double *x, int *nobs, int *vars, double *weights, double *mean, double *cov)
{ /* wrapper for 'FM_mean_and_Scatter' */
  int n = *nobs, p = *vars;

  FM_center_and_Scatter(x, n, p, weights, mean, cov);
}

void
cov_MSSD(double *x, int *nobs, int *vars, double *mean, double *cov)
{ /* wrapper for 'FM_cov_MSSD' */
  int n = *nobs, p = *vars;

  FM_cov_MSSD(x, n, p, mean, cov);
}

void
cov4th(double *x, int *nobs, int *vars, double *mean, double *cov)
{ /* wrapper for 'FM_cov4th' */
  int n = *nobs, p = *vars;

  FM_cov4th(x, n, p, mean, cov);
}

void
Psi2Q(double *Psi, double *s, int *p)
{ /* scales the 'Psi' matrix related to the log-transformation
   * (see Section 5 of Harris (1985). Biometrika 72, 103-107) */
  for (int i = 0; i < *p; i++) {
    Psi[i * (*p + 1)] /= SQR(s[i]);
    for (int j = i + 1; j < *p; j++) {
      *(Psi + i + j * *p) = *(Psi + i + j * *p) / (s[i] * s[j]);
      *(Psi + j + i * *p) = *(Psi + i + j * *p);
    }
  }
}

void
FM_cov2cor(double *cov, int p)
{ /* scales a 'covariance' matrix into the corresponding correlation matrix */
  double *s;

  s = (double *) Calloc(p, double);

  for (int i = 0; i < p; i++)
    s[i] = cov[i * (p + 1)];

  for (int i = 0; i < p; i++) {
    cov[i * (p + 1)] = 1.0;
    for (int j = i + 1; j < p; j++) {
      *(cov + i + j * p) = *(cov + i + j * p) / sqrt(s[i] * s[j]);
      *(cov + j + i * p) = *(cov + i + j * p);
    }
  }
  Free(s);
}
