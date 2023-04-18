/* ID: cor_struct.c, last updated 2022-07-06, F.Osorio */

#include "fastmatrix.h"

void
cor_AR1(double *cor, int *p, double *rho)
{ /* wrapper to 'FM_cor_AR1'*/
  FM_cor_AR1(cor, *p, *rho);
}

void
FM_cor_AR1(double *cor, int p, double rho)
{ /* AR(1) correlation matrix */
  int pow;

  /* fast return if possible */
  if (rho == 0.0) {
    for (int i = 0; i < p; i++)
      cor[i * (p + 1)] = 1.0;
    return;
  }

  /* autoregressive correlation structure */
  for (int i = 0; i < p; i++) {
    cor[i * (p + 1)] = 1.0;
    for (int j = i + 1; j < p; j++) {
      pow = (int) abs(i - j);
      *(cor + i + j * p) = R_pow_di(rho, pow);
      *(cor + j + i * p) = *(cor + i + j * p);
    }
  }
}

void
cor_CS(double *cor, int *p, double *rho)
{ /* wrapper to 'FM_cor_CS'*/
  FM_cor_CS(cor, *p, *rho);
}

void
FM_cor_CS(double *cor, int p, double rho)
{ /* compound symmetry correlation matrix */

  /* fast return if possible */
  if (rho == 0.0) {
    for (int i = 0; i < p; i++)
      cor[i * (p + 1)] = 1.0;
    return;
  }

  /* compound symmetry correlation structure */
  for (int i = 0; i < p; i++) {
    cor[i * (p + 1)] = 1.0;
    for (int j = i + 1; j < p; j++) {
      *(cor + i + j * p) = rho;
      *(cor + j + i * p) = *(cor + i + j * p);
    }
  }
}
