/* $ID: wilson_hilferty.c, last updated 2020-11-24, F.Osorio */

#include "fastmatrix.h"

/* Wilson-Hilferty transformation */

void
wilson_hilferty_chisq(double *distances, int *n, int *p, double *z)
{ /* Wilson-Hilferty transformation for chi-squared variables (to be called by R) */
  FM_WH_chisq(distances, *n, *p, z);
}

void
FM_WH_chisq(double *distances, int n, int p, double *z)
{ /* Wilson-Hilferty transformation for chi-squared variables */
  double f, q = 1. / 3., s = 2. / 9.;

  for (int i = 0; i < n; i++) {
    f = *distances++ / p;
    *z++ = (R_pow(f, q) - (1. - s / p)) / sqrt(s / p);
  }
}

void
FM_WH_F(double *distances, int n, int p, double eta, double *z)
{ /* Wilson-Hilferty transformation for F variables */
  double f, q = 1. / 3., r = 2. / 3., s = 2. / 9.;

  for (int i = 0; i < n; i++) {
    f = *distances++ / p;
    f /= 1. - 2. * eta;
    *z++ = ((1. - s * eta) * R_pow(f, q) - (1. - s / p)) / sqrt(s * eta * R_pow(f, r) + s / p);
  }
}
